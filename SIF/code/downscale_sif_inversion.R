# Downscale OCO-2 SIF_757nm from satellite footprints to the sub_road_grid
# cells around Lamont, OK, with a regularized (Bayesian / Tikhonov) linear
# inversion, one inversion per day with target-mode observations.
#
# Forward model:  y = K x + e
#   y  SIF_757nm of each footprint                         (n_obs)
#   x  SIF_757nm of each grid cell (the unknowns)          (n_cells)
#   K  Jacobian: K[i, j] = area(footprint i in cell j) / area(footprint i in grid),
#      so every row sums to 1 and is the fraction of footprint i in each cell
#   e  retrieval error, Se = diag(SIF_Uncertainty_757nm^2 + model_err^2)
#
# Regularization: Gaussian prior x ~ N(xa, Sa) with xa = the day's mean SIF in
# every cell and Sa[j, k] = prior_sd^2 * exp(-d_jk / corr_length_km), d = distance
# between cell centroids. Solution (Rodgers 2000, optimal estimation):
#   S_hat = (K' Se^-1 K + Sa^-1)^-1
#   x_hat = xa + S_hat K' Se^-1 (y - K xa)
#   A     = S_hat K' Se^-1 K          (averaging kernel; DOFS = trace(A))
#
# Effective resolution: row j of A is how the estimate in cell j responds to the
# true field, so its width is cell j's resolution. Per cell, with the row as a
# response density A[j, k] / area_k over cells within res_radius_km:
#   res_area_pr_km2   participation area (sum A)^2 / sum(A^2 / area), positive
#                     part of the row; equals R for a flat response over area R
#                     and the cell's own area for a perfect response. The
#                     effective resolution used in the paper.
#   res_area_km2      area where the density is at least half its peak (FWHM
#                     area); on this grid it rarely resolves below one cell
#   ak_row_sum        sum of row j: near 1 where the estimate comes from data
#                     (possibly from neighbouring cells), near 0 where it is prior
# Both areas include the smoothing of the prior correlation (corr_length_km).
# They are meaningful only where ak_row_sum is well above 0.
#
# Input:  oco2_sif_raster_extent.gpkg (from subset_oco2_to_extent.R),
#         sub_road_grid.shp
# Output (in out_dir):
#   downscaled_sif757.gpkg   layer "cells": one row per cell per day with
#                            sif757, sif757_sd, avg_kernel, n_footprints,
#                            res_area_km2, res_area_pr_km2, ak_row_sum
#   A/<date>.rds             full averaging kernel for the days in save_A_days
#   daily_summary.csv        one row per day: n_obs, prior, DOFS, fit stats
#   maps/<date>.png          footprints vs downscaled SIF vs uncertainty

library(sf)
library(data.table)
library(Matrix)

# ---- Settings ---------------------------------------------------------------
proj_dir  <- path.expand("~/Documents/cropland_downscaling_paper")
oco_file  <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
out_dir   <- file.path(proj_dir, "downscaled_sif757")

target_mode    <- 2       # Metadata_MeasurementMode for target observations
quality_flags  <- c(0, 1) # Quality_Flag values to keep
min_coverage   <- 0.5     # drop footprints with less than half their area in the grid
min_obs        <- 50      # skip days with fewer usable footprints

corr_length_km <- 1.0     # prior correlation length between cells: sqrt of the area-weighted mean
                          # area of contiguous same-crop patches, CDL 2014-2023 (1.03 km2; crop_patch_sizes.R)
prior_sd       <- NULL    # prior SD of cell SIF; NULL = estimate per day (below)
min_prior_sd   <- 0.05    # floor for the per-day estimate
model_err      <- 0       # extra error added in quadrature to each footprint (W/m2/sr/um)
make_maps      <- TRUE
res_radius_km  <- 6       # search radius for the width of each averaging-kernel row
save_A_days    <- as.Date("2019-04-19")   # days whose full averaging kernel is saved (Figure 2)

dir.create(file.path(out_dir, "maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "A"), showWarnings = FALSE)

# ---- Grid -------------------------------------------------------------------
grid <- st_read(grid_file, quiet = TRUE)
# GeoPackage reserves "FID" for its own feature ID, so keep the original as grid_fid
names(grid)[toupper(names(grid)) == "FID"] <- "grid_fid"
grid$cell_id <- seq_len(nrow(grid))
grid_crs <- st_crs(grid)                      # NAD83 / Oklahoma North (m)
n_cells <- nrow(grid)
domain <- st_union(st_geometry(grid))

# Prior correlation between cells (fixed for all days)
cent <- st_coordinates(st_centroid(st_geometry(grid)))
dist_km <- as.matrix(dist(cent)) / 1000
corr <- exp(-dist_km / corr_length_km)
corr_inv <- solve(corr)
cell_area <- as.numeric(st_area(grid)) / 1e6  # km2
near_mask <- dist_km <= res_radius_km

# ---- Width of each averaging-kernel row (effective resolution) -----------------
resolution_metrics <- function(A) {
  dens <- sweep(A, 2, cell_area, "/") * near_mask        # response per km2, local only
  peak <- apply(dens, 1, max)
  half <- dens >= 0.5 * peak & near_mask
  Ap <- pmax(A, 0) * near_mask
  list(res_area_km2 = as.vector(half %*% cell_area),
       res_area_pr_km2 = rowSums(Ap)^2 / rowSums(sweep(Ap^2, 2, cell_area, "/")),
       ak_row_sum = rowSums(A))
}

# ---- OCO-2 target observations ----------------------------------------------
oco <- st_read(oco_file, quiet = TRUE,
               query = paste("SELECT * FROM oco2_sif WHERE Metadata_MeasurementMode =",
                             target_mode))
# Normalise SIF (and its uncertainty) to a nadir view: vza_correction.R
source(file.path(proj_dir, "vza_correction.R"))
oco <- correct_vza(oco)
oco <- oco[oco$Quality_Flag %in% quality_flags &
           !is.na(oco$Science_SIF_757nm) &
           !is.na(oco$Science_SIF_Uncertainty_757nm), ]
oco$day <- as.Date(oco$time_utc)

# Footprints are convex quadrilaterals; some have their corners stored out of
# order (self-intersecting), which the convex hull of the 4 corners fixes
oco <- st_transform(oco, grid_crs)
st_geometry(oco) <- st_convex_hull(st_geometry(oco))

days <- sort(unique(oco$day))
cat(length(days), "days with target observations\n")

# ---- Jacobian for one day ---------------------------------------------------
build_jacobian <- function(fp) {
  fp$obs_id <- seq_len(nrow(fp))
  pieces <- suppressWarnings(st_intersection(fp[, "obs_id"], grid[, "cell_id"]))
  pieces$area <- as.numeric(st_area(pieces))
  pieces <- pieces[pieces$area > 0, ]

  in_grid <- tapply(pieces$area, pieces$obs_id, sum)
  coverage <- numeric(nrow(fp))
  coverage[as.integer(names(in_grid))] <- in_grid / as.numeric(st_area(fp))[as.integer(names(in_grid))]
  keep <- which(coverage >= min_coverage)

  pieces <- pieces[pieces$obs_id %in% keep, ]
  row <- match(pieces$obs_id, keep)
  K <- sparseMatrix(i = row, j = pieces$cell_id, x = pieces$area,
                    dims = c(length(keep), n_cells))
  K <- Diagonal(x = 1 / rowSums(K)) %*% K       # each row sums to 1
  list(K = K, keep = keep, coverage = coverage[keep])
}

# ---- Inversion for one day --------------------------------------------------
invert_day <- function(d) {
  fp <- oco[oco$day == d, ]
  jac <- build_jacobian(fp)
  if (length(jac$keep) < min_obs) return(NULL)
  K <- jac$K
  y <- fp$Science_SIF_757nm[jac$keep]
  se2 <- fp$Science_SIF_Uncertainty_757nm[jac$keep]^2 + model_err^2

  # Prior: the day's error-weighted mean everywhere; spread = the part of the
  # footprint-to-footprint variance not explained by retrieval noise
  xa_val <- sum(y / se2) / sum(1 / se2)
  sd_a <- if (is.null(prior_sd)) max(sqrt(max(var(y) - mean(se2), 0)), min_prior_sd) else prior_sd
  xa <- rep(xa_val, n_cells)
  Sa_inv <- corr_inv / sd_a^2

  KtSe <- t(K) %*% Diagonal(x = 1 / se2)
  H <- as.matrix(KtSe %*% K)                    # K' Se^-1 K   (n_cells x n_cells)
  S_hat <- solve(H + Sa_inv)
  x_hat <- as.vector(xa + S_hat %*% as.vector(KtSe %*% (y - as.vector(K %*% xa))))
  A <- S_hat %*% H                              # full averaging kernel
  A_diag <- diag(A)
  rm_ <- resolution_metrics(A)
  if (d %in% save_A_days) saveRDS(A, file.path(out_dir, "A", paste0(d, ".rds")))

  resid <- y - as.vector(K %*% x_hat)
  cells <- data.table(cell_id = grid$cell_id, date = d,
                      sif757 = x_hat, sif757_sd = sqrt(diag(S_hat)),
                      avg_kernel = A_diag,
                      n_footprints = as.vector(colSums(K > 0)),
                      res_area_km2 = rm_$res_area_km2, res_area_pr_km2 = rm_$res_area_pr_km2,
                      ak_row_sum = rm_$ak_row_sum)
  ok <- A_diag > 0.5
  summary <- data.table(date = d, n_obs = length(y),
                        n_dropped_coverage = nrow(fp) - length(y),
                        obs_mean = mean(y), obs_sd = sd(y),
                        median_obs_unc = median(sqrt(se2)),
                        prior_mean = xa_val, prior_sd = sd_a,
                        dofs = sum(A_diag),
                        chi2_reduced = sum(resid^2 / se2) / length(y),
                        rmse_fit = sqrt(mean(resid^2)),
                        rmse_prior = sqrt(mean((y - xa_val)^2)),
                        cells_constrained = sum(ok),
                        res_area_med_constrained = if (any(ok)) median(rm_$res_area_pr_km2[ok]) else NA_real_,
                        area_finer_than_footprint = sum(cell_area[rm_$res_area_pr_km2 < 2.9 & rm_$ak_row_sum > 0.5]))
  list(cells = cells, summary = summary, fp = fp[jac$keep, ], y = y)
}

# ---- Maps -------------------------------------------------------------------
map_day <- function(res) {
  d <- res$summary$date
  g <- merge(grid[, "cell_id"], res$cells, by = "cell_id")
  lim <- quantile(c(res$y, g$sif757), c(0.02, 0.98), na.rm = TRUE)
  pal <- hcl.colors(50, "YlGn", rev = TRUE)
  col_of <- function(v) pal[cut(pmin(pmax(v, lim[1]), lim[2]),
                                seq(lim[1], lim[2], length.out = 51), include.lowest = TRUE)]
  png(file.path(out_dir, "maps", paste0(d, ".png")), width = 2700, height = 1000, res = 200)
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
  o <- order(res$y)
  plot(domain, border = "grey40", main = sprintf("%s: %d footprints (SIF757)", d, length(res$y)))
  plot(st_geometry(res$fp)[o], col = adjustcolor(col_of(res$y[o]), 0.35), border = NA, add = TRUE)
  plot(st_geometry(g), col = col_of(g$sif757), border = NA,
       main = sprintf("Downscaled SIF757 (DOFS %.0f)", res$summary$dofs))
  sd_pal <- hcl.colors(50, "Purples", rev = TRUE)
  plot(st_geometry(g), col = sd_pal[cut(g$sif757_sd, 50)], border = NA,
       main = sprintf("Posterior SD (%.2f-%.2f)", min(g$sif757_sd), max(g$sif757_sd)))
  dev.off()
}

# ---- Run --------------------------------------------------------------------
all_cells <- list(); all_summary <- list()
for (k in seq_along(days)) {
  res <- tryCatch(invert_day(days[k]), error = function(e) {
    message(days[k], ": failed - ", conditionMessage(e)); NULL
  })
  if (is.null(res)) next
  all_cells[[k]] <- res$cells
  all_summary[[k]] <- res$summary
  if (make_maps) map_day(res)
  s <- res$summary
  message(sprintf("%s  n=%4d  prior %.2f +/- %.2f  DOFS %5.1f  chi2 %.2f",
                  s$date, s$n_obs, s$prior_mean, s$prior_sd, s$dofs, s$chi2_reduced))
}

cells <- rbindlist(all_cells)
summary <- rbindlist(all_summary)
fwrite(summary, file.path(out_dir, "daily_summary.csv"))

out <- merge(grid[, c("cell_id", setdiff(names(grid), c("cell_id", attr(grid, "sf_column"))))],
             cells, by = "cell_id")
out <- out[order(out$date, out$cell_id), ]
st_write(out, file.path(out_dir, "downscaled_sif757.gpkg"), layer = "cells",
         delete_dsn = TRUE, quiet = TRUE)
cat("\n", nrow(summary), "days downscaled;", nrow(out), "cell-days written to",
    file.path(out_dir, "downscaled_sif757.gpkg"), "\n")

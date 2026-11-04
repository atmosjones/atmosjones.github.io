# Estimate OCO-2 SIF_757nm for each crop class around Lamont, OK, on each day
# with target observations: a regularized inversion like
# downscale_sif_inversion.R, but with one unknown per crop class instead of
# one per grid cell, so the data, not the prior, determine the answer.
#
# Forward model:  y = Kc x + e
#   y   SIF_757nm of each footprint                              (n_obs)
#   x   SIF_757nm of each crop class that day                    (n_classes)
#   Kc  = K M: K[i, j] = fraction of footprint i (within the grid) in cell j
#       (rows sum to 1); M[j, c] = 1 if cell j's majority crop is in class c.
#       So Kc[i, c] = fraction of footprint i in class c; rows still sum to 1
#   e   retrieval error, Se = diag(SIF_Uncertainty_757nm^2 + model_err^2)
# Prior: x ~ N(day mean, prior_sd^2 I), weak, to keep sparsely observed
# classes stable. Posterior as in Rodgers (2000):
#   S_hat = (Kc' Se^-1 Kc + Sa^-1)^-1,  x_hat = xa + S_hat Kc' Se^-1 (y - Kc xa)
#
# Each cell counts entirely as its majority crop (sub_road_crop maj_frac is
# often well under 1), so class values are for "cells where this crop
# dominates", not pure fields.
#
# Input:  oco2_sif_raster_extent.gpkg, sub_road_grid.shp, and the crop columns
#         in downscaled_sif757/downscaled_sif757.gpkg (add_crop_to_sif.R; that
#         file's crop_year already handles fall-planted winter wheat)
# Output (in out_dir):
#   sif757_by_crop.csv         one row per day and crop class
#   sif757_by_crop_days.csv    one row per day: n_obs, DOFS, fit statistics
#   sif757_by_crop_posterior.rds   posterior covariance per day (classes are
#                                  spatially mixed, so their errors correlate)

library(sf)
library(data.table)
library(Matrix)

# ---- Settings ---------------------------------------------------------------
proj_dir  <- path.expand("~/Documents/cropland_downscaling_paper")
oco_file  <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
crop_file <- file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg")
out_dir   <- file.path(proj_dir, "downscaled_sif757")

target_mode   <- 2
quality_flags <- c(0, 1)
min_coverage  <- 0.5    # drop footprints with less than half their area in the grid
min_obs       <- 50
prior_sd      <- 1.0    # weak prior on each class (W/m2/sr/um)
model_err     <- 0      # extra footprint error added in quadrature

# Crop classes to estimate; CDL names not listed here go to "Other"
crop_groups <- list(
  "Winter Wheat"             = "Winter Wheat",
  "Dbl Crop WinWht/Soybeans" = "Dbl Crop WinWht/Soybeans",
  "Soybeans"                 = "Soybeans",
  "Corn"                     = "Corn",
  "Sorghum"                  = "Sorghum",
  "Cotton"                   = "Cotton",
  "Grassland/Pasture"        = "Grassland/Pasture",
  "Fallow/Idle Cropland"     = "Fallow/Idle Cropland",
  "Developed"                = c("Developed/Open Space", "Developed/Low Intensity",
                                 "Developed/Med Intensity", "Developed/High Intensity"),
  "Trees/Wetlands"           = c("Deciduous Forest", "Woody Wetlands",
                                 "Herbaceous Wetlands")
)
class_levels <- c(names(crop_groups), "Other")
group_of <- function(crop_name) {
  g <- rep("Other", length(crop_name))
  for (k in names(crop_groups)) g[crop_name %in% crop_groups[[k]]] <- k
  factor(g, levels = class_levels)
}

# ---- Grid and crop classes per day ------------------------------------------
grid <- st_read(grid_file, quiet = TRUE)
grid$cell_id <- seq_len(nrow(grid))
grid_crs <- st_crs(grid)
n_cells <- nrow(grid)
cell_area <- as.numeric(st_area(grid))

cell_crops <- st_read(crop_file, query = "SELECT cell_id, date, crop_name FROM cells",
                      quiet = TRUE)
cell_crops <- as.data.table(if (inherits(cell_crops, "sf")) st_drop_geometry(cell_crops) else cell_crops)
if (!"crop_name" %in% names(cell_crops)) {
  stop("Crop columns missing from ", crop_file, "; run add_crop_to_sif.R first")
}
cell_crops[, date := as.Date(date)]
# Days without a crop map (years with no CDL yet) can't be split by crop
no_crop <- cell_crops[, .(missing = anyNA(crop_name)), by = date][missing == TRUE, date]
if (length(no_crop)) message("Skipping ", length(no_crop), " days without crop labels")
cell_crops <- cell_crops[!date %in% no_crop]
cell_crops[, crop_class := group_of(crop_name)]

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
oco <- st_transform(oco, grid_crs)
st_geometry(oco) <- st_convex_hull(st_geometry(oco))  # fix out-of-order corners
days <- sort(intersect(unique(oco$day), unique(cell_crops$date)))
cat(length(days), "days with target observations and crop maps\n")

# ---- Footprint-to-cell Jacobian (as in downscale_sif_inversion.R) -----------
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
  K <- sparseMatrix(i = match(pieces$obs_id, keep), j = pieces$cell_id,
                    x = pieces$area, dims = c(length(keep), n_cells))
  K <- Diagonal(x = 1 / rowSums(K)) %*% K       # each row sums to 1
  list(K = K, keep = keep)
}

# ---- Inversion for one day --------------------------------------------------
invert_day <- function(d) {
  fp <- oco[oco$day == d, ]
  jac <- build_jacobian(fp)
  if (length(jac$keep) < min_obs) return(NULL)
  y <- fp$Science_SIF_757nm[jac$keep]
  se2 <- fp$Science_SIF_Uncertainty_757nm[jac$keep]^2 + model_err^2

  # Cell -> class aggregation for this day's crop map
  cls <- cell_crops[date == d][order(cell_id), crop_class]
  M <- sparseMatrix(i = seq_len(n_cells), j = as.integer(cls),
                    dims = c(n_cells, length(class_levels)))
  Kc <- as.matrix(jac$K %*% M)                    # rows sum to 1
  colnames(Kc) <- class_levels

  # Only classes some footprint actually sees
  seen <- colSums(Kc) > 0
  Kc <- Kc[, seen, drop = FALSE]
  p <- ncol(Kc)

  xa_val <- sum(y / se2) / sum(1 / se2)
  xa <- rep(xa_val, p)
  KtSe <- t(Kc / se2)                             # Kc' Se^-1
  H <- KtSe %*% Kc
  S_hat <- solve(H + diag(1 / prior_sd^2, p))
  x_hat <- as.vector(xa + S_hat %*% (KtSe %*% (y - Kc %*% xa)))
  A_diag <- diag(S_hat %*% H)
  resid <- y - as.vector(Kc %*% x_hat)
  chi2 <- sum(resid^2 / se2) / (length(y) - p)
  dimnames(S_hat) <- list(colnames(Kc), colnames(Kc))

  class_area <- tapply(cell_area, cls, sum)[colnames(Kc)] / 1e6
  class_cells <- table(cls)[colnames(Kc)]
  est <- data.table(date = d, crop_class = colnames(Kc),
                    sif757 = x_hat,
                    sif757_sd = sqrt(diag(S_hat)),
                    # SD scaled up by the misfit, for when chi2 > 1
                    sif757_sd_scaled = sqrt(diag(S_hat) * max(chi2, 1)),
                    avg_kernel = A_diag,
                    footprint_weight = colSums(Kc),  # effective no. of footprints
                    n_cells = as.integer(class_cells),
                    area_km2 = as.numeric(class_area))
  days_row <- data.table(date = d, n_obs = length(y), n_classes = p,
                         prior_mean = xa_val, dofs = sum(A_diag),
                         chi2_reduced = chi2,
                         rmse_fit = sqrt(mean(resid^2)),
                         rmse_prior = sqrt(mean((y - xa_val)^2)))
  list(est = est, day = days_row, S_hat = S_hat)
}

# ---- Run --------------------------------------------------------------------
est <- list(); day_rows <- list(); posterior <- list()
for (k in seq_along(days)) {
  res <- tryCatch(invert_day(days[k]), error = function(e) {
    message(days[k], ": failed - ", conditionMessage(e)); NULL
  })
  if (is.null(res)) next
  est[[k]] <- res$est
  day_rows[[k]] <- res$day
  posterior[[as.character(days[k])]] <- res$S_hat
  message(sprintf("%s  n=%4d  classes %2d  DOFS %4.1f  chi2 %.2f",
                  days[k], res$day$n_obs, res$day$n_classes, res$day$dofs,
                  res$day$chi2_reduced))
}
est <- rbindlist(est)
day_rows <- rbindlist(day_rows)

fwrite(est, file.path(out_dir, "sif757_by_crop.csv"))
fwrite(day_rows, file.path(out_dir, "sif757_by_crop_days.csv"))
saveRDS(posterior, file.path(out_dir, "sif757_by_crop_posterior.rds"))
cat("\n", nrow(day_rows), "days;", nrow(est), "class-days written to", out_dir, "\n")

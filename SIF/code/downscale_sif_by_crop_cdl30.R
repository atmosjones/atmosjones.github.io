# Estimate OCO-2 SIF_757nm for each crop class around Lamont, OK, on each day
# with target observations, using the 30 m USDA Cropland Data Layer directly:
# the same inversion as downscale_sif_by_crop.R, but each footprint's row of
# the Jacobian is the exact fraction of its area in each crop class, from the
# 30 m pixels, instead of assigning whole grid cells to their majority crop.
#
# Forward model:  y = Kc x + e
#   y   SIF_757nm of each footprint                               (n_obs)
#   x   SIF_757nm of each crop class that day                     (n_classes)
#   Kc  Kc[i, c] = fraction of footprint i covered by class c (30 m CDL);
#       rows sum to 1 over the whole footprint
#   e   retrieval error, Se = diag(SIF_Uncertainty_757nm^2 + model_err^2)
# Prior x ~ N(day mean, prior_sd^2 I), weak. Posterior as in Rodgers (2000).
#
# Crop classes and the crop-year rule (fall-planted winter wheat) are in
# crop_classes.R.
#
# Usage:  Rscript downscale_sif_by_crop_cdl30.R [757|740]   (default 757)
#   757 = Science/SIF_757nm; 740 = SIF_740nm (Lite product, from 757 + 771),
#   the wavelength to compare with TROPOMI-based products such as SIFNet.
#
# Input:  oco2_sif_raster_extent.gpkg, cdl_30m/CDL_<year>.tif (CropScape)
# Output (in out_dir; <band> = 757 or 740):
#   sif<band>_by_crop_cdl30.csv           one row per day and crop class
#   sif<band>_by_crop_cdl30_days.csv      one row per day: n_obs, DOFS, fit stats
#   sif<band>_by_crop_cdl30_posterior.rds posterior covariance per day
#   sif<band>_by_crop_cdl30_design.rds    per day: sounding IDs, y, error
#                                         variance and the Jacobian Kc

library(terra)
library(sf)
library(data.table)
library(exactextractr)

# ---- Settings ---------------------------------------------------------------
proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
oco_file <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
cdl_dir  <- file.path(proj_dir, "cdl_30m")
out_dir  <- file.path(proj_dir, "downscaled_sif757")

args <- commandArgs(trailingOnly = TRUE)
band <- if (length(args) >= 1) args[1] else "757"
sif_var <- c("757" = "Science_SIF_757nm", "740" = "SIF_740nm")[[band]]
unc_var <- c("757" = "Science_SIF_Uncertainty_757nm", "740" = "SIF_Uncertainty_740nm")[[band]]
tag <- paste0("sif", band)

target_mode   <- 2
quality_flags <- c(0, 1)
min_cdl_cover <- 0.95   # drop footprints less than 95% covered by CDL pixels
min_obs       <- 50
prior_sd      <- 1.0    # weak prior on each class (W/m2/sr/um)
model_err     <- 0      # extra footprint error added in quadrature

# Crop classes, CDL rasters and class_raster(date), incl. the fall winter
# wheat rule
source(file.path(proj_dir, "crop_classes.R"))

# ---- OCO-2 target observations ----------------------------------------------
oco <- st_read(oco_file, quiet = TRUE,
               query = paste("SELECT * FROM oco2_sif WHERE Metadata_MeasurementMode =",
                             target_mode))
# Normalise SIF (and its uncertainty) to a nadir view: vza_correction.R
source(file.path(proj_dir, "vza_correction.R"))
oco <- correct_vza(oco)
oco <- oco[oco$Quality_Flag %in% quality_flags &
           !is.na(oco[[sif_var]]) & !is.na(oco[[unc_var]]), ]
oco$day <- as.Date(oco$time_utc)
oco <- st_transform(oco, cdl_crs)
st_geometry(oco) <- st_convex_hull(st_geometry(oco))  # fix out-of-order corners
# A few footprints have collapsed corners (a line or point, no area); drop them
degenerate <- !st_geometry_type(oco) %in% "POLYGON" | as.numeric(st_area(oco)) <= 0
if (any(degenerate)) message("Dropping ", sum(degenerate), " zero-area footprints")
oco <- oco[!degenerate, ]
days <- sort(unique(oco$day))
no_cdl <- !format(days, "%Y") %in% names(cdl)
if (any(no_cdl)) message("Skipping ", sum(no_cdl), " days in years without a CDL yet (",
                         paste(unique(format(days[no_cdl], "%Y")), collapse = ", "), ")")
days <- days[!no_cdl]
cat(length(days), "days with target observations and a CDL\n")

# ---- Inversion for one day --------------------------------------------------
invert_day <- function(d) {
  fp <- oco[oco$day == d, ]
  cls <- class_raster(d)

  # Fraction of each footprint in each class, from the 30 m pixels
  # (coverage_fraction = share of each pixel inside the footprint)
  px <- exact_extract(cls, fp, progress = FALSE)
  fr <- t(vapply(px, function(p) {
    ok <- !is.na(p$value)
    w <- tapply(p$coverage_fraction[ok],
                factor(p$value[ok], levels = seq_along(class_levels)), sum)
    w[is.na(w)] <- 0
    c(w, cover = sum(p$coverage_fraction[ok]) / sum(p$coverage_fraction))
  }, numeric(length(class_levels) + 1)))
  keep <- which(fr[, length(class_levels) + 1] >= min_cdl_cover)
  if (length(keep) < min_obs) return(NULL)
  Kc <- fr[keep, seq_along(class_levels), drop = FALSE]
  Kc <- Kc / rowSums(Kc)                          # rows sum to 1
  colnames(Kc) <- class_levels

  y <- fp[[sif_var]][keep]
  se2 <- fp[[unc_var]][keep]^2 + model_err^2

  seen <- colSums(Kc) > 0
  Kc <- Kc[, seen, drop = FALSE]
  p <- ncol(Kc)

  xa_val <- sum(y / se2) / sum(1 / se2)
  xa <- rep(xa_val, p)
  KtSe <- t(Kc / se2)
  H <- KtSe %*% Kc
  S_hat <- solve(H + diag(1 / prior_sd^2, p))
  x_hat <- as.vector(xa + S_hat %*% (KtSe %*% (y - Kc %*% xa)))
  A_diag <- diag(S_hat %*% H)
  resid <- y - as.vector(Kc %*% x_hat)
  chi2 <- sum(resid^2 / se2) / (length(y) - p)
  dimnames(S_hat) <- list(colnames(Kc), colnames(Kc))

  # Class area within the area the footprints cover that day
  covered <- vect(st_union(st_geometry(fp[keep, ])))
  area_in <- freq(mask(crop(cls, covered), covered))
  area_km2 <- setNames(area_in$count * prod(res(cls)) / 1e6,
                       class_levels[area_in$value])[colnames(Kc)]

  est <- data.table(date = d, crop_class = colnames(Kc),
                    sif757 = x_hat,
                    sif757_sd = sqrt(diag(S_hat)),
                    sif757_sd_scaled = sqrt(diag(S_hat) * max(chi2, 1)),
                    avg_kernel = A_diag,
                    footprint_weight = colSums(Kc),
                    mean_frac_in_footprints = colMeans(Kc),
                    area_km2 = as.numeric(area_km2))
  days_row <- data.table(date = d,
                         lite_build = if ("lite_build" %in% names(fp)) paste(unique(fp$lite_build), collapse = "+") else NA_character_,
                         n_obs = length(y),
                         n_dropped = nrow(fp) - length(y), n_classes = p,
                         prior_mean = xa_val, dofs = sum(A_diag),
                         chi2_reduced = chi2,
                         rmse_fit = sqrt(mean(resid^2)),
                         rmse_prior = sqrt(mean((y - xa_val)^2)))
  setnames(est, c("sif757", "sif757_sd", "sif757_sd_scaled"),
           c(tag, paste0(tag, "_sd"), paste0(tag, "_sd_scaled")))
  design <- list(sounding_id = fp$Metadata_SoundingId[keep], y = y, se2 = se2,
                 Kc = Kc)
  list(est = est, day = days_row, S_hat = S_hat, design = design)
}

# ---- Run --------------------------------------------------------------------
est <- list(); day_rows <- list(); posterior <- list(); designs <- list()
for (k in seq_along(days)) {
  res <- tryCatch(invert_day(days[k]), error = function(e) {
    message(days[k], ": failed - ", conditionMessage(e)); NULL
  })
  if (is.null(res)) next
  est[[k]] <- res$est
  day_rows[[k]] <- res$day
  posterior[[as.character(days[k])]] <- res$S_hat
  designs[[as.character(days[k])]] <- res$design
  message(sprintf("%s  n=%4d  classes %2d  DOFS %4.1f  chi2 %.2f",
                  days[k], res$day$n_obs, res$day$n_classes, res$day$dofs,
                  res$day$chi2_reduced))
}
est <- rbindlist(est)
day_rows <- rbindlist(day_rows)

fwrite(est, file.path(out_dir, paste0(tag, "_by_crop_cdl30.csv")))
fwrite(day_rows, file.path(out_dir, paste0(tag, "_by_crop_cdl30_days.csv")))
saveRDS(posterior, file.path(out_dir, paste0(tag, "_by_crop_cdl30_posterior.rds")))
saveRDS(designs, file.path(out_dir, paste0(tag, "_by_crop_cdl30_design.rds")))
cat("\n", nrow(day_rows), "days;", nrow(est), "class-days written to", out_dir, "\n")

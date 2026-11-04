# Compare the OCO-2 SIF downscaling with HLS NDVI (30 m, Landsat + Sentinel-2).
#
# 1. Footprints: how well NDVI explains the OCO-2 footprint SIF on each target
#    day, against the crop-class model and the two combined. All predictions
#    are cross-validated (fit on a random half of the day's footprints,
#    predict the other half). NDVI comes from the clearest HLS scene within
#    max_gap_days of the OCO-2 day.
# 2. Crop classes: crop-class SIF (downscale_sif_by_crop_cdl30.R) against
#    crop-class NDVI interpolated in time to the OCO-2 day.
# 3. Drought in NDVI: the same seasonal-curve / residual analysis as
#    drought_analysis.R, on every clear HLS scene (many more dates than
#    OCO-2), tested at the level of scenes.
# 4. SIF beyond NDVI: on OCO-2 days, are SIF residuals related to drought once
#    the NDVI residual (greenness anomaly) is accounted for?
#
# Needs: hls_ndvi_lamont.R, ndvi_by_crop.R, downscale_sif_by_crop_cdl30.R 757,
#        drought_analysis.R (residuals and the study-area drought index)
# Output in ndvi_comparison/:
#   footprint_stats.csv        per OCO-2 day: correlation of each predictor with OCO-2
#   class_comparison.csv       crop-class SIF and NDVI on OCO-2 days
#   ndvi_residuals.csv         every clear class-scene: NDVI, seasonal fit, residual, drought
#   drought_stats.csv          drought tests for NDVI and for SIF given NDVI

library(terra)
library(sf)
library(data.table)
library(exactextractr)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
ndvi_dir <- file.path(proj_dir, "hls_ndvi")
res_dir  <- file.path(proj_dir, "downscaled_sif757")
out_dir  <- file.path(proj_dir, "ndvi_comparison")
oco_file <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")

max_gap_days    <- 5     # HLS scene must be within this many days of the OCO-2 day
min_scene_clear <- 0.8   # ...and at least this clear over the study area
min_fp_clear    <- 0.7   # footprints need this much clear NDVI
min_class_clear <- 0.7   # class-scenes need this much of the class clear
max_interp_days <- 16    # interpolate class NDVI only between scenes this close
analysis_classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
                      "Grassland/Pasture", "Sorghum", "Cotton")
loess_span <- 0.4
prior_sd   <- 1.0
cv_seed    <- 1

dir.create(out_dir, showWarnings = FALSE)
season_of <- function(m) fcase(m %in% c(12, 1, 2), "Winter", m %in% 3:5, "Spring",
                               m %in% 6:8, "Summer", default = "Autumn")
spearman <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.table(n = sum(ok), spearman_r = unname(t$estimate), spearman_p = t$p.value)
}

scenes  <- fread(file.path(ndvi_dir, "scenes.csv"))[, date := as.Date(date)]
by_crop <- fread(file.path(ndvi_dir, "ndvi_by_crop.csv"))[, date := as.Date(date)]
drought <- fread(file.path(proj_dir, "drought_index_study_area.csv"))[, date := as.Date(date)]
dsci_at  <- function(d) approx(as.numeric(drought$date), drought$dsci, as.numeric(d), rule = 2)$y
terciles <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))   # drought_analysis.R
dsci_class <- function(x) factor(fcase(x < terciles$low_upper, terciles$label_low,
                                       x < terciles$high_lower, terciles$label_mid,
                                       default = terciles$label_high),
                                 levels = unlist(terciles[, .(label_low, label_mid, label_high)]))

# ---- 1. Footprints ----------------------------------------------------------
design <- readRDS(file.path(res_dir, "sif757_by_crop_cdl30_design.rds"))
fp_all <- st_read(oco_file, quiet = TRUE,
                  query = "SELECT Metadata_SoundingId, geom FROM oco2_sif WHERE Metadata_MeasurementMode = 2")
fp_all$sid <- as.character(fp_all$Metadata_SoundingId)
st_geometry(fp_all) <- st_convex_hull(st_geometry(fp_all))
scene_clear <- by_crop[crop_class == "All", .(granule, date, clear_frac)]

# Weighted least squares with a weak prior (as in the crop-class inversion);
# X = crop-class fractions and/or centred NDVI
fit_predict <- function(X_tr, y_tr, se2_tr, X_te) {
  xa <- sum(y_tr / se2_tr) / sum(1 / se2_tr)
  XtSe <- t(X_tr / se2_tr)
  S <- solve(XtSe %*% X_tr + diag(1 / prior_sd^2, ncol(X_tr)))
  b <- S %*% (XtSe %*% (y_tr - xa))
  as.vector(xa + X_te %*% b)
}

set.seed(cv_seed)
fp_stats <- rbindlist(lapply(names(design), function(dn) {
  d <- as.Date(dn); ds <- design[[dn]]
  cand <- scene_clear[abs(as.numeric(date - d)) <= max_gap_days & clear_frac >= min_scene_clear]
  if (nrow(cand) == 0) return(NULL)
  sc <- cand[order(-clear_frac, abs(as.numeric(date - d)))][1]
  r <- (rast(scenes[granule == sc$granule, file]) - 100) / 100

  fp <- fp_all[match(ds$sounding_id, fp_all$sid), ]
  ex <- exact_extract(r, st_transform(fp, crs(r)), function(v, w) {
    c(mean = sum(v * w, na.rm = TRUE) / sum(w[!is.na(v)]), clear = sum(w[!is.na(v)]) / sum(w))
  }, progress = FALSE)
  ndvi <- ex[1, ]; ok <- ex[2, ] >= min_fp_clear & !is.na(ndvi)
  if (sum(ok) < 50) return(NULL)

  y <- ds$y[ok]; se2 <- ds$se2[ok]
  Kc <- ds$Kc[ok, colSums(ds$Kc[ok, , drop = FALSE]) > 0, drop = FALSE]
  N  <- matrix(ndvi[ok] - mean(ndvi[ok]), ncol = 1)
  designs <- list(crop = Kc, ndvi = N, crop_ndvi = cbind(Kc, N))
  fold <- sample(rep(1:2, length.out = length(y)))
  preds <- lapply(designs, function(X) {
    p <- numeric(length(y))
    for (k in 1:2) {
      tr <- fold != k
      p[!tr] <- fit_predict(X[tr, , drop = FALSE], y[tr], se2[tr], X[!tr, , drop = FALSE])
    }
    p
  })
  data.table(date = d, hls_date = sc$date, hls_sensor = substr(sc$granule, 5, 7),
             n = length(y), r_crop_cv = cor(y, preds$crop), r_ndvi_cv = cor(y, preds$ndvi),
             r_crop_ndvi_cv = cor(y, preds$crop_ndvi),
             noise_ceiling = sqrt(max(0, 1 - mean(se2) / var(y))))
}))
fwrite(fp_stats, file.path(out_dir, "footprint_stats.csv"))

# ---- 2. Crop classes on OCO-2 days ------------------------------------------
clear_class <- by_crop[crop_class %in% analysis_classes & clear_frac >= min_class_clear &
                       !is.na(ndvi)]
# Average same-day scenes (Landsat and Sentinel-2 on one day)
clear_class <- clear_class[, .(ndvi = mean(ndvi)), by = .(crop_class, date)]
setorder(clear_class, crop_class, date)

interp_to <- function(cls, d, value_col = "ndvi") {
  s <- clear_class[crop_class == cls]
  before <- s[date <= d][.N]; after <- s[date >= d][1]
  if (nrow(before) == 0 || nrow(after) == 0 || is.na(before$date) || is.na(after$date)) return(NA_real_)
  if (as.numeric(d - before$date) > max_interp_days || as.numeric(after$date - d) > max_interp_days)
    return(NA_real_)
  if (before$date == after$date) return(before[[value_col]])
  approx(as.numeric(c(before$date, after$date)), c(before[[value_col]], after[[value_col]]),
         as.numeric(d))$y
}

sif_cls <- fread(file.path(res_dir, "sif757_by_crop_cdl30.csv"))[
  crop_class %in% analysis_classes & avg_kernel >= 0.9][, date := as.Date(date)]
sif_cls[, ndvi := mapply(interp_to, crop_class, date)]
fwrite(sif_cls, file.path(out_dir, "class_comparison.csv"))

# ---- 3. Drought in NDVI -----------------------------------------------------
nd <- copy(clear_class)
nd[, `:=`(year = as.integer(format(date, "%Y")), doy = as.integer(format(date, "%j")),
          season = season_of(as.integer(format(date, "%m"))))]
seasonal_fit <- function(fit_data, new_doy) {
  wrap <- rbind(fit_data, copy(fit_data)[, doy := doy - 365], copy(fit_data)[, doy := doy + 365])
  predict(loess(ndvi ~ doy, data = wrap, span = loess_span / 3, degree = 2),
          newdata = data.frame(doy = new_doy))
}
nd[, ndvi_fit := {                       # leave-one-year-out, as for SIF
  out <- numeric(.N)
  for (y in unique(year)) out[year == y] <- seasonal_fit(.SD[year != y], doy[year == y])
  out
}, by = crop_class, .SDcols = c("doy", "ndvi", "year")]
nd[, residual := ndvi - ndvi_fit]
nd[, dsci := dsci_at(date)]
nd[, drought := dsci_class(dsci)]   # DSCI classes (drought_analysis.R)
fwrite(nd, file.path(out_dir, "ndvi_residuals.csv"))

# Scene-level tests (residuals of crops on the same date share weather and
# viewing effects, so test the mean across crops per date)
nd_day <- nd[, .(residual = mean(residual), dsci = dsci[1], year = year[1],
                 season = season[1]), by = date]
nd[, residual_vs_crops := residual - mean(residual), by = date]

# ---- 4. SIF beyond NDVI -----------------------------------------------------
# NDVI residual at each OCO-2 class-day, interpolated like NDVI itself
sif_res <- fread(file.path(proj_dir, "drought_residuals.csv"))[, date := as.Date(date)]
clear_class <- merge(clear_class, nd[, .(crop_class, date, residual)], by = c("crop_class", "date"))
sif_res[, ndvi_residual := mapply(interp_to, as.character(crop_class), date,
                                  MoreArgs = list(value_col = "residual"))]
sif_day <- sif_res[!is.na(ndvi_residual), .(sif_residual = mean(residual),
                                            ndvi_residual = mean(ndvi_residual),
                                            dsci = dsci[1]), by = date]
# SIF residual left over after the NDVI anomaly: residual of a linear fit
sif_day[, sif_beyond_ndvi := resid(lm(sif_residual ~ ndvi_residual))]
fwrite(sif_day, file.path(out_dir, "sif_vs_ndvi_residuals_by_day.csv"))

stats <- rbind(
  cbind(test = "NDVI: scene-mean residual vs drought index", spearman(nd_day$residual, nd_day$dsci)),
  cbind(test = "NDVI: same, without 2022-2023", nd_day[!year %in% 2022:2023, spearman(residual, dsci)]),
  nd_day[, cbind(test = paste("NDVI, scene-mean residual vs drought:", season[1]),
                 spearman(residual, dsci)), by = season][, season := NULL],
  nd[, cbind(test = paste("NDVI relative to other crops:", crop_class[1]),
             spearman(residual_vs_crops, dsci)), by = crop_class][, crop_class := NULL],
  cbind(test = "SIF residual vs NDVI residual (OCO-2 days)",
        spearman(sif_day$sif_residual, sif_day$ndvi_residual)),
  cbind(test = "SIF residual vs drought (OCO-2 days with NDVI)",
        spearman(sif_day$sif_residual, sif_day$dsci)),
  cbind(test = "SIF residual beyond NDVI vs drought",
        spearman(sif_day$sif_beyond_ndvi, sif_day$dsci)))
fwrite(stats, file.path(out_dir, "drought_stats.csv"))

# ---- Summary ----------------------------------------------------------------
cat("\nFootprints (", nrow(fp_stats), "OCO-2 days with a clear HLS scene; medians):\n")
print(fp_stats[, lapply(.SD, median), .SDcols = c("r_crop_cv", "r_ndvi_cv", "r_crop_ndvi_cv",
                                                  "noise_ceiling")], digits = 3)
cc <- sif_cls[!is.na(ndvi)]
cc[, `:=`(sif_a = sif757 - mean(sif757), ndvi_a = ndvi - mean(ndvi)), by = date]
cat("\nCrop classes:", nrow(cc), "class-days with NDVI; r(SIF, NDVI) =",
    round(cor(cc$sif757, cc$ndvi), 3), "; within-day contrasts r =",
    round(cor(cc$sif_a, cc$ndvi_a), 3), "\n")
cat("\nNDVI drought analysis:", nrow(nd_day), "clear scene dates\n")
print(stats, digits = 2)

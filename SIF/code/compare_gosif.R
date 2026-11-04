# Compare the OCO-2 crop-class SIF downscaling with GOSIF v2 (Li & Xiao 2019:
# OCO-2 SIF_757nm extended in space and time with MODIS EVI and MERRA-2
# meteorology; 8-day composites at 0.05 deg), alongside the SIFNet comparison.
#
# Both are compared at 757 nm, against the 757 nm crop-class inversion.
# GOSIF is a DAILY-MEAN SIF, while the inversion fits OCO-2's instantaneous
# Science_SIF_757nm. GOSIF is put on the instantaneous scale with OCO-2's own
# per-sounding Science_daily_correction_factor (Daily_SIF_757nm =
# Science_SIF_757nm x factor): footprint predictions use each footprint's
# factor, crop-class means use the day's median factor. Correlations and
# anomaly skill don't depend on this scaling; bias and ratio do.
# Each OCO-2 target day is matched to the GOSIF 8-day interval containing it.
#
# GOSIF pixels (~5.5 km) are larger than OCO-2 footprints (~1.3 x 2.2 km), and
# the study extent covers only about 6 x 4 of them, so GOSIF can resolve
# little of the field-scale variation the crop-class model targets. Note too
# that GOSIF was trained on OCO-2 SIF, so it is not independent of OCO-2.
#
# 1. Footprints: for every OCO-2 footprint, the SIF each product predicts
#    against the SIF OCO-2 measured (GOSIF: area-weighted mean over the
#    footprint). The crop-class prediction is two-fold cross-validated as in
#    compare_sifnet.R; a day-mean prediction is the no-skill baseline.
# 2. Crop classes: GOSIF averaged over each crop class (30 m CDL, within
#    raster_extent.shp) for every 8-day interval, against the OCO-2 crop-class
#    estimates on matched days.
# 3. Summary table of the footprint skill of GOSIF and SIFNet, on each
#    product's own days and on the days both were compared.
#
# Needs: downscale_sif_by_crop_cdl30.R run with band 757; compare_sifnet.R for
# the joint summary (skipped if its output is missing).
#
# Output (in out_dir):
#   footprint_predictions.csv  one row per footprint on matched days
#   footprint_stats.csv        one row per day: bias and skill of each product
#   gosif_class_means.csv      GOSIF mean per crop class per 8-day interval
#                              (daily-mean units)
#   class_comparison.csv       OCO-2 vs GOSIF per crop class on matched days
#                              (instantaneous units)
#   product_skill_summary.csv  GOSIF and SIFNet footprint skill side by side

library(terra)
library(sf)
library(data.table)
library(exactextractr)

# ---- Settings ---------------------------------------------------------------
proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
gosif_dir  <- file.path(proj_dir, "GOSIF")
gosif_scale <- 1e-4          # stored as integers x 10000 (W m-2 um-1 sr-1)
gosif_fill  <- 32766         # 32766 / 32767 = snow / water / no data
interval_days <- 8
oco_file   <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
extent_shp <- file.path(proj_dir, "raster_extent.shp")
res_dir    <- file.path(proj_dir, "downscaled_sif757")
design_file <- file.path(res_dir, "sif757_by_crop_cdl30_design.rds")
class_file  <- file.path(res_dir, "sif757_by_crop_cdl30.csv")
cdl_dir    <- file.path(proj_dir, "cdl_30m")
sifnet_stats <- file.path(proj_dir, "sifnet_comparison", "footprint_stats.csv")
out_dir    <- file.path(proj_dir, "gosif_comparison")
prior_sd   <- 1.0     # as in downscale_sif_by_crop_cdl30.R
cv_seed    <- 1

dir.create(out_dir, showWarnings = FALSE)
source(file.path(proj_dir, "crop_classes.R"))   # class_levels, class_raster()

# ---- GOSIF intervals ----------------------------------------------------------
gosif_files <- sort(list.files(gosif_dir, "^GOSIF_[0-9]{7}[.]tif([.]gz)?$", full.names = TRUE))
gosif <- data.table(file = gosif_files,
                    start = as.Date(sub("^GOSIF_([0-9]{7}).*", "\\1", basename(gosif_files)), "%Y%j"))
setorder(gosif, start)
# 8-day composites; the last one of each year runs to 31 Dec
gosif[, end := pmin(start + interval_days - 1, as.Date(paste0(format(start, "%Y"), "-12-31")))]
interval_of <- function(d) gosif[start <= d & end >= d]

extent <- vect(st_read(extent_shp, quiet = TRUE))
ext_ll <- ext(project(extent, "EPSG:4326")) + 0.1     # pad so edge footprints are covered
read_gosif <- function(f) {                            # daily-mean SIF757, W m-2 um-1 sr-1
  src <- if (grepl("[.]gz$", f)) paste0("/vsigzip/", normalizePath(f)) else f
  r <- crop(rast(src), ext_ll)
  r <- ifel(r >= gosif_fill, NA, r) * gosif_scale
  names(r) <- "gosif"
  r
}

# ---- Matched OCO-2 days -----------------------------------------------------
design <- readRDS(design_file)
days <- as.Date(names(design))
days <- days[vapply(days, function(d) nrow(interval_of(d)) == 1, logical(1))]
cat(length(days), "OCO-2 target days fall within a GOSIF interval\n")

fp_all <- st_read(oco_file, quiet = TRUE,
                  query = paste("SELECT Metadata_SoundingId, Science_daily_correction_factor, geom",
                                "FROM oco2_sif WHERE Metadata_MeasurementMode = 2"))
fp_all$sid <- as.character(fp_all$Metadata_SoundingId)
st_geometry(fp_all) <- st_convex_hull(st_geometry(fp_all))

# Bayesian crop-class fit (same as downscale_sif_by_crop_cdl30.R)
fit_classes <- function(Kc, y, se2) {
  xa <- sum(y / se2) / sum(1 / se2)
  KtSe <- t(Kc / se2)
  S_hat <- solve(KtSe %*% Kc + diag(1 / prior_sd^2, ncol(Kc)))
  list(x = as.vector(xa + S_hat %*% (KtSe %*% (y - Kc %*% rep(xa, ncol(Kc))))),
       mean = xa)
}

# ---- 1. Footprint comparison ------------------------------------------------
set.seed(cv_seed)
fp_rows <- rbindlist(lapply(days, function(d) {
  ds <- design[[as.character(d)]]
  keep_cols <- colSums(ds$Kc) > 0
  Kc <- ds$Kc[, keep_cols, drop = FALSE]
  n <- length(ds$y)

  fold <- sample(rep(1:2, length.out = n))
  pred_crop <- pred_mean <- numeric(n)
  for (k in 1:2) {
    tr <- fold != k
    fit <- fit_classes(Kc[tr, , drop = FALSE], ds$y[tr], ds$se2[tr])
    pred_crop[!tr] <- as.vector(Kc[!tr, , drop = FALSE] %*% fit$x)
    pred_mean[!tr] <- fit$mean
  }

  fp <- fp_all[match(ds$sounding_id, fp_all$sid), ]
  iv <- interval_of(d)
  r <- read_gosif(iv$file)
  fp_ll <- st_transform(fp, crs(r))
  gosif_daily <- exact_extract(r, fp_ll, "mean", progress = FALSE)
  xy <- st_coordinates(st_centroid(st_geometry(fp_ll)))

  data.table(date = d, gosif_start = iv$start,
             sounding_id = ds$sounding_id, obs = ds$y, obs_unc = sqrt(ds$se2),
             daily_factor = fp$Science_daily_correction_factor,
             pred_crop_cv = pred_crop, pred_day_mean = pred_mean,
             gosif_cell = cellFromXY(r, xy), gosif_daily = gosif_daily,
             pred_gosif = gosif_daily / fp$Science_daily_correction_factor)
}))
fwrite(fp_rows, file.path(out_dir, "footprint_predictions.csv"))

# Same statistics as compare_sifnet.R. n_gosif_px = GOSIF pixels containing that
# day's footprint centres (how much spatial detail GOSIF can offer at all).
footprint_stats <- fp_rows[!is.na(pred_gosif), {
  anom_mse <- function(p) mean(((obs - mean(obs)) - (p - mean(p)))^2)
  .(n = .N,
    n_gosif_px = uniqueN(gosif_cell),
    obs_mean = mean(obs), gosif_mean = mean(pred_gosif),
    gosif_bias = mean(pred_gosif - obs),
    gosif_ratio = mean(pred_gosif) / mean(obs),
    r_gosif = if (sd(pred_gosif) > 0) cor(obs, pred_gosif) else NA_real_,
    r_crop_cv = cor(obs, pred_crop_cv),
    noise_ceiling = sqrt(max(0, 1 - mean(obs_unc^2) / var(obs))),
    skill_gosif = 1 - anom_mse(pred_gosif) / anom_mse(pred_day_mean),
    skill_crop_cv = 1 - anom_mse(pred_crop_cv) / anom_mse(pred_day_mean))
}, by = .(date, gosif_start)]
fwrite(footprint_stats, file.path(out_dir, "footprint_stats.csv"))

# ---- 2. Crop-class means of GOSIF -------------------------------------------
# Every 8-day interval that has a CDL year (for the time-series figure).
gosif_class <- rbindlist(lapply(seq_len(nrow(gosif)), function(k) {
  mid <- gosif$start[k] + interval_days %/% 2
  cls <- tryCatch(class_raster(mid), error = function(e) NULL)
  if (is.null(cls)) return(NULL)
  ext_cdl <- project(extent, crs(cls))
  cls <- mask(crop(cls, ext_cdl), ext_cdl)
  # GOSIF on the 30 m CDL grid (nearest neighbour), then mean per class:
  # an area-weighted mean of GOSIF over each class's pixels
  g30 <- project(read_gosif(gosif$file[k]), cls, method = "near")
  z <- zonal(g30, cls, fun = "mean", na.rm = TRUE)
  data.table(gosif_start = gosif$start[k], mid_date = mid,
             crop_class = class_levels[z[[1]]], gosif_daily_sif757 = z[[2]])
}))
fwrite(gosif_class, file.path(out_dir, "gosif_class_means.csv"))

# Matched days, GOSIF scaled to instantaneous with the day's median factor
day_factor <- fp_rows[, .(daily_factor = median(daily_factor, na.rm = TRUE)), by = date]
oco_class <- fread(class_file)[, date := as.Date(date)][date %in% days]
oco_class[, gosif_start := as.Date(vapply(date, function(d) as.numeric(interval_of(d)$start), numeric(1)))]
class_comp <- merge(oco_class, gosif_class, by = c("gosif_start", "crop_class"))
class_comp <- merge(class_comp, day_factor, by = "date")
class_comp[, gosif_sif757 := gosif_daily_sif757 / daily_factor]
fwrite(class_comp, file.path(out_dir, "class_comparison.csv"))

# ---- 3. GOSIF vs SIFNet footprint skill ---------------------------------------
med <- function(dt, r, skill, ratio, label) {
  dt[, .(product = label, days = .N, median_r = median(get(r), na.rm = TRUE),
         median_r_crop_cv = median(r_crop_cv), median_noise_ceiling = median(noise_ceiling),
         median_skill = median(get(skill), na.rm = TRUE),
         median_skill_crop_cv = median(skill_crop_cv),
         median_ratio = median(get(ratio), na.rm = TRUE))]
}
summ <- list(med(footprint_stats, "r_gosif", "skill_gosif", "gosif_ratio", "GOSIF (757 nm), all days"))
if (file.exists(sifnet_stats)) {
  sn <- fread(sifnet_stats)[, date := as.Date(date)]
  summ <- c(summ, list(med(sn, "r_sifnet", "skill_sifnet", "sifnet_ratio", "SIFNet (740 nm), all days")))
  common <- intersect(footprint_stats$date, sn$date)
  if (length(common)) {
    summ <- c(summ, list(
      med(footprint_stats[date %in% common], "r_gosif", "skill_gosif", "gosif_ratio",
          "GOSIF (757 nm), days compared with both"),
      med(sn[date %in% common], "r_sifnet", "skill_sifnet", "sifnet_ratio",
          "SIFNet (740 nm), days compared with both")))
  }
} else {
  cat("(no", sifnet_stats, "- run compare_sifnet.R for the joint summary)\n")
}
summ <- rbindlist(summ)
fwrite(summ, file.path(out_dir, "product_skill_summary.csv"))

# ---- Summary ----------------------------------------------------------------
cat("\nFootprint level (", nrow(footprint_stats), "days; medians):\n")
print(footprint_stats[, .(r_gosif = median(r_gosif, na.rm = TRUE), r_crop_cv = median(r_crop_cv),
                          noise_ceiling = median(noise_ceiling),
                          skill_gosif = median(skill_gosif), skill_crop_cv = median(skill_crop_cv),
                          gosif_ratio = median(gosif_ratio),
                          gosif_pixels_per_day = median(n_gosif_px))], digits = 3)
cat("Days where GOSIF correlates better with OCO-2 than the crop-class model:",
    sum(footprint_stats$r_gosif > footprint_stats$r_crop_cv, na.rm = TRUE), "of",
    nrow(footprint_stats), "\n")

cc <- class_comp[avg_kernel >= 0.9 & !crop_class %in% c("Developed", "Other")]
cat("\nCrop classes on matched days (avg_kernel >= 0.9, excl. Developed/Other):\n")
cat("  r(OCO-2, GOSIF) across all class-days:", round(cor(cc$sif757, cc$gosif_sif757), 3), "\n")
cc[, `:=`(oco_anom = sif757 - mean(sif757), gosif_anom = gosif_sif757 - mean(gosif_sif757)), by = date]
cat("  r of crop contrasts within a day (day means removed):",
    round(cor(cc$oco_anom, cc$gosif_anom), 3), "\n")
cat("  spread of crop contrasts (SD of class minus day mean): OCO-2",
    round(sd(cc$oco_anom), 3), " GOSIF", round(sd(cc$gosif_anom), 3), "\n")
cat("  GOSIF / OCO-2 ratio (median):", round(median(cc$gosif_sif757 / cc$sif757), 2), "\n")

cat("\nGOSIF vs SIFNet (footprint level, medians):\n")
print(summ, digits = 3)

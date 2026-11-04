# Compare the OCO-2 crop-class SIF downscaling with SIFNet (TROPOMI SIF at
# 740 nm downscaled with a CNN and NDVI; 16-day composites at 0.005 deg).
#
# Both are compared at 740 nm: OCO-2 SIF_740nm from the Lite product against
# SIFNet "sif" (overpass-time SIF; "sif_dc" is the length-of-day corrected
# daily value). OCO-2 and TROPOMI both pass over near 13:30 local time.
# Each OCO-2 target day is matched to the SIFNet interval that contains it.
#
# 1. Footprints: for every OCO-2 footprint, the SIF each product predicts
#    against the SIF OCO-2 measured. SIFNet's prediction is its area-weighted
#    mean over the footprint. The crop-class model's prediction is
#    cross-validated: fitted to a random half of the day's footprints, then
#    used to predict the other half (and vice versa). A day-mean prediction is
#    the no-skill baseline.
# 2. Crop classes: SIFNet averaged over each crop class (30 m CDL, within
#    raster_extent.shp) for every 16-day interval, against the OCO-2
#    crop-class estimates on matched days.
#
# Needs: downscale_sif_by_crop_cdl30.R run with band 740 (its *_design.rds
# holds each day's footprints and Jacobian).
#
# Output (in out_dir):
#   footprint_predictions.csv  one row per footprint on matched days
#   footprint_stats.csv        one row per day: bias and skill of each product
#   sifnet_class_means.csv     SIFNet mean per crop class per 16-day interval
#   class_comparison.csv       OCO-2 vs SIFNet per crop class on matched days

library(terra)
library(sf)
library(data.table)
library(exactextractr)

# ---- Settings ---------------------------------------------------------------
proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
sifnet_dir <- file.path(proj_dir, "SIFNet_OK")
sifnet_var <- "sif"
interval_days <- 16
oco_file   <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
extent_shp <- file.path(proj_dir, "raster_extent.shp")
res_dir    <- file.path(proj_dir, "downscaled_sif757")
design_file <- file.path(res_dir, "sif740_by_crop_cdl30_design.rds")
class_file  <- file.path(res_dir, "sif740_by_crop_cdl30.csv")
cdl_dir    <- file.path(proj_dir, "cdl_30m")
out_dir    <- file.path(proj_dir, "sifnet_comparison")
prior_sd   <- 1.0     # as in downscale_sif_by_crop_cdl30.R
cv_seed    <- 1

dir.create(out_dir, showWarnings = FALSE)
source(file.path(proj_dir, "crop_classes.R"))   # class_levels, class_raster()

# ---- SIFNet intervals -------------------------------------------------------
sifnet_files <- sort(list.files(sifnet_dir, "^SIFnet_OK_SIF_.*[.]nc$",
                                recursive = TRUE, full.names = TRUE))
sifnet <- data.table(file = sifnet_files,
                     start = as.Date(sub(".*_([0-9-]{10})[.]nc$", "\\1", sifnet_files)))
setorder(sifnet, start)
sifnet[, end := start + interval_days - 1]
read_sifnet <- function(f) rast(f, subds = sifnet_var)
interval_of <- function(d) sifnet[start <= d & end >= d]

# ---- Matched OCO-2 days -----------------------------------------------------
design <- readRDS(design_file)
days <- as.Date(names(design))
days <- days[vapply(days, function(d) nrow(interval_of(d)) == 1, logical(1))]
cat(length(days), "OCO-2 target days fall within a SIFNet interval\n")

fp_all <- st_read(oco_file, quiet = TRUE,
                  query = "SELECT Metadata_SoundingId, geom FROM oco2_sif WHERE Metadata_MeasurementMode = 2")
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

  # Two-fold cross-validation of the crop-class model
  fold <- sample(rep(1:2, length.out = n))
  pred_crop <- pred_mean <- numeric(n)
  for (k in 1:2) {
    tr <- fold != k
    fit <- fit_classes(Kc[tr, , drop = FALSE], ds$y[tr], ds$se2[tr])
    pred_crop[!tr] <- as.vector(Kc[!tr, , drop = FALSE] %*% fit$x)
    pred_mean[!tr] <- fit$mean
  }

  # SIFNet averaged over each footprint
  fp <- fp_all[match(ds$sounding_id, fp_all$sid), ]
  r <- read_sifnet(interval_of(d)$file)
  pred_sifnet <- exact_extract(r, st_transform(fp, crs(r)), "mean", progress = FALSE)

  data.table(date = d, sifnet_start = interval_of(d)$start,
             sounding_id = ds$sounding_id, obs = ds$y, obs_unc = sqrt(ds$se2),
             pred_crop_cv = pred_crop, pred_day_mean = pred_mean,
             pred_sifnet = pred_sifnet)
}))
fwrite(fp_rows, file.path(out_dir, "footprint_predictions.csv"))

# Skill per day. Correlation and the anomaly MSE ignore each product's overall
# level, so they measure how well it places SIF in space; "bias" is the level
# difference. noise_ceiling = the best correlation possible given the OCO-2
# retrieval noise.
footprint_stats <- fp_rows[!is.na(pred_sifnet), {
  anom_mse <- function(p) mean(((obs - mean(obs)) - (p - mean(p)))^2)
  .(n = .N,
    obs_mean = mean(obs), sifnet_mean = mean(pred_sifnet),
    sifnet_bias = mean(pred_sifnet - obs),
    sifnet_ratio = mean(pred_sifnet) / mean(obs),
    r_sifnet = cor(obs, pred_sifnet),
    r_crop_cv = cor(obs, pred_crop_cv),
    noise_ceiling = sqrt(max(0, 1 - mean(obs_unc^2) / var(obs))),
    skill_sifnet = 1 - anom_mse(pred_sifnet) / anom_mse(pred_day_mean),
    skill_crop_cv = 1 - anom_mse(pred_crop_cv) / anom_mse(pred_day_mean))
}, by = .(date, sifnet_start)]
fwrite(footprint_stats, file.path(out_dir, "footprint_stats.csv"))

# ---- 2. Crop-class means of SIFNet ------------------------------------------
extent <- vect(st_read(extent_shp, quiet = TRUE))
sifnet_class <- rbindlist(lapply(seq_len(nrow(sifnet)), function(k) {
  mid <- sifnet$start[k] + interval_days %/% 2
  cls <- tryCatch(class_raster(mid), error = function(e) NULL)
  if (is.null(cls)) return(NULL)
  ext_cdl <- project(extent, crs(cls))
  cls <- mask(crop(cls, ext_cdl), ext_cdl)
  # SIFNet on the 30 m CDL grid (nearest neighbour), then mean per class:
  # an area-weighted mean of SIFNet over each class's pixels
  s30 <- project(read_sifnet(sifnet$file[k]), cls, method = "near")
  z <- zonal(s30, cls, fun = "mean", na.rm = TRUE)
  data.table(sifnet_start = sifnet$start[k], mid_date = mid,
             crop_class = class_levels[z[[1]]], sifnet_sif740 = z[[2]])
}))
fwrite(sifnet_class, file.path(out_dir, "sifnet_class_means.csv"))

oco_class <- fread(class_file)[, date := as.Date(date)][date %in% days]
oco_class[, sifnet_start := as.Date(vapply(date, function(d) as.numeric(interval_of(d)$start), numeric(1)))]
class_comp <- merge(oco_class, sifnet_class, by = c("sifnet_start", "crop_class"))
fwrite(class_comp, file.path(out_dir, "class_comparison.csv"))

# ---- Summary ----------------------------------------------------------------
cat("\nFootprint level (", nrow(footprint_stats), "days; medians):\n")
print(footprint_stats[, .(r_sifnet = median(r_sifnet), r_crop_cv = median(r_crop_cv),
                          noise_ceiling = median(noise_ceiling),
                          skill_sifnet = median(skill_sifnet),
                          skill_crop_cv = median(skill_crop_cv),
                          sifnet_ratio = median(sifnet_ratio))], digits = 3)
cat("Days where SIFNet correlates better with OCO-2 than the crop-class model:",
    sum(footprint_stats$r_sifnet > footprint_stats$r_crop_cv), "of",
    nrow(footprint_stats), "\n")

cc <- class_comp[avg_kernel >= 0.9 & !crop_class %in% c("Developed", "Other")]
cat("\nCrop classes on matched days (avg_kernel >= 0.9, excl. Developed/Other):\n")
cat("  r(OCO-2, SIFNet) across all class-days:", round(cor(cc$sif740, cc$sifnet_sif740), 3), "\n")
cc[, `:=`(oco_anom = sif740 - mean(sif740), sifnet_anom = sifnet_sif740 - mean(sifnet_sif740)), by = date]
cat("  r of crop contrasts within a day (day means removed):",
    round(cor(cc$oco_anom, cc$sifnet_anom), 3), "\n")
cat("  SIFNet / OCO-2 ratio (median):", round(median(cc$sifnet_sif740 / cc$sif740), 2), "\n")

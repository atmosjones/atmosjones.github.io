# Mean HLS NDVI per crop class in the study area (raster_extent.shp) for
# every HLS scene, using the same 30 m CDL crop classes as the SIF inversion
# (crop_classes.R, including the fall winter-wheat rule).
#
# Input:  hls_ndvi/scenes.csv and scenes/*.tif (hls_ndvi_lamont.R), cdl_30m/
# Output: hls_ndvi/ndvi_by_crop.csv, one row per scene and class:
#           ndvi        mean NDVI over the class's clear pixels
#           clear_frac  fraction of the class's pixels that were clear
#         crop_class "All" = the whole study area (also for years without a CDL)

library(terra)
library(sf)
library(data.table)
library(parallel)

proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
cdl_dir    <- file.path(proj_dir, "cdl_30m")
ndvi_dir   <- file.path(proj_dir, "hls_ndvi")
extent_shp <- file.path(proj_dir, "raster_extent.shp")
n_workers  <- 4

source(file.path(proj_dir, "crop_classes.R"))   # class_levels, class_raster(), cdl

scenes <- fread(file.path(ndvi_dir, "scenes.csv"))[, date := as.Date(date)]
study <- project(vect(st_read(extent_shp, quiet = TRUE)), cdl_crs)
template <- mask(crop(cdl[[length(cdl)]], study), study)   # 30 m grid of the study area
in_study <- !is.na(values(template)[, 1])

scene_stats <- function(k) {
  s <- scenes[k]
  # NDVI on the CDL grid (nearest neighbour; both are 30 m)
  ndvi <- project((rast(s$file) - 100) / 100, template, method = "near")
  v <- values(ndvi)[, 1][in_study]
  all <- data.table(crop_class = "All", ndvi = mean(v, na.rm = TRUE),
                    clear_frac = mean(!is.na(v)), n_pixels = length(v))
  cls <- tryCatch(values(mask(crop(class_raster(s$date), study), study))[, 1][in_study],
                  error = function(e) NULL)       # no CDL for this year yet
  by_class <- if (is.null(cls)) NULL else
    data.table(c = cls, v = v)[!is.na(c), .(ndvi = mean(v, na.rm = TRUE),
                                            clear_frac = mean(!is.na(v)), n_pixels = .N), by = c][
      , .(crop_class = class_levels[c], ndvi, clear_frac, n_pixels)]
  cbind(granule = s$granule, sensor = s$sensor, date = s$date, rbind(all, by_class))
}

out <- rbindlist(mclapply(seq_len(nrow(scenes)), scene_stats, mc.cores = n_workers))
out[is.nan(ndvi), ndvi := NA]
setorder(out, date, granule, crop_class)
fwrite(out, file.path(ndvi_dir, "ndvi_by_crop.csv"))
cat(nrow(scenes), "scenes;", out[crop_class == "All" & clear_frac >= 0.7, .N],
    "with at least 70% of the study area clear\n")

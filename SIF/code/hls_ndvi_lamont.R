# NDVI over the study area from Harmonized Landsat Sentinel-2 (HLS v2.0,
# 30 m), every scene from 2014 to now, for comparison with the OCO-2 SIF.
#
# Each scene's red and NIR bands and quality mask are read remotely for the
# study-area window only; NDVI is masked for cloud, cloud shadow, adjacent
# cloud, snow, water and high aerosol, and saved as one small GeoTIFF per
# scene (bytes: NDVI x 100 + 100, 255 = masked). Reruns skip scenes already
# saved.
#   Red / NIR: Landsat 8-9 (L30) B04 / B05, Sentinel-2 (S30) B04 / B8A.
#
# Needs a NASA Earthdata login in ~/.netrc.
#
# Output: hls_ndvi/scenes/<granule>.tif, hls_ndvi/scenes.csv (scene list)

library(terra)
library(sf)
library(httr)
library(jsonlite)
library(data.table)
library(parallel)

proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir    <- file.path(proj_dir, "hls_ndvi")
scene_dir  <- file.path(out_dir, "scenes")
cdl_file   <- file.path(proj_dir, "cdl_30m", "CDL_2019.tif")  # sets the window
start_date <- "2014-01-01"
end_date   <- format(Sys.Date())
max_cloud  <- 80     # skip scenes with more cloud than this over the whole tile (%)
n_workers  <- 4

dir.create(scene_dir, recursive = TRUE, showWarnings = FALSE)
cookie <- file.path(out_dir, ".earthdata_cookies")
setGDALconfig(c("GDAL_HTTP_NETRC=YES",
                paste0("GDAL_HTTP_COOKIEFILE=", cookie),
                paste0("GDAL_HTTP_COOKIEJAR=", cookie),
                "GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR",
                "GDAL_HTTP_MAX_RETRY=5", "GDAL_HTTP_RETRY_DELAY=5"))

# Window: the CDL extent (all OCO-2 footprints plus ~2 km), in lon/lat for
# the search and in the tile's UTM projection for reading
win_5070 <- ext(rast(cdl_file))
e <- as.vector(win_5070)   # xmin, xmax, ymin, ymax
win_ll <- st_bbox(st_transform(st_as_sfc(st_bbox(
  c(xmin = e[[1]], ymin = e[[3]], xmax = e[[2]], ymax = e[[4]]), crs = st_crs(5070))), 4326))

# ---- Find scenes ------------------------------------------------------------
search_hls <- function(short_name) {
  out <- list(); page <- 1
  repeat {
    r <- GET("https://cmr.earthdata.nasa.gov/search/granules.json",
             query = list(short_name = short_name, version = "2.0",
                          bounding_box = paste(win_ll, collapse = ","),
                          temporal = paste0(start_date, "T00:00:00Z,", end_date, "T23:59:59Z"),
                          page_size = 2000, page_num = page))
    stop_for_status(r)
    e <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)$feed$entry
    if (length(e) == 0) break
    out[[page]] <- rbindlist(lapply(e, function(x) {
      urls <- vapply(x$links, `[[`, "", "href")
      data.table(granule = x$producer_granule_id, cloud = as.numeric(x$cloud_cover),
                 urls = list(urls[grepl("^https.*[.]tif$", urls)]))
    }))
    if (length(e) < 2000) break
    page <- page + 1
  }
  rbindlist(out)
}
scenes <- rbind(search_hls("HLSL30"), search_hls("HLSS30"))
scenes[, sensor := substr(granule, 5, 7)]
scenes[, tile := sub(".*[.](T[0-9A-Z]{5})[.].*", "\\1", granule)]
scenes[, date := as.Date(substr(sub(".*[.]T[0-9A-Z]{5}[.]", "", granule), 1, 7), "%Y%j")]
scenes <- scenes[cloud <= max_cloud]
scenes[, file := file.path(scene_dir, paste0(granule, ".tif"))]
cat(nrow(scenes), "scenes (", paste(scenes[, .N, by = sensor][, paste(sensor, N)], collapse = ", "),
    "); tiles:", paste(unique(scenes$tile), collapse = ", "), "\n")

# ---- NDVI per scene ---------------------------------------------------------
# Fmask bits: 1 cloud, 2 adjacent to cloud/shadow, 3 shadow, 4 snow/ice,
# 5 water, 6-7 aerosol level (3 = high)
scene_ndvi <- function(i) {
  s <- scenes[i]
  if (file.exists(s$file)) return(TRUE)
  urls <- s$urls[[1]]
  pick <- function(b) paste0("/vsicurl/", urls[grepl(paste0("[.]", b, "[.]tif$"), urls)])
  bands <- if (s$sensor == "S30") c("B04", "B8A") else c("B04", "B05")
  tryCatch({
    red_full <- rast(pick(bands[1]))
    win <- ext(project(as.polygons(rast(win_5070, crs = "EPSG:5070")), crs(red_full)))
    red <- crop(red_full, win)
    nir <- crop(rast(pick(bands[2])), win)
    fm  <- as.integer(values(crop(rast(pick("Fmask")), win)))
    r <- values(red)[, 1]; n <- values(nir)[, 1]
    ndvi <- (n - r) / (n + r)
    bad <- is.na(fm) | bitwAnd(fm, 62L) > 0 | bitwAnd(bitwShiftR(fm, 6L), 3L) == 3L |
           is.na(ndvi) | r <= 0 | n <= 0
    out <- round(pmax(pmin(ndvi, 1), -1) * 100) + 100
    out[bad] <- 255
    tmp <- paste0(s$file, ".part")
    writeRaster(setValues(rast(red), out), tmp, datatype = "INT1U", NAflag = 255,
                overwrite = TRUE, filetype = "GTiff",
                gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    file.rename(tmp, s$file)
  }, error = function(e) { message(s$granule, " failed: ", conditionMessage(e)); FALSE })
}

todo <- which(!file.exists(scenes$file))
cat("Processing", length(todo), "scenes\n")
t0 <- Sys.time()
ok <- unlist(mclapply(todo, scene_ndvi, mc.cores = n_workers, mc.preschedule = FALSE))
cat(sum(ok), "of", length(todo), "done in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins"))), "min\n")

fwrite(scenes[file.exists(file), .(granule, sensor, tile, date, cloud, file)],
       file.path(out_dir, "scenes.csv"))

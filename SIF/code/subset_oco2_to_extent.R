# Pull every OCO-2 Lite SIF sounding whose footprint intersects
# raster_extent.shp out of the daily Lite SIF files and save them as one
# GeoPackage. Two sources:
#   - v11r (Sep 2014 - Mar 2024): gs://YOUR-BUCKET/data/oco2_sif/all/
#   - v11.2r (Apr 2024 on): NASA GES DISC, found through the CMR search API.
#     v11.2r continues the record (NASA has not reprocessed earlier years in
#     it); the lite_build column (e.g. B11012Ar, B11218Ar) records which
#     version each sounding came from.
#
# The files (~70 GB) are downloaded in batches to a temporary folder,
# subset, and deleted, so they never all need to fit on disk. Per-file
# results are kept in parts_dir so an interrupted run picks up where it
# stopped. All soundings are kept (every MeasurementMode) with every
# per-sounding variable in the file; footprints are the 4-corner polygons.
#
# Needs: gcloud CLI logged in (gcloud auth login) for the bucket, and a NASA
# Earthdata login in ~/.netrc for GES DISC. Bucket files are fetched with
# curl from the Cloud Storage JSON API: `gcloud storage cp` of many files
# hangs on macOS after downloading, leaving *.gstmp files behind.

library(ncdf4)
library(sf)
library(data.table)
library(parallel)
library(httr)
library(jsonlite)

# ---- Settings ---------------------------------------------------------------
proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
extent_shp <- file.path(proj_dir, "raster_extent.shp")
bucket_dir <- "gs://YOUR-BUCKET/data/oco2_sif/all/"
gesdisc_short_name <- "OCO2_L2_Lite_SIF"
gesdisc_version    <- "11.2r"
gesdisc_start      <- "2024-04-01"   # bucket (v11r) covers the record before this
out_file   <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
parts_dir  <- file.path(proj_dir, ".oco2_subset_parts")
batch_size <- 50
n_workers  <- 4   # parallel subsetting
n_download <- 3   # parallel downloads; more streams were slower on this network

# Planar geometry: some footprints' corners are stored in an order that makes
# the polygon cross itself, which the spherical (s2) engine rejects; over a
# ~25 km box treating lon/lat as flat is accurate enough
sf_use_s2(FALSE)

dir.create(parts_dir, showWarnings = FALSE)
dl_dir <- file.path(tempdir(), "oco2_batch")

extent <- st_transform(st_read(extent_shp, quiet = TRUE), 4326)
ext_geom <- st_union(st_geometry(extent))
bb <- st_bbox(ext_geom)

# ---- Subset one file --------------------------------------------------------
subset_file <- function(nc_path) {
  nc <- nc_open(nc_path, suppress_dimvals = TRUE)
  on.exit(nc_close(nc))

  lat_c <- ncvar_get(nc, "Geolocation/footprint_latitude_vertices")
  lon_c <- ncvar_get(nc, "Geolocation/footprint_longitude_vertices")
  n <- ncol(lat_c)

  # Cheap bounding-box test on the corners before building polygons
  cand <- which(colSums(lat_c >= bb$ymin) > 0 & colSums(lat_c <= bb$ymax) > 0 &
                colSums(lon_c >= bb$xmin) > 0 & colSums(lon_c <= bb$xmax) > 0)
  cand <- cand[colSums(is.na(lat_c[, cand, drop = FALSE])) == 0 &
               colSums(is.na(lon_c[, cand, drop = FALSE])) == 0]
  if (length(cand) == 0) return(NULL)

  polys <- st_sfc(lapply(cand, function(i) {
    st_polygon(list(cbind(c(lon_c[, i], lon_c[1, i]), c(lat_c[, i], lat_c[1, i]))))
  }), crs = 4326)
  hit <- lengths(st_intersects(polys, ext_geom)) > 0
  if (!any(hit)) return(NULL)
  idx <- cand[hit]

  # Every variable with one value per sounding
  per_sounding <- Filter(function(v) {
    d <- nc$var[[v]]$varsize
    length(d) == 1 && d == n
  }, names(nc$var))
  cols <- lapply(per_sounding, function(v) {
    x <- ncvar_get(nc, v, start = min(idx), count = max(idx) - min(idx) + 1)
    x[idx - min(idx) + 1]
  })
  names(cols) <- gsub("/", "_", per_sounding)
  d <- as.data.table(cols)

  # Exact UTC from the sounding ID (YYYYMMDDhhmmss + two digits); time_tai93
  # counts from 1993-01-01 and includes leap seconds
  sid <- format(d$Metadata_SoundingId, scientific = FALSE, trim = TRUE)
  d[, time_utc := as.POSIXct(substr(sid, 1, 14), format = "%Y%m%d%H%M%S", tz = "UTC")]
  d[, Metadata_SoundingId := sid]
  d[, source_file := basename(nc_path)]

  st_sf(d, geometry = polys[hit])
}

# ---- Batches ----------------------------------------------------------------
listing <- suppressWarnings(system2("gcloud", c("storage", "ls", bucket_dir),
                                    stdout = TRUE, stderr = TRUE))
files <- grep("^gs://.*[.]nc4$", listing, value = TRUE)
if (length(files) == 0) {
  stop("Could not list ", bucket_dir, " - is the gcloud login still valid? ",
       "Run `gcloud auth login` and try again.\n", paste(listing, collapse = "\n"))
}
# v11.2r daily files at GES DISC, listed through NASA's CMR search API
cmr_files <- function() {
  out <- character(); page <- 1
  repeat {
    r <- GET("https://cmr.earthdata.nasa.gov/search/granules.json",
             query = list(short_name = gesdisc_short_name, version = gesdisc_version,
                          temporal = paste0(gesdisc_start, "T00:00:00Z,"),
                          page_size = 2000, page_num = page))
    stop_for_status(r)
    e <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)$feed$entry
    if (length(e) == 0) break
    out <- c(out, vapply(e, function(g) {
      h <- vapply(g$links, `[[`, "", "href")
      h[grepl("^https://.*[.]nc4$", h)][1]
    }, ""))
    if (length(e) < 2000) break
    page <- page + 1
  }
  out[!is.na(out)]
}
files <- c(files, cmr_files())
cat(sum(grepl("^gs://", files)), "files in the bucket,",
    sum(grepl("^https://", files)), "at GES DISC\n")

part_of <- function(f) file.path(parts_dir, sub("[.]nc4$", ".rds", basename(f)))
todo <- files[!file.exists(part_of(files))]
cat(length(files), "files;", length(todo), "left to process\n")

# Download one gs:// object with curl; returns TRUE on success
gcs_get <- function(url, dest_dir, token) {
  bucket <- sub("^gs://([^/]+)/.*", "\\1", url)
  object <- utils::URLencode(sub("^gs://[^/]+/", "", url), reserved = TRUE)
  dest <- file.path(dest_dir, basename(url))
  status <- system2("curl", c("-sf", "--retry", "3", "--max-time", "900",
                              "-H", shQuote(paste("Authorization: Bearer", token)),
                              "-o", shQuote(paste0(dest, ".part")),
                              shQuote(sprintf("https://storage.googleapis.com/storage/v1/b/%s/o/%s?alt=media",
                                              bucket, object))))
  ok <- status == 0 && file.rename(paste0(dest, ".part"), dest)
  if (!ok) unlink(paste0(dest, ".part"))
  ok
}

# Download one file from NASA GES DISC (Earthdata login from ~/.netrc)
cookie_file <- file.path(tempdir(), "earthdata_cookies")
nasa_get <- function(url, dest_dir) {
  dest <- file.path(dest_dir, basename(url))
  status <- system2("curl", c("-sfL", "-n", "--retry", "3", "--max-time", "900",
                              "-c", shQuote(cookie_file), "-b", shQuote(cookie_file),
                              "-o", shQuote(paste0(dest, ".part")), shQuote(url)))
  ok <- status == 0 && file.rename(paste0(dest, ".part"), dest)
  if (!ok) unlink(paste0(dest, ".part"))
  ok
}

batches <- split(todo, ceiling(seq_along(todo) / batch_size))
t0 <- Sys.time()
for (b in seq_along(batches)) {
  unlink(dl_dir, recursive = TRUE); dir.create(dl_dir)
  token <- NULL
  if (any(grepl("^gs://", batches[[b]]))) {
    token <- suppressWarnings(system2("gcloud", c("auth", "print-access-token"),
                                      stdout = TRUE, stderr = FALSE))  # ~1 h lifetime
    if (length(token) != 1 || !grepl("^ya29[.]", token)) {
      stop("gcloud could not provide an access token (login expired?). ",
           "Run `gcloud auth login`, then rerun this script to continue.")
    }
  }
  got <- unlist(mclapply(batches[[b]], function(u) {
    if (grepl("^gs://", u)) gcs_get(u, dl_dir, token) else nasa_get(u, dl_dir)
  }, mc.cores = n_download, mc.preschedule = FALSE))
  if (!all(got)) message("  ", sum(!got), " downloads failed in batch ", b, " (retried on rerun)")
  if (!any(got)) stop("Every download in batch ", b, " failed; stopping. ",
                      "Check the network, `gcloud auth login` and ~/.netrc, then rerun to continue.")

  local <- file.path(dl_dir, basename(batches[[b]]))
  res <- mclapply(local, function(f) {
    if (!file.exists(f)) return("missing")
    out <- tryCatch(subset_file(f), error = function(e) conditionMessage(e))
    if (is.character(out)) return(out)
    saveRDS(out, part_of(f))   # NULL (no soundings) is saved too, so it's skipped next time
    "ok"
  }, mc.cores = n_workers)
  bad <- which(unlist(res) != "ok")
  for (i in bad) message("  ", basename(local[i]), ": ", res[[i]])

  done <- sum(file.exists(part_of(files)))
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("batch %d/%d: %d/%d files done, %.0f min elapsed",
                  b, length(batches), done, length(files), el))
}
unlink(dl_dir, recursive = TRUE)

# ---- Combine ----------------------------------------------------------------
missing <- files[!file.exists(part_of(files))]
if (length(missing)) {
  warning(length(missing), " files could not be processed; rerun to retry:\n",
          paste(basename(missing), collapse = "\n"))
}
parts <- Filter(Negate(is.null), lapply(part_of(files[!files %in% missing]), readRDS))
# v11.2r adds a few quality-flag variables and renames one, so columns differ
# between versions; missing ones are filled with NA
oco <- rbindlist(lapply(parts, as.data.table), fill = TRUE)
oco[, lite_build := sub(".*_(B[0-9]+[A-Za-z]*)_.*", "\\1", source_file)]
oco <- st_sf(oco[order(time_utc)], sf_column_name = "geometry")

st_write(oco, out_file, layer = "oco2_sif", delete_dsn = file.exists(out_file),
         quiet = TRUE)
cat("\n", nrow(oco), "soundings from", length(unique(oco$source_file)), "of",
    length(files), "files written to", out_file, "\n")
print(table(MeasurementMode = oco$Metadata_MeasurementMode))
print(range(oco$time_utc))

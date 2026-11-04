# Build sub_road_crop/sub_road_crop<year>.shp (majority CDL crop per
# sub_road_grid cell) from the 30 m CDL, for CDL years that don't have one yet.
# Existing files are left alone.
#
# Fields, as in the existing files:
#   FID       grid cell id (= sub_road_grid.shp FID)
#   crop      CDL code covering the largest area of the cell
#   maj_frac  fraction of the cell's area covered by that code
#
# Note: the earlier files (2014-2024) agree with this on `crop` for ~98% of
# cells, but their maj_frac is typically about half of the area fraction
# computed here, so maj_frac is not comparable between those years and
# years made with this script.

library(terra)
library(sf)
library(exactextractr)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
cdl_dir  <- file.path(proj_dir, "cdl_30m")
crop_dir <- file.path(proj_dir, "sub_road_crop")
grid     <- st_read(file.path(proj_dir, "sub_road_grid.shp"), quiet = TRUE)

cdl_years <- sub("^CDL_([0-9]{4})[.]tif$", "\\1",
                 list.files(cdl_dir, "^CDL_[0-9]{4}[.]tif$"))
for (y in cdl_years) {
  out <- file.path(crop_dir, paste0("sub_road_crop", y, ".shp"))
  if (file.exists(out)) next
  r <- rast(file.path(cdl_dir, paste0("CDL_", y, ".tif")))
  maj <- exact_extract(r, st_transform(grid, crs(r)), function(v, w) {
    ok <- !is.na(v) & v != 0
    if (!any(ok)) return(c(NA, NA))
    a <- tapply(w[ok], v[ok], sum)
    c(as.numeric(names(a)[which.max(a)]), max(a) / sum(w))
  }, progress = FALSE)
  crops <- st_sf(FID = grid$FID, crop = maj[1, ], maj_frac = maj[2, ],
                 geometry = st_geometry(grid))
  st_write(crops, out, quiet = TRUE)
  cat("Wrote", out, "\n")
}

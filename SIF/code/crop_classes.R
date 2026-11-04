# Crop classes from the 30 m USDA Cropland Data Layer, shared by
# downscale_sif_by_crop_cdl30.R and compare_sifnet.R.
#
# Before sourcing, set cdl_dir (folder with CDL_<year>.tif). Defines:
#   class_levels      class names, in order (index = class raster value)
#   code_to_class()   CDL codes -> class index
#   cdl               list of CDL rasters by year; cdl_crs their CRS
#   class_raster(d)   class raster (values 1..length(class_levels)) for date d
#
# Crop year: a CDL year labels the crop harvested that year, so from
# fall_start on, pixels that are winter wheat (or a double crop starting with
# it) in next year's CDL take next year's class; all others use the
# calendar year.

fall_start  <- "09-01"
wheat_codes <- c(24, 26, 225, 236)   # Winter Wheat and WinWht double crops

# Crop classes (CDL codes); codes not listed go to "Other". Same classes as
# downscale_sif_by_crop.R so the two can be compared.
crop_groups <- list(
  "Winter Wheat"             = 24,
  "Dbl Crop WinWht/Soybeans" = 26,
  "Soybeans"                 = 5,
  "Corn"                     = 1,
  "Sorghum"                  = 4,
  "Cotton"                   = 2,
  "Grassland/Pasture"        = 176,
  "Fallow/Idle Cropland"     = 61,
  "Developed"                = 121:124,
  "Trees/Wetlands"           = c(141, 142, 143, 190, 195)
)
class_levels <- c(names(crop_groups), "Other")
code_to_class <- function(codes) {
  out <- rep(length(class_levels), length(codes))      # "Other"
  for (k in seq_along(crop_groups)) out[codes %in% crop_groups[[k]]] <- k
  out[is.na(codes) | codes == 0] <- NA                 # 0 = no data
  out
}

cdl_files <- list.files(cdl_dir, "^CDL_[0-9]{4}[.]tif$", full.names = TRUE)
cdl <- setNames(lapply(cdl_files, terra::rast),
                sub("^CDL_([0-9]{4})[.]tif$", "\\1", basename(cdl_files)))
cdl_crs <- terra::crs(cdl[[1]])

class_raster <- function(d) {
  yr <- format(d, "%Y")
  if (!yr %in% names(cdl)) stop("No CDL for ", yr)
  codes <- cdl[[yr]]
  nxt <- as.character(as.integer(yr) + 1)
  if (format(d, "%m-%d") >= fall_start && nxt %in% names(cdl)) {
    codes <- terra::ifel(cdl[[nxt]] %in% wheat_codes, cdl[[nxt]], codes)
  }
  terra::app(codes, code_to_class)
}

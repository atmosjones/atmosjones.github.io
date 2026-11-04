# Add the dominant crop type of each sub_road grid cell (USDA Cropland Data
# Layer, one shapefile per year in sub_road_crop/) to the downscaled SIF in
# downscaled_sif757/downscaled_sif757.gpkg, matching each date to the crop
# map for its calendar year - except in the fall: a CDL year labels the crop
# harvested that year, so winter wheat sown in Sep-Oct appears in the NEXT
# year's map. From fall_start on, a cell that is winter wheat (or a double
# crop starting with it) in next year's map takes next year's crop.
#
# Adds to layer "cells":
#   crop_year   CDL year used
#   crop_code   CDL class code of the dominant crop
#   crop_name   CDL class name
#   crop_frac   fraction of the cell covered by that crop (maj_frac)
#
# Run after downscale_sif_inversion.R; rerunning the inversion rewrites the
# GeoPackage without these columns.

library(sf)
library(data.table)

proj_dir  <- path.expand("~/Documents/cropland_downscaling_paper")
crop_dir  <- file.path(proj_dir, "sub_road_crop")
sif_file  <- file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg")
sif_layer <- "cells"

fall_start  <- "09-01"                 # month-day from which fall planting applies
wheat_codes <- c(24, 26, 225, 236)     # Winter Wheat and WinWht double crops

# CDL class names for the codes in this area
cdl_names <- c(
  "1" = "Corn", "2" = "Cotton", "4" = "Sorghum", "5" = "Soybeans",
  "21" = "Barley", "24" = "Winter Wheat",
  "26" = "Dbl Crop WinWht/Soybeans", "28" = "Oats", "31" = "Canola",
  "36" = "Alfalfa", "37" = "Other Hay/Non Alfalfa", "42" = "Dry Beans",
  "57" = "Herbs", "59" = "Sod/Grass Seed", "61" = "Fallow/Idle Cropland",
  "111" = "Open Water", "121" = "Developed/Open Space",
  "122" = "Developed/Low Intensity", "123" = "Developed/Med Intensity",
  "124" = "Developed/High Intensity", "141" = "Deciduous Forest",
  "176" = "Grassland/Pasture", "190" = "Woody Wetlands",
  "195" = "Herbaceous Wetlands", "205" = "Triticale",
  "225" = "Dbl Crop WinWht/Corn", "236" = "Dbl Crop WinWht/Sorghum",
  "240" = "Dbl Crop Soybeans/Oats", "254" = "Dbl Crop Barley/Soybeans")

# ---- Crop maps, all years ---------------------------------------------------
crop_files <- sort(list.files(crop_dir, "^sub_road_crop[0-9]{4}[.]shp$", full.names = TRUE))
crops <- rbindlist(lapply(crop_files, function(f) {
  d <- st_drop_geometry(st_read(f, quiet = TRUE))
  data.table(grid_fid  = d$FID,
             crop_year = as.integer(sub("\\D*(\\d{4})[.]shp$", "\\1", basename(f))),
             crop_code = as.integer(d$crop),
             crop_frac = d$maj_frac)
}))
crops[, crop_name := unname(cdl_names[as.character(crop_code)])]
unknown <- crops[is.na(crop_name), unique(crop_code)]
if (length(unknown)) warning("No name for CDL codes: ", paste(unknown, collapse = ", "))
cat("Crop maps for", paste(range(crops$crop_year), collapse = "-"), "\n")

# ---- Join to the downscaled SIF ---------------------------------------------
sif <- st_read(sif_file, layer = sif_layer, quiet = TRUE)
sif <- sif[, setdiff(names(sif), c("crop_year", "crop_code", "crop_name", "crop_frac"))]
d <- as.Date(sif$date)
sif$crop_year <- as.integer(format(d, "%Y"))

# Fall dates: use next year's map where it shows winter wheat
wheat_next <- crops[crop_code %in% wheat_codes, .(grid_fid, crop_year = crop_year - 1L)]
is_fall <- format(d, "%m-%d") >= fall_start
to_next <- is_fall & paste(sif$grid_fid, sif$crop_year) %in%
  paste(wheat_next$grid_fid, wheat_next$crop_year)
sif$crop_year[to_next] <- sif$crop_year[to_next] + 1L
cat(sum(to_next), "fall cell-days use next year's map (winter wheat planted);",
    sum(is_fall & !to_next), "fall cell-days keep the calendar year\n")

n_before <- nrow(sif)
sif <- merge(sif, crops, by = c("grid_fid", "crop_year"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(sif) == n_before)
sif <- sif[order(sif$date, sif$cell_id), ]

missing_years <- setdiff(unique(sif$crop_year), unique(crops$crop_year))
if (length(missing_years)) {
  warning("No crop map for year(s) ", paste(missing_years, collapse = ", "),
          "; crop columns are NA there")
}

# Write to a temporary file first so the original survives a failed write
tmp <- sub("[.]gpkg$", "_tmp.gpkg", sif_file)
st_write(sif, tmp, layer = sif_layer, delete_dsn = TRUE, quiet = TRUE)
invisible(file.rename(tmp, sif_file))

cat(nrow(sif), "cell-days updated in", sif_file, "\n")
print(as.data.table(st_drop_geometry(sif))[, .(cell_days = .N,
                                              mean_sif757 = round(mean(sif757), 3)),
                                          keyby = crop_name][order(-cell_days)])

# Crop area and crop rotation in the study area (raster_extent.shp) from the
# 30 m USDA Cropland Data Layer, 2014-2024.
#
# Area: share of the study area in each crop class (classes from
# crop_classes.R) for every CDL year.
#
# Rotation rate: of the pixels that are cropland in two consecutive years,
# the fraction whose CDL crop code changed. Cropland = CDL crop codes
# (1-61, 66-80, 204-254; includes fallow, excludes grassland/pasture,
# hay is included as a crop). Codes are compared at full CDL detail, so e.g.
# winter wheat -> double crop wheat/soybeans counts as a change.
# CDL classification errors also show up as "changes", so this is an upper
# bound on true rotation. rotation_rate_interior repeats it for field
# interiors only (the pixel's 3x3 neighbourhood is one crop in both years),
# which removes most edge and speckle errors; true rotation is likely
# between the two.
#
# Outputs (in proj_dir):
#   crop_area_by_year.csv         year, crop class, area (km2), share of study area
#   crop_rotation_by_year.csv     year pair, cropland pixels, fraction changed
#   crop_transitions.csv          year pair, from crop, to crop, pixels (cropland)

library(terra)
library(sf)
library(data.table)

proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
cdl_dir    <- file.path(proj_dir, "cdl_30m")
extent_shp <- file.path(proj_dir, "raster_extent.shp")

source(file.path(proj_dir, "crop_classes.R"))   # cdl, class_levels, code_to_class()

cropland_codes <- c(1:61, 66:80, 204:254)

cdl_names <- c(
  "1" = "Corn", "2" = "Cotton", "4" = "Sorghum", "5" = "Soybeans", "21" = "Barley",
  "24" = "Winter Wheat", "26" = "Dbl Crop WinWht/Soybeans", "27" = "Rye", "28" = "Oats",
  "31" = "Canola", "36" = "Alfalfa", "37" = "Other Hay/Non Alfalfa", "42" = "Dry Beans",
  "44" = "Other Crops", "57" = "Herbs", "58" = "Clover/Wildflowers", "59" = "Sod/Grass Seed",
  "61" = "Fallow/Idle Cropland", "205" = "Triticale", "225" = "Dbl Crop WinWht/Corn",
  "236" = "Dbl Crop WinWht/Sorghum", "238" = "Dbl Crop WinWht/Cotton",
  "240" = "Dbl Crop Soybeans/Oats", "254" = "Dbl Crop Barley/Soybeans")
name_of <- function(code) {
  n <- cdl_names[as.character(code)]
  ifelse(is.na(n), paste("CDL", code), n)
}

# ---- CDL clipped to the study area ------------------------------------------
study <- project(vect(st_read(extent_shp, quiet = TRUE)), cdl_crs)
years <- sort(as.integer(names(cdl)))
codes <- rast(lapply(as.character(years), function(y) mask(crop(cdl[[y]], study), study)))
names(codes) <- years
pix_km2 <- prod(res(codes)) / 1e6

# ---- Area by crop class -----------------------------------------------------
area <- rbindlist(lapply(seq_along(years), function(k) {
  v <- values(codes[[k]])[, 1]
  v <- v[!is.na(v) & v != 0]
  cls <- factor(class_levels[code_to_class(v)], levels = class_levels)
  n <- table(cls)
  data.table(year = years[k], crop_class = names(n), pixels = as.integer(n),
             area_km2 = as.integer(n) * pix_km2, share = as.integer(n) / length(v))
}))
fwrite(area, file.path(proj_dir, "crop_area_by_year.csv"))

# ---- Rotation between consecutive years -------------------------------------
# 3x3 uniformity: the pixel's code equals the min and max of its neighbourhood
uniform <- function(r) {
  (focal(r, w = 3, fun = "min", na.rm = FALSE) == r) &
    (focal(r, w = 3, fun = "max", na.rm = FALSE) == r)
}

rotation <- list(); transitions <- list()
for (k in seq_along(years)[-1]) {
  a <- codes[[k - 1]]; b <- codes[[k]]
  ok <- values((a %in% cropland_codes) & (b %in% cropland_codes))[, 1] == 1
  interior <- ok & values(uniform(a) & uniform(b))[, 1] %in% 1
  va <- values(a)[, 1]; vb <- values(b)[, 1]
  changed <- va != vb
  rotation[[k]] <- data.table(from_year = years[k - 1], to_year = years[k],
                              cropland_km2 = sum(ok) * pix_km2,
                              rotation_rate = mean(changed[ok]),
                              interior_km2 = sum(interior) * pix_km2,
                              rotation_rate_interior = mean(changed[interior]))
  keep <- which(ok)
  va <- va[keep]; vb <- vb[keep]
  tr <- data.table(from = va, to = vb)[, .(pixels = .N), by = .(from, to)]
  tr[, `:=`(from_year = years[k - 1], to_year = years[k],
            from_crop = name_of(from), to_crop = name_of(to))]
  transitions[[k]] <- tr
}
rotation <- rbindlist(rotation)
transitions <- rbindlist(transitions)[order(from_year, -pixels)]
transitions[, share_of_cropland := pixels / sum(pixels), by = from_year]
fwrite(rotation, file.path(proj_dir, "crop_rotation_by_year.csv"))
fwrite(transitions, file.path(proj_dir, "crop_transitions.csv"))

# ---- Summary ----------------------------------------------------------------
cat("Rotation rate (fraction of cropland pixels changing crop):\n")
print(rotation[, .(years = paste(from_year, to_year, sep = "-"),
                   cropland_km2 = round(cropland_km2), rotation_rate = round(rotation_rate, 3),
                   interior_km2 = round(interior_km2),
                   rotation_rate_interior = round(rotation_rate_interior, 3))])
cat(sprintf("Mean %.3f (field interiors %.3f)\n", mean(rotation$rotation_rate),
            mean(rotation$rotation_rate_interior)))
cat("\nMost common changes, all years (share of cropland pixel-years):\n")
print(transitions[from != to, .(pixels = sum(pixels)), by = .(from_crop, to_crop)][
  , share := round(pixels / sum(transitions$pixels), 3)][order(-pixels)][1:8])

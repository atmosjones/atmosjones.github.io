# Figure 1 (alternative): the cropland mosaic and its patch sizes.
#
#   a  2024 Cropland Data Layer at 10 m (CroplandCROS export), grouped into the
#      paper's crop classes, with the section-road grid used for downscaling
#   b  Share of area in contiguous patches at least a given size, stacked by
#      class (CDL 2014-2023), with the resolution of the grid and of each SIF
#      product marked (from crop_patch_sizes.R)
#
# Panel a uses 2024 because it is the only 10 m CDL year; panel b pools the
# 30 m years 2014-2023 (2024 patch statistics aren't comparable, see
# crop_patch_sizes.R).
#
# Needs: crop_patch_sizes.R run first (patch_sizes/plot_patch_size_stacked.rds).
# Output: fig1_cropland_patches.png / .pdf, and cdl2024_classes_10m.tif (the classed
# map in panel a, NAD83 / Oklahoma North, with the figure's colors and class names).

library(sf)
library(terra)
library(ggplot2)
library(patchwork)

# ---- Settings ---------------------------------------------------------------
proj_dir  <- path.expand("~/Documents/cropland_downscaling_paper")
map_year  <- 2024
cros_file <- path.expand(sprintf("~/Documents/downscaling/CROS/%d_CROS/clipped.mrf", map_year))
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
stack_rds <- file.path(proj_dir, "patch_sizes", "plot_patch_size_stacked.rds")
cdl_dir   <- file.path(proj_dir, "cdl_30m")                      # needed by crop_classes.R
site_ll   <- c(lon = -97.486, lat = 36.604)                      # Lamont TCCON / ARM SGP CF

source(file.path(proj_dir, "crop_classes.R"))   # class_levels, code_to_class()
crop_class_colors <- c(                          # as in make_plots.R / crop_patch_sizes.R
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d",
  "Soybeans" = "#1b9e77", "Corn" = "#d95f02", "Sorghum" = "#e7298a",
  "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e",
  "Fallow/Idle Cropland" = "#bdb76b", "Trees/Wetlands" = "#1f4e3d",
  "Developed" = "grey55", "Other" = "grey80")

# ---- a. 10 m crop map + grid ------------------------------------------------
grid <- st_read(grid_file, quiet = TRUE)
crs_g <- st_crs(grid)
domain <- st_union(st_geometry(grid))
site <- st_transform(st_sfc(st_point(site_ll), crs = 4326), crs_g)

r <- rast(cros_file)
levels(r) <- NULL                                         # raw CDL codes
dom_r <- vect(st_transform(st_sf(geometry = domain), crs(r)))
r <- crop(r, dom_r, snap = "out")
r <- project(r, crs_g$wkt, method = "near", res = res(rast(cros_file))[1] * cos(36.6 * pi / 180))
r <- mask(r, vect(st_sf(geometry = domain)))
cls <- app(r, code_to_class)

# Save the classed map as a GeoTIFF with the figure's colors and class names
# embedded (opens in QGIS / ArcGIS with the same palette as the figure)
out_tif <- file.path(proj_dir, sprintf("cdl%d_classes_%dm.tif", map_year, round(res(rast(cros_file))[1])))
cls_out <- cls
names(cls_out) <- "crop_class"
levels(cls_out) <- data.frame(value = seq_along(class_levels), crop_class = class_levels)
rgb <- col2rgb(crop_class_colors[class_levels])
coltab(cls_out) <- data.frame(value = seq_along(class_levels), red = rgb[1, ], green = rgb[2, ],
                              blue = rgb[3, ], alpha = 255)
writeRaster(cls_out, out_tif, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))
cat(sprintf("%d CDL, %.1f m pixels, %d x %d\n", map_year, res(cls)[1], nrow(cls), ncol(cls)))

# draw the raster as an image (millions of pixels: much lighter than geom_raster)
v <- as.matrix(cls, wide = TRUE)
img <- matrix(c(crop_class_colors[class_levels], NA)[ifelse(is.na(v), length(class_levels) + 1, v)],
              nrow = nrow(v))
img[is.na(img)] <- "#FFFFFF00"
e <- as.vector(ext(cls))
bb <- st_bbox(domain)
present <- class_levels[sort(unique(as.vector(v[!is.na(v)])))]
key <- data.frame(x = bb[["xmin"]], y = bb[["ymin"]],
                  crop_class = factor(present, levels = names(crop_class_colors)))
km <- 5; x0 <- bb[["xmin"]] + 600; y0 <- bb[["ymin"]] - 1300; h <- 300
pa <- ggplot() +
  annotation_raster(img, xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], interpolate = FALSE) +
  geom_tile(data = key, aes(x = x, y = y, fill = crop_class), width = 0, height = 0) +   # legend only
  geom_sf(data = grid, fill = NA, color = "white", linewidth = 0.12) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = site, shape = 24, fill = "black", color = "white", size = 2) +
  annotate("rect", xmin = x0, xmax = x0 + km * 500, ymin = y0, ymax = y0 + h, fill = "black") +
  annotate("rect", xmin = x0 + km * 500, xmax = x0 + km * 1000, ymin = y0, ymax = y0 + h,
           fill = "white", color = "black", linewidth = 0.3) +
  annotate("text", x = c(x0, x0 + km * 1000), y = y0 + h + 550, label = c("0", paste(km, "km")), size = 2.6) +
  scale_fill_manual(values = crop_class_colors, drop = TRUE, name = NULL,
                    guide = guide_legend(ncol = 3, byrow = TRUE,
                                         override.aes = list(width = 1, height = 1))) +
  coord_sf(crs = crs_g, xlim = c(bb[["xmin"]] - 800, bb[["xmax"]] + 800),
           ylim = c(bb[["ymin"]] - 1900, bb[["ymax"]] + 800), expand = FALSE, datum = NA) +
  labs(title = sprintf("%d Cropland Data Layer (10 m) and the %s-cell section-road grid",
                       map_year, format(nrow(grid), big.mark = ","))) +
  theme_void(base_size = 9) +
  theme(panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.38, "cm"), legend.text = element_text(size = 7.5),
        legend.position = "bottom", legend.justification = "left", plot.title.position = "plot",
        plot.margin = margin(2, 4, 2, 2))

# Inset: where the study area is in Oklahoma (placed beside the legend, so it hides no map)
states <- st_read(file.path(proj_dir, "context_data", "cb_2023_us_state_20m.shp"), quiet = TRUE)
ok_aea <- st_transform(states[states$STUSPS == "OK", ], 5070)
nb_aea <- st_transform(states[states$STUSPS %in% c("TX", "KS", "CO", "NM", "MO", "AR"), ], 5070)
dom_aea <- st_transform(st_sf(geometry = domain), 5070)
okb <- st_bbox(ok_aea)
inset <- ggplot() +
  geom_sf(data = nb_aea, fill = "grey96", color = "grey75", linewidth = 0.15) +
  geom_sf(data = ok_aea, fill = "grey85", color = "grey30", linewidth = 0.3) +
  geom_sf(data = st_buffer(st_centroid(dom_aea), 18000), fill = NA, color = "#c0392b", linewidth = 0.7) +
  geom_sf(data = dom_aea, fill = "#c0392b", color = "#c0392b", linewidth = 0.3) +
  annotate("text", x = (okb[["xmin"]] + okb[["xmax"]]) / 2 + 60000, y = okb[["ymin"]] + 90000,
           label = "Oklahoma", size = 2.4, fontface = "italic", color = "grey25") +
  coord_sf(crs = 5070, xlim = c(okb[["xmin"]] - 15000, okb[["xmax"]] + 15000),
           ylim = c(okb[["ymin"]] - 15000, okb[["ymax"]] + 15000), expand = FALSE, datum = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3))
pa <- pa + inset_element(inset, left = 0.755, bottom = 0.0, right = 1.0, top = 0.165,
                         align_to = "full", ignore_tag = TRUE)

# ---- b. stacked patch-size distribution -------------------------------------
pb <- readRDS(stack_rds) + guides(fill = "none") +          # panel a's legend covers both
  labs(title = "Area in contiguous patches at least a given size (CDL 2014-2023)") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9), axis.text = element_text(size = 8))

# ---- Assemble ---------------------------------------------------------------
# wrap_elements: panel a is drawn as one self-contained graphic (title, map, legend), so
# it isn't shrunk to match b's plot area, which b's tall top-axis labels push down
fig <- (wrap_elements(full = pa) | pb) + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
saveRDS(list(a = pa, b = pb), file.path(proj_dir, "fig1_panels_patches.rds"))   # for make_fig1_combined.R
ggsave(file.path(proj_dir, "fig1_cropland_patches.png"), fig, width = 13, height = 6.6, dpi = 300, bg = "white")
ggsave(file.path(proj_dir, "fig1_cropland_patches.pdf"), fig, width = 13, height = 6.6, bg = "white")

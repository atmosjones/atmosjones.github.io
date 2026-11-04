# Distribution of contiguous crop areas (patches) in the study area, from the
# CroplandCROS exports of the USDA Cropland Data Layer, 2014-2024.
#
# A patch is a set of same-class pixels joined edge-to-edge (rook / 4-neighbour
# connectivity; set `directions` to 8 to also join diagonal neighbours). CDL
# codes are grouped into the paper's crop classes (crop_classes.R), so e.g. all
# forest and wetland codes form one "Trees/Wetlands" class.
#
# Resolution: the CDL is 10 m from 2024 on (30 m before); those years are
# aggregated to 30 m (most common class) so patch topology compares across years.
#
# Pixel areas: the CROS exports are 30 m in Web Mercator (EPSG:3857), which at
# 36.6 N is really ~24 m on the ground, so areas come from terra::cellSize()
# (true geodesic area), not 900 m2 per pixel.
#
# Patches touching the edge of raster_extent.shp are cut off by it; they are
# flagged (touches_edge) and the summary reports their share of each class.
#
# The patch-count distribution is dominated by CDL speckle (1-2 pixel patches),
# so the summary centres on AREA-WEIGHTED statistics: the size of the patch an
# average pixel of the class belongs to. Reference scales (ref_scales): the
# section-road grid cell and the resolution of each SIF product compared.
#
# Output (in out_dir):
#   patches.csv               one row per patch: year, class, area, edge flag
#   patch_size_summary.csv    one row per year x class
#   patch_size_by_class.csv   pooled over all years
#   plot_patch_size_distribution.png
#   plot_patch_size_stacked.png (+ .rds, the ggplot object, for make_fig1_patches.R)

library(terra)
library(sf)
library(data.table)
library(ggplot2)

# ---- Settings ---------------------------------------------------------------
proj_dir   <- path.expand("~/Documents/cropland_downscaling_paper")
cros_dir   <- path.expand("~/Documents/downscaling/CROS")        # <year>_CROS/clipped.mrf
extent_shp <- file.path(proj_dir, "raster_extent.shp")
grid_file  <- file.path(proj_dir, "sub_road_grid.shp")
cdl_dir    <- file.path(proj_dir, "cdl_30m")                      # needed by crop_classes.R
out_dir    <- file.path(proj_dir, "patch_sizes")
directions <- 4
target_res <- 30           # CDL 2024 onward is 10 m; aggregate (modal) to 30 m so years compare
pool_years <- 2014:2023    # pooled stats + figure. The 10 m CDL (2024-) still gives much larger
                           # patches after aggregation (a different classification), so it is kept
                           # in the per-year table only
min_patch_ha <- 1          # patch counts in the summary also given for patches >= this size

dir.create(out_dir, showWarnings = FALSE)
source(file.path(proj_dir, "crop_classes.R"))   # class_levels, code_to_class()

# ---- Reference scales (km2) -------------------------------------------------
# Lat/lon grids converted at the site latitude (36.6 N). OCO-2 is the nominal
# 1.29 x 2.25 km footprint (target-mode footprints here are smaller, median
# ~0.5 km2, because pointing at the site cancels most of the along-track smear).
grid <- st_read(grid_file, quiet = TRUE)
cell_km2 <- median(as.numeric(st_area(grid))) / 1e6
deg2_km2 <- function(d, lat = 36.6) (d * 111.32 * cos(lat * pi / 180)) * (d * 110.95)
ref_scales <- c(FLEX =.3*.3, sifnet = deg2_km2(0.005), grid_cell = cell_km2, oco2 = 1.29 * 2.25,
                tropomi = 3.5 * 5.5, gosif = deg2_km2(0.05))
# plotmath labels (degree / times signs render in any locale)
# TROPOMI and GOSIF are too close on a log axis for two labels: one shared label,
# both ticks/lines kept
ref_labels <- c("FLEX", "SIFNet~(0.005*degree)", "Grid~cell", "OCO*'-'*2~(1.3%*%2.25~km)",
                "TROPOMI~(3.5%*%5.5~km)~'/'~GOSIF~(0.05*degree)", "''")
print(round(ref_scales, 3))

# ---- Patches, one year at a time --------------------------------------------
years <- sort(as.integer(sub("_CROS$", "", list.files(cros_dir, "^[0-9]{4}_CROS$"))))
extent <- vect(st_read(extent_shp, quiet = TRUE))

patch_year <- function(yr) {
  r <- rast(file.path(cros_dir, paste0(yr, "_CROS"), "clipped.mrf"))
  levels(r) <- NULL                                  # raw CDL codes, not the category labels
  f <- round(target_res / res(r)[1])
  if (f > 1) r <- aggregate(r, f, fun = "modal", na.rm = TRUE)
  ext_r <- project(extent, crs(r))
  r <- mask(crop(r, ext_r), ext_r)
  cls <- app(r, code_to_class)
  area <- cellSize(cls, unit = "km")                 # true area of each pixel
  # pixels inside the extent with a neighbour outside it (or on the raster edge)
  inside <- !is.na(r)
  nb_out <- focal(inside, w = matrix(c(NA, 1, NA, 1, NA, 1, NA, 1, NA), 3), fun = "min",
                  fillvalue = 0)                     # NA weights = ignored (diagonals, centre)
  edge <- inside & (nb_out == 0)

  rbindlist(lapply(seq_along(class_levels), function(k) {
    x <- ifel(cls == k, 1, NA)
    if (is.na(global(x, "max", na.rm = TRUE)[[1]])) return(NULL)
    p <- patches(x, directions = directions, zeroAsNA = TRUE)
    a <- zonal(c(area, edge), p, fun = c("sum"), na.rm = TRUE)
    n <- zonal(x, p, fun = "sum", na.rm = TRUE)
    data.table(year = yr, crop_class = class_levels[k], patch = a[[1]],
               n_px = n[[2]], area_km2 = a[[2]], touches_edge = a[[3]] > 0)
  }))
}

patches_dt <- rbindlist(lapply(years, function(yr) {
  t0 <- Sys.time(); out <- patch_year(yr)
  cat(yr, ":", nrow(out), "patches,", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  out
}))
fwrite(patches_dt, file.path(out_dir, "patches.csv"))

# ---- Summaries --------------------------------------------------------------
# area-weighted quantile: the patch size below which a given share of the
# class's AREA lies
wq <- function(a, p) { o <- order(a); a <- a[o]; cw <- cumsum(a) / sum(a); a[which(cw >= p)[1]] }
summarise <- function(d) d[, .(
  class_area_km2   = sum(area_km2),
  n_patches        = .N,
  n_patches_ge_min = sum(area_km2 * 100 >= min_patch_ha),
  median_patch_km2 = median(area_km2),
  aw_mean_km2      = sum(area_km2^2) / sum(area_km2),       # area-weighted mean patch size
  aw_p25_km2       = wq(area_km2, 0.25),
  aw_median_km2    = wq(area_km2, 0.50),
  aw_p75_km2       = wq(area_km2, 0.75),
  largest_km2      = max(area_km2),
  share_gt_sifnet    = sum(area_km2[area_km2 > ref_scales[[1]]]) / sum(area_km2),
  share_gt_grid_cell = sum(area_km2[area_km2 > ref_scales[[2]]]) / sum(area_km2),
  share_gt_oco2      = sum(area_km2[area_km2 > ref_scales[[3]]]) / sum(area_km2),
  share_gt_tropomi   = sum(area_km2[area_km2 > ref_scales[[4]]]) / sum(area_km2),
  share_gt_gosif     = sum(area_km2[area_km2 > ref_scales[[5]]]) / sum(area_km2),
  share_in_edge_patches = sum(area_km2[touches_edge]) / sum(area_km2))]
by_year <- patches_dt[, summarise(.SD), by = .(year, crop_class)]
pooled  <- patches_dt[year %in% pool_years, summarise(.SD), by = crop_class]
setorder(pooled, -class_area_km2)
fwrite(by_year, file.path(out_dir, "patch_size_summary.csv"))
fwrite(pooled,  file.path(out_dir, "patch_size_by_class.csv"))

# ---- Figure: share of class area in patches at least a given size -------------
# (area-weighted exceedance curve: read off at the footprint / grid-cell size,
# it is the share of the class that sits in patches the downscaling can resolve)
main <- pooled[!crop_class %in% c("Developed", "Other") & class_area_km2 > 5, crop_class]
cd <- patches_dt[year %in% pool_years & crop_class %in% main][order(-area_km2)]
cd[, share_ge := cumsum(area_km2) / sum(area_km2), by = crop_class]   # share in patches >= this size
cd[, crop_class := factor(crop_class, levels = main)]
ggplot(cd, aes(x = area_km2, y = share_ge, color = crop_class)) +
  geom_step(linewidth = 0.6, direction = "vh") +
  geom_vline(xintercept = ref_scales, color = "grey35", linewidth = 0.35,
             linetype = c(1, 1, 1, 1, 1, 1)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = c("0.001", "0.01", "0.1", "1", "10", "100"),
                limits = c(NA, 40),
                sec.axis = dup_axis(breaks = unname(ref_scales), labels = parse(text = ref_labels),
                                    name = NULL)) +
  scale_color_manual(values = c(
    "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
    "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3",
    "Grassland/Pasture" = "#66a61e", "Fallow/Idle Cropland" = "#bdb76b",
    "Trees/Wetlands" = "#1f4e3d")) +
  labs(x = expression(Contiguous~patch~area~(km^2)), y = "Share of class area in patches this size or larger",
       color = NULL, title = sprintf("Contiguous crop areas, CDL %d-%d (%d-neighbour patches)",
                                     min(pool_years), max(pool_years), directions)) +
  theme_bw() +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 7.5),
        axis.ticks.x.top = element_line(color = "grey35"), plot.margin = margin(5.5, 75, 5.5, 5.5))
ggsave(file.path(out_dir, "plot_patch_size_distribution.png"), g, width = 8, height = 5.4, dpi = 300)

# ---- Figure: the same, stacked, as a share of all the classes shown ---------------
# Each band = that class's area in patches >= x, divided by the total area of
# the classes shown; at the left edge the bands are the classes' area shares.
x_min <- 4e-4                              # just below a single 30 m pixel (~0.00058 km2 here), so
xs <- 10^seq(log10(x_min), log10(40), length.out = 400)   # the plot starts at 100% and shows the 1-px step
total_area <- sum(pooled[crop_class %in% main, class_area_km2])
st <- cd[, {
  a <- sort(area_km2, decreasing = TRUE)
  ge <- vapply(xs, function(x) sum(a[a >= x]), numeric(1))
  .(x = xs, share = ge / total_area)
}, by = crop_class]
st[, crop_class := factor(crop_class, levels = rev(main))]     # largest class at the bottom
ggplot(st, aes(x = x, y = share, fill = crop_class)) +
  geom_area(position = "stack", outline.type = "full", color = "white", linewidth = 0.1) +
  #geom_vline(xintercept = ref_scales, color = "grey15", linewidth = 0.35,
  #           linetype = c(1, 1, 1, 1, 1, 1)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = c("0.001", "0.01", "0.1", "1", "10", "100"),
                limits = c(x_min, 40), expand = c(0, 0),
                sec.axis = dup_axis(breaks = unname(ref_scales), labels = parse(text = ref_labels),
                                    name = NULL)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +   # no hard limit: rounding can put the stack a hair over 1
  scale_fill_manual(values = c(
    "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
    "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3",
    "Grassland/Pasture" = "#66a61e", "Fallow/Idle Cropland" = "#bdb76b",
    "Trees/Wetlands" = "#1f4e3d"), breaks = levels(st$crop_class)) +   # first level is stacked on top
  labs(x = expression(Contiguous~patch~area~(km^2)),
       y = "Share of area in patches this size or larger", fill = NULL,
       title = sprintf("Contiguous crop areas by class, CDL %d-%d (%d-neighbour patches)",
                       min(pool_years), max(pool_years), directions)) +
  theme_classic()
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 7.5),
        axis.ticks.x.top = element_line(color = "grey35"), panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 75, 5.5, 5.5))
ggsave(file.path(out_dir, "plot_patch_size_stacked.png"), gs, width = 8, height = 5.4, dpi = 300)
saveRDS(gs, file.path(out_dir, "plot_patch_size_stacked.rds"))   # reused as Figure 1b (make_fig1_patches.R)
at_ref <- st[, .(x, share)][, .(total = sum(share)), by = x]
cat("\nStacked plot: share of the classes shown in patches >= each reference scale\n")
print(round(vapply(ref_scales, function(r) at_ref[which.min(abs(log10(x) - log10(r))), total], numeric(1)), 3))

cat("\nPooled over", length(pool_years), "years (area-weighted patch sizes, km2; class area = mean per year):\n")
print(pooled[, .(crop_class, class_area_km2 = class_area_km2 / length(pool_years),
                 aw_median_km2, largest_km2, share_gt_sifnet, share_gt_grid_cell, share_gt_oco2,
                 share_gt_tropomi, share_gt_gosif, share_edge = share_in_edge_patches)],
      digits = 2)

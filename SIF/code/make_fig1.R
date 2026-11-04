# Figure 1: the study domain and how overlapping OCO-2 target-mode soundings
# are inverted onto the section-road grid (downscale_sif_inversion.R).
#
#   a  Location: north-central Oklahoma, the four counties, the domain and the
#      Lamont TCCON / ARM SGP Central Facility that OCO-2 targets
#   b  Domain: 30 m CDL crop classes for the example day's crop year, with the
#      1,058 grid cells bounded by section roads
#   c  One target day: every footprint (grey outlines, showing the overlap)
#      and its measured SIF757 (points at footprint centres)
#   d  One footprint over the grid: the fraction of it in each cell = its row
#      of the forward operator A in y = A x + e (the manuscript's notation)
#   e  Downscaled SIF757 per grid cell that day (posterior mean)
#   f  Averaging kernel per cell (1 = set by the observations, 0 = prior only)
#   g  Overlapping soundings per cell (saved for Figure 2 only)
#   h  Effective resolution per cell (saved for Figure 2 only)
#
# Output: fig1_domain_inversion.png / .pdf in the project folder (the earlier QGIS fig1.* is left alone).

library(sf)
library(terra)
library(data.table)
library(ggplot2)
library(patchwork)

# ---- Settings ---------------------------------------------------------------
proj_dir  <- path.expand("~/Documents/cropland_downscaling_paper")
day       <- as.Date("2019-04-19")      # example target day (high DOFS, spring wheat contrast)
oco_file  <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
cell_file <- file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg")
states_shp <- file.path(proj_dir, "context_data", "cb_2023_us_state_20m.shp")
county_file <- file.path(proj_dir, "lamont_counties.gpkg")
cdl_dir   <- file.path(proj_dir, "cdl_30m")
site_ll   <- c(lon = -97.486, lat = 36.604)   # Lamont TCCON / ARM SGP Central Facility
target_mode <- 2; quality_flags <- c(0, 1); min_coverage <- 0.5   # as in downscale_sif_inversion.R

source(file.path(proj_dir, "crop_classes.R"))   # class_levels, class_raster()
crop_class_colors <- c(                          # as in make_plots.R
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d",
  "Soybeans" = "#1b9e77", "Corn" = "#d95f02", "Sorghum" = "#e7298a",
  "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e",
  "Fallow/Idle Cropland" = "#bdb76b", "Trees/Wetlands" = "#1f4e3d",
  "Developed" = "grey55", "Other" = "grey80")
sif_label <- expression(atop(SIF[757], (mW~m^-2~sr^-1~nm^-1)))

# ---- Data -------------------------------------------------------------------
grid <- st_read(grid_file, quiet = TRUE)
grid$cell_id <- seq_len(nrow(grid))
crs_g <- st_crs(grid)                       # NAD83 / Oklahoma North (m)
domain <- st_union(st_geometry(grid))
site <- st_transform(st_sfc(st_point(site_ll), crs = 4326), crs_g)

cells <- st_read(cell_file, quiet = TRUE,
                 query = sprintf("SELECT * FROM cells WHERE date = '%s'", day))
if (!nrow(cells)) stop("no downscaled cells for ", day)

fp <- st_read(oco_file, quiet = TRUE,
              query = paste("SELECT * FROM oco2_sif WHERE Metadata_MeasurementMode =", target_mode))
# Normalise SIF (and its uncertainty) to a nadir view: vza_correction.R
source(file.path(proj_dir, "vza_correction.R"))
fp <- correct_vza(fp)
fp <- fp[fp$Quality_Flag %in% quality_flags & !is.na(fp$Science_SIF_757nm) &
         !is.na(fp$Science_SIF_Uncertainty_757nm) & as.Date(fp$time_utc) == day, ]
fp <- st_transform(fp, crs_g)
st_geometry(fp) <- st_convex_hull(st_geometry(fp))
fp$fid_ <- seq_len(nrow(fp))
ins <- suppressWarnings(st_intersection(fp[, "fid_"], st_sf(geometry = domain)))
in_area <- tapply(as.numeric(st_area(ins)), ins$fid_, sum)
cov <- in_area[as.character(fp$fid_)] / as.numeric(st_area(fp))
fp <- fp[!is.na(cov) & cov >= min_coverage, ]
fp <- fp[order(fp$Science_SIF_757nm), ]
cat(nrow(fp), "footprints on", format(day), "\n")

# shared SIF colour range for c and e
sif_lim <- quantile(c(fp$Science_SIF_757nm, cells$sif757), c(0.02, 0.98))
sif_scale <- function(aes_name) {
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = sif_lim, oob = scales::squish,
                       name = sif_label, aesthetics = aes_name)
}

# ---- Helpers ----------------------------------------------------------------
bb <- st_bbox(domain)
scalebar <- function(km = 5, x0 = bb$xmin + 800, y0 = bb$ymin - 1500, h = 350) {
  m <- km * 1000
  list(annotate("rect", xmin = x0, xmax = x0 + m / 2, ymin = y0, ymax = y0 + h, fill = "black"),
       annotate("rect", xmin = x0 + m / 2, xmax = x0 + m, ymin = y0, ymax = y0 + h,
                fill = "white", color = "black", linewidth = 0.3),
       annotate("text", x = c(x0, x0 + m), y = y0 + h + 600, size = 2.5,
                label = c("0", paste(km, "km"))))
}
map_theme <- theme_void(base_size = 9) +
  theme(panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.4),
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 9, face = "bold"))
dom_lims <- coord_sf(crs = crs_g, xlim = c(bb$xmin - 1200, bb$xmax + 1200),
                     ylim = c(bb$ymin - 2200, bb$ymax + 1200), expand = FALSE, datum = NA)

# ---- a. Location ------------------------------------------------------------
states <- st_read(states_shp, quiet = TRUE)
states <- states[!states$STUSPS %in% c("AK", "HI", "PR"), ]
counties <- st_read(county_file, quiet = TRUE)
ok <- st_transform(states, 5070)
ok_bb <- st_bbox(st_transform(states[states$STUSPS == "OK", ], 5070))
dom_5070 <- st_transform(st_as_sfc(st_bbox(domain)), 5070)
pa <- ggplot() +
  geom_sf(data = ok, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = ok[ok$STUSPS == "OK", ], fill = "grey88", color = "grey30", linewidth = 0.4) +
  geom_sf(data = st_transform(counties, 5070), fill = NA, color = "grey35", linewidth = 0.25) +
  geom_sf(data = dom_5070, fill = NA, color = "#c0392b", linewidth = 0.6) +
  geom_sf(data = st_transform(site, 5070), shape = 24, fill = "black", size = 1.6) +
  annotate("text", x = (ok_bb$xmin + ok_bb$xmax) / 2 + 60000, y = (ok_bb$ymin + ok_bb$ymax) / 2 - 40000,
           label = "Oklahoma", size = 3.2, fontface = "italic", color = "grey30") +
  coord_sf(crs = 5070, xlim = c(ok_bb$xmin - 20000, ok_bb$xmax + 20000),
           ylim = c(ok_bb$ymin - 20000, ok_bb$ymax + 60000), expand = FALSE) +
  labs(title = "Location: Grant, Kay, Garfield and Noble counties") +
  theme_bw(base_size = 8) +
  theme(axis.title = element_blank(), panel.grid = element_line(color = "grey92", linewidth = 0.2),
        plot.title = element_text(size = 9, face = "bold"))
conus <- ggplot() +
  geom_sf(data = ok, fill = "grey90", color = "grey60", linewidth = 0.1) +
  geom_sf(data = ok[ok$STUSPS == "OK", ], fill = "grey45", color = NA) +
  coord_sf(crs = 5070, datum = NA) + theme_void() +
  theme(panel.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3))
pa <- pa + inset_element(conus, left = 0.0, bottom = 0.0, right = 0.36, top = 0.34, align_to = "panel",
                         ignore_tag = TRUE)

# ---- b. Domain: crop classes + grid -----------------------------------------
cls <- class_raster(day)
dom_v <- vect(st_transform(st_sf(geometry = domain), crs(cls)))
cls <- mask(crop(cls, dom_v), dom_v)
cls <- project(cls, crs_g$wkt, method = "near", res = 30)
cdf <- as.data.frame(cls, xy = TRUE, na.rm = TRUE)
names(cdf)[3] <- "cls"
cdf$crop_class <- factor(class_levels[cdf$cls], levels = names(crop_class_colors))
pb <- ggplot() +
  geom_raster(data = cdf, aes(x = x, y = y, fill = crop_class)) +
  geom_sf(data = grid, fill = NA, color = "white", linewidth = 0.15) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = site, shape = 24, fill = "black", color = "white", size = 2) +
  scale_fill_manual(values = crop_class_colors, drop = TRUE, name = NULL) +
  scalebar() + dom_lims +
  labs(title = sprintf("Domain: %d section-road cells, %s CDL", nrow(grid),
                       format(day, "%Y"))) +
  map_theme + theme(legend.key.height = unit(0.32, "cm"), legend.text = element_text(size = 6.5))

# ---- d. one footprint's row of K --------------------------------------------
# a footprint near the domain centre that is split over several cells
ctr <- st_centroid(domain)
near <- fp[as.numeric(st_distance(st_centroid(st_geometry(fp)), ctr)) < 3000, ]
pieces_of <- function(f) {
  p <- suppressWarnings(st_intersection(st_sf(geometry = st_geometry(f)), grid[, "cell_id"]))
  p$frac <- as.numeric(st_area(p)) / sum(as.numeric(st_area(p)))
  p
}
score <- vapply(seq_len(nrow(near)), function(i) {
  p <- pieces_of(near[i, ]); sum(p$frac > 0.08) - 0.01 * sum(p$frac <= 0.08)
}, numeric(1))
one <- near[which.max(score), ]
pc <- pieces_of(one)
zb <- st_bbox(st_buffer(st_geometry(one), 900))
pd <- ggplot() +
  geom_sf(data = pc, aes(fill = frac), color = NA) +
  geom_sf(data = grid, fill = NA, color = "grey45", linewidth = 0.35) +
  geom_sf(data = one, fill = NA, color = "black", linewidth = 0.8) +
  geom_sf_text(data = pc[pc$frac >= 0.03, ], aes(label = sprintf("%.2f", frac)), size = 2.6) +
  scale_fill_gradient(low = "#eef4fa", high = "#4a7fb5", limits = c(0, NA), guide = "none") +
  annotate("label", x = zb$xmin + 80, y = zb$ymin + 80, hjust = 0, vjust = 0, size = 2.8,
           label.size = 0, fill = alpha("white", 0.85), parse = TRUE,
           label = "y[i] == sum(A[ij]~x[j], j) + epsilon[i]") +
  annotate("label", x = zb$xmin + 80, y = zb$ymax - 80, hjust = 0, vjust = 1, size = 2.4,
           label.size = 0, fill = alpha("white", 0.85),
           label = sprintf("footprint i: y = %.2f\nA[ij] = fraction of i in cell j", one$Science_SIF_757nm)) +
  coord_sf(crs = crs_g, xlim = c(zb$xmin, zb$xmax), ylim = c(zb$ymin, zb$ymax), expand = FALSE, datum = NA) +
  labs(title = "One footprint = one row of A") +
  map_theme

# ---- c. footprints of the day -----------------------------------------------
fpc <- st_centroid(st_geometry(fp))
pc_dat <- st_sf(sif = fp$Science_SIF_757nm, geometry = fpc)
pcp <- ggplot() +
  geom_sf(data = st_sf(geometry = domain), fill = "grey97", color = NA) +
  geom_sf(data = fp, fill = NA, color = alpha("grey30", 0.06), linewidth = 0.15) +
  geom_sf(data = pc_dat, aes(color = sif), size = 0.55, stroke = 0) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = st_as_sfc(zb), fill = NA, color = "#c0392b", linewidth = 0.5) +
  geom_sf(data = one, fill = NA, color = "black", linewidth = 0.5) +
  sif_scale("colour") + guides(colour = "none") +
  scalebar() + dom_lims +
  labs(title = sprintf("%s: %s OCO-2 footprints", format(day, "%d %b %Y"),
                       format(nrow(fp), big.mark = ","))) +
  map_theme

# ---- e. downscaled SIF; f. averaging kernel ---------------------------------
dofs <- sum(cells$avg_kernel)
pe <- ggplot() +
  geom_sf(data = cells, aes(fill = sif757), color = NA) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  sif_scale("fill") + scalebar() + dom_lims +
  labs(title = "Inverted SIF per cell (posterior mean)") +
  map_theme
# averaging kernel: diverging about 0.5, the point where data and prior contribute equally
pf <- ggplot() +
  geom_sf(data = cells, aes(fill = avg_kernel), color = NA) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       name = "Averaging\nkernel") +
  scalebar() + dom_lims +
  labs(title = sprintf("Averaging kernel (DOFS %.0f; %d cells > 0.5)",
                       dofs, sum(cells$avg_kernel > 0.5))) +
  map_theme

# ---- g. overlapping soundings per cell; h. effective resolution ---------------
pg <- ggplot() +
  geom_sf(data = cells, aes(fill = n_footprints), color = NA) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis_c(option = "rocket", direction = -1, trans = "sqrt", name = "Soundings\nper cell") +
  scalebar() + dom_lims +
  labs(title = sprintf("Overlapping soundings per cell (max %d)", max(cells$n_footprints))) +
  map_theme
# width of each averaging-kernel row (participation area, downscale_sif_inversion.R);
# grey = estimate mostly prior (row sum < 0.5); outline = finer than a nominal
# OCO-2 footprint, away from the grid edge where rows are truncated (effective_resolution.R)
nominal_fp <- 2.9
edge_line <- st_cast(domain, "MULTILINESTRING")
cells$interior <- as.numeric(st_distance(st_centroid(st_geometry(cells)), edge_line)) > 1500
cells$res_plot <- ifelse(cells$ak_row_sum > 0.5, cells$res_area_pr_km2, NA)
fine <- st_union(cells[!is.na(cells$res_plot) & cells$res_plot < nominal_fp & cells$interior, ])
ph <- ggplot() +
  geom_sf(data = cells, aes(fill = res_plot), color = NA) +
  geom_sf(data = st_sf(geometry = fine), fill = NA, color = "black", linewidth = 0.45) +
  geom_sf(data = st_sf(geometry = domain), fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis_c(option = "mako", trans = "log10", limits = c(0.64, 12), oob = scales::squish,
                       breaks = c(0.64, 1, 2.9, 10), labels = c("0.64 (cell)", "1", "2.9 (OCO-2)", "10"),
                       na.value = "grey85", name = expression(km^2)) +
  scalebar() + dom_lims +
  labs(title = bquote(bold("Effective resolution (median" ~ .(sprintf("%.1f", median(cells$res_area_pr_km2[cells$avg_kernel > 0.5]))) ~
                                     km^2 ~ "where AK > 0.5)"))) +
  map_theme

# ---- Assemble ---------------------------------------------------------------
fig <- (pa + pb) / (pcp + pd) / (pe + pf) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
saveRDS(list(c = pcp, d = pd, e = pe, f = pf, g = pg, h = ph), file.path(proj_dir, "fig1_panels_inversion.rds"))   # for make_fig2.R
ggsave(file.path(proj_dir, "fig1_domain_inversion.png"), fig, width = 8, height = 9.2, dpi = 300, bg = "white")
ggsave(file.path(proj_dir, "fig1_domain_inversion.pdf"), fig, width = 8, height = 9.2, bg = "white")
cat("example footprint", one$Metadata_SoundingId, "split over", nrow(pc), "cells;",
    "DOFS", round(dofs), "\n")

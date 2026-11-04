# Figure 1 for the manuscript, in the layout of fig1_v3.pptx, drawn entirely in R
# so it is sharp at any size (vector except the 10 m crop raster, which is
# embedded at full resolution).
#
#   a  top:    Oklahoma, its counties and neighbouring states, with the study area
#      bottom: 2024 Cropland Data Layer (10 m, cdl2024_classes_10m.tif from
#              make_fig1_patches.R), the section-road grid, county lines, the
#              Lamont TCCON site (star) and a 25 km scale bar
#   b  share of area in contiguous same-class patches at least a given size,
#      stacked by class (CDL 2014-2023, patch_sizes/patches.csv from
#      crop_patch_sizes.R), with the resolution of FLEX, SIFnet, the downscaling
#      grid, OCO-2 and TROPOMI marked; bands and legend boxes are shares of the whole
#      study area (developed and other classes in the legend only)
#
# Output: fig1_manuscript.png (600 dpi, at the printed width of 6.5 in)

library(sf)
library(terra)
library(data.table)
library(ggplot2)
library(patchwork)
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "fig_fonts.R"))   # Barlow, save_fig()

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
ctx <- file.path(proj_dir, "context_data")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
site_ll <- c(lon = -97.486, lat = 36.604)                   # Lamont TCCON / ARM SGP CF
pool_years <- 2014:2023                                     # as in crop_patch_sizes.R
x_min <- 4e-4
stem_col <- "#1f5f7a"
crop_class_colors <- c(
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
  "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e",
  "Fallow/Idle Cropland" = "#bdb76b", "Trees/Wetlands" = "#1f4e3d", "Developed" = "#8c8c8c", "Other" = "#cccccc")
legend_names <- c("Winter Wheat" = "Winter wheat", "Dbl Crop WinWht/Soybeans" = "Wheat/soy (double)", "Fallow/Idle Cropland" = "Fallow/idle",
                  "Grassland/Pasture" = "Grassland/pasture", "Trees/Wetlands" = "Trees/wetlands")
base <- 9.5                                                 # base font size (pt) at the printed size (6.5 in wide)

grid <- st_read(grid_file, quiet = TRUE)
crs_g <- st_crs(grid)
domain <- st_union(st_geometry(grid))

# ---- a, top: locator map --------------------------------------------------------
states <- st_read(file.path(ctx, "cb_2023_us_state_20m.shp"), quiet = TRUE)
counties <- st_read(file.path(ctx, "cb_2023_us_county_20m.shp"), quiet = TRUE)
ok_cty <- st_transform(counties[counties$STUSPS == "OK", ], 4326)
states <- st_transform(states, 4326)
dom_ll <- st_transform(domain, 4326)
lab <- data.frame(label = c("CO", "KS", "MO", "NM", "OK", "AR", "TX", "LA"),
                  x = c(-102.6, -98.2, -93.95, -103.45, -97.5, -93.9, -99.9, -93.45),
                  y = c(37.45, 37.45, 37.2, 35.0, 35.55, 35.0, 34.0, 32.85))
halo <- function(d, size) {      # white halo behind dark text
  o <- expand.grid(dx = c(-1, 0, 1) * 0.035, dy = c(-1, 0, 1) * 0.03)
  h <- do.call(rbind, lapply(seq_len(nrow(o)), function(i) transform(d, x = x + o$dx[i], y = y + o$dy[i])))
  list(geom_text(data = h, aes(x, y, label = label), color = "white", size = size),
       geom_text(data = d, aes(x, y, label = label), color = "grey15", size = size))
}
pa1 <- ggplot() +
  geom_sf(data = states, fill = "grey76", color = "grey10", linewidth = 0.2) +
  geom_sf(data = ok_cty, fill = NA, color = "grey13", linewidth = 0.45) +
  geom_sf(data = states[states$STUSPS == "OK", ], fill = NA, color = "grey10", linewidth = 0.6) +
  geom_sf(data = st_as_sfc(st_bbox(dom_ll)), fill = "red", color = "red", linewidth = 0.3) +
  halo(lab, size = (base - 2.5) / .pt) +
  coord_sf(crs = 4326, xlim = c(-103.85, -93.25), ylim = c(32.6, 37.85), expand = FALSE,
           label_axes = list(top = "E", bottom = "E", left = "N", right = "N")) +
  scale_x_continuous(breaks = seq(-102, -94, 2)) +
  scale_y_continuous(breaks = c(34, 36)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = base, base_family = fig_font) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey76"),
        panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.3),
        axis.text = element_text(size = base - 2, color = "grey10"),
        axis.ticks = element_line(color = "grey20", linewidth = 0.3), axis.ticks.length = unit(1.5, "pt"),
        plot.margin = margin(0, 0, 4, 0))

# ---- a, bottom: crop map ----------------------------------------------------------
cls <- rast(file.path(proj_dir, "cdl2024_classes_10m.tif"))
v <- as.matrix(cls, wide = TRUE)
lv <- cats(cls)[[1]]
pal <- crop_class_colors[lv$crop_class]
img <- matrix(c(pal, NA)[ifelse(is.na(v), length(pal) + 1, match(v, lv$value))], nrow = nrow(v))
img[is.na(img)] <- "#FFFFFF00"
e <- as.vector(ext(cls))
bb <- st_bbox(domain)
view <- c(bb[["xmin"]] - 700, bb[["xmax"]] + 700, bb[["ymin"]] - 5600, bb[["ymax"]] + 500)
cty4 <- st_transform(st_read(file.path(proj_dir, "lamont_counties.gpkg"), quiet = TRUE), crs_g)
cty_lines <- st_intersection(st_cast(st_geometry(cty4), "MULTILINESTRING"), domain)   # county lines inside the domain only
site <- st_coordinates(st_transform(st_sfc(st_point(site_ll), crs = 4326), crs_g))
star <- function(x, y, r_out, r_in = r_out * 0.42) {       # five-pointed star polygon
  a <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 11)[1:10]
  r <- rep(c(r_out, r_in), 5)
  data.frame(x = x + r * cos(a), y = y + r * sin(a))
}
sb_x0 <- bb[["xmin"]]; sb_y <- bb[["ymin"]] - 1900; sb_h <- 350
sb <- data.frame(xmin = sb_x0 + (0:1) * 5000, xmax = sb_x0 + (1:2) * 5000, fill = c("black", "white"))
pa2 <- ggplot() +
  annotation_raster(img, xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], interpolate = FALSE) +
  geom_sf(data = grid, fill = NA, color = "grey40", linewidth = 0.12) +
  geom_sf(data = cty_lines, color = "black", linewidth = 0.35) +
  geom_polygon(data = star(site[1], site[2], 650), aes(x, y), fill = "black", color = "white", linewidth = 0.4) +
  geom_rect(data = sb, aes(xmin = xmin, xmax = xmax, ymin = sb_y, ymax = sb_y + sb_h, fill = fill),
            color = "black", linewidth = 0.25) +
  annotate("text", x = sb_x0 + 10000 + 500, y = sb_y + sb_h / 2, label = "10 km", hjust = 0, vjust = 0.5,
           size = fs_tick / .pt, color = "grey10") +
  scale_fill_identity() +
  coord_sf(crs = crs_g, xlim = view[1:2], ylim = view[3:4], expand = FALSE, datum = NA) +
  theme_void(base_size = base, base_family = fig_font) + theme(plot.margin = margin(0, 0, 0, 0))

# ---- b: stacked patch-size curve ---------------------------------------------------
p <- fread(file.path(proj_dir, "patch_sizes", "patches.csv"))[year %in% pool_years]
cls_area <- p[, .(area = sum(area_km2)), by = crop_class]
main <- cls_area[!crop_class %in% c("Developed", "Other") & area / length(pool_years) > 5][order(-area), crop_class]
total_area <- cls_area[, sum(area)]      # all classes: bands and legend are shares of the whole study area
xs <- 10^seq(log10(x_min), log10(40), length.out = 400)
st <- p[crop_class %in% main, {
  a <- sort(area_km2, decreasing = TRUE); cs <- cumsum(a)
  k <- vapply(xs, function(x) sum(a >= x), integer(1))
  .(x = xs, share = ifelse(k > 0, cs[pmax(k, 1)], 0) / total_area)
}, by = crop_class]
st[, crop_class := factor(crop_class, levels = rev(main))]
tot <- st[, .(total = sum(share)), by = x]

deg2_km2 <- function(d, lat = 36.6) (d * 111.32 * cos(lat * pi / 180)) * (d * 110.95)
cell_km2 <- median(as.numeric(st_area(grid))) / 1e6
refs <- data.table(label = c("FLEX", "SIFNet", "Grid", "OCO-2", "TROPOMI"),
                   x = c(0.3 * 0.3, deg2_km2(0.005), cell_km2, 1.29 * 2.25, 3.5 * 5.5),
                   y_user = c(0.95, 0.82, 0.63, 0.35, 0.10))         # staggered as in fig1_v3.pptx
refs[, curve := vapply(x, function(r) tot[which.min(abs(log10(x) - log10(r))), total], numeric(1))]
refs[, y := pmax(y_user, curve + 0.06)]                              # stem always ends above the stack
pb1 <- ggplot(st, aes(x, share, fill = crop_class)) +
  geom_area(position = "stack", outline.type = "full", color = "white", linewidth = 0.1) +
  geom_segment(data = refs, aes(x = x, xend = x, y = 0, yend = y), inherit.aes = FALSE,
               color = stem_col, linewidth = 1.1, lineend = "butt") +
  geom_point(data = refs, aes(x, y), inherit.aes = FALSE, color = stem_col, size = 2) +
  geom_text(data = refs, aes(x, y, label = label), inherit.aes = FALSE, vjust = -0.75, size = base / .pt,
            color = "grey5") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10), labels = c("0.001", "0.01", "0.1", "1", "10"),
                limits = c(x_min, 40), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = crop_class_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1.07), clip = "off") +
  labs(x = expression(Contiguous~patch~area~(km^2)), y = "Share of study area in\npatches this size or larger") +
  theme_classic(base_size = base, base_family = fig_font) +
  theme(axis.text = element_text(color = "grey10", size = base - 1), plot.margin = margin(14, 6, 4, 2))

# legend: two columns of colored boxes with each class's share of the study area. Developed and
# Other are one legend-only entry: they are not plotted, because roads join into very large
# connected "patches" that would misrepresent field sizes
leg <- cls_area[crop_class %in% main][order(-area)][, `:=`(pct = sprintf("%.0f%%", 100 * area / total_area),
                                                          name = ifelse(crop_class %in% names(legend_names),
                                                                        legend_names[crop_class], crop_class))]
leg <- rbind(leg, data.table(crop_class = "Other", area = cls_area[!crop_class %in% main, sum(area)],
                             name = "Developed/other"), fill = TRUE)
leg[crop_class == "Other", pct := sprintf("%.0f%%", 100 * area / total_area)]
leg[, txt_col := ifelse(crop_class == "Other", "grey10", "white")]
n_left <- ceiling(nrow(leg) / 2)
leg[, `:=`(col = ifelse(.I <= n_left, 0, 1), row = ifelse(.I <= n_left, .I, .I - n_left))]
box_w <- 0.2; col_x <- c(0, 1.1)
pb2 <- ggplot(leg) +
  geom_tile(aes(x = col_x[col + 1] + box_w / 2, y = -row, fill = crop_class), width = box_w, height = 0.78,
            color = "#1d2b36", linewidth = 0.35) +
  geom_text(aes(x = col_x[col + 1] + box_w / 2, y = -row, label = pct, color = txt_col), size = (fs_tick - 1) / .pt) +
  scale_color_identity() +
  geom_text(aes(x = col_x[col + 1] + box_w + 0.05, y = -row, label = name), hjust = 0, size = fs_tick / .pt,
            color = "grey5") +
  scale_fill_manual(values = crop_class_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.02, 2.12), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-n_left - 0.5, -0.5), expand = c(0, 0)) +
  theme_void(base_size = base, base_family = fig_font) + theme(plot.margin = margin(6, 0, 0, 20))

# ---- assemble ---------------------------------------------------------------------
# a: crop map, with a small Oklahoma locator (counties, study area in red) set in below the domain;
# b: patch-area curve; one shared legend across the bottom (the colours serve both panels)
ok_state <- states[states$STUSPS == "OK", ]
locator <- ggplot() +
  geom_sf(data = ok_cty, fill = "grey88", color = "grey55", linewidth = 0.12) +
  geom_sf(data = ok_state, fill = NA, color = "grey15", linewidth = 0.35) +
  geom_sf(data = st_as_sfc(st_bbox(dom_ll)), fill = "red", color = "red", linewidth = 0.2) +
  coord_sf(crs = 4326, expand = FALSE, datum = NA) +
  theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
pa <- pa2 + labs(tag = "a.") + inset_element(locator, left = 0.56, right = 1, bottom = 0, top = 0.15, align_to = "panel",
                          ignore_tag = TRUE)
legend_strip <- local({
  ncol <- 4; d <- copy(leg)
  d[, `:=`(col = (.I - 1) %/% 3, row = (.I - 1) %% 3)]
  bw <- 0.2                                                   # box width, in column units
  ggplot(d) +
    geom_tile(aes(x = col + bw / 2, y = -row, fill = crop_class), width = bw, height = 0.78,
              color = "#1d2b36", linewidth = 0.3) +
    geom_text(aes(x = col + bw / 2, y = -row, label = pct, color = txt_col), size = (fs_tick - 1) / .pt,
              family = fig_font) +
    geom_text(aes(x = col + bw + 0.05, y = -row, label = name), hjust = 0, size = fs_tick / .pt,
              family = fig_font, color = "grey5") +
    scale_color_identity() +
    scale_fill_manual(values = crop_class_colors, guide = "none") +
    scale_x_continuous(limits = c(-0.02, ncol), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-2.5, 0.5), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(6, 4, 0, 4))
})
fig <- (pa | pb1 + labs(tag = "b.")) / legend_strip +
  plot_layout(heights = c(1, 0.19)) &
  theme(plot.tag = element_text(size = base + 0.5, face = "bold", family = fig_font),
        plot.tag.position = c(0, 1)) & font_theme
fig[[1]] <- fig[[1]] + plot_layout(widths = c(0.86, 1))
save_fig(fig, "fig1_manuscript", height = 3.6)

cat("legend shares (% of the study area):\n"); print(leg[, .(crop_class, pct)])
cat("reference scales (km2) and stem tops:\n"); print(refs[, .(label, x = signif(x, 3), curve = round(curve, 2), y)])

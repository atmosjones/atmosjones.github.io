# Effective resolution of the grid inversion across all targets, from the widths
# of the averaging-kernel rows that downscale_sif_inversion.R saves per cell
# (res_area_pr_km2: participation area of row j; see that script's header).
#
# Per target: median resolution of the constrained cells (averaging-kernel
# diagonal > 0.5) and the grid area resolved finer than a nominal OCO-2
# footprint (2.9 km2), counting only cells whose estimate comes from data
# (row sum > 0.5) and lying more than edge_km from the grid edge, where rows are
# truncated by the boundary and look narrower than they are. Per cell: median resolution over targets, against the median
# number of overlapping soundings.
# Output (effective_resolution/): resolution_by_target.csv, resolution_by_cell.csv,
#   plot_effective_resolution.png, plot_resolution_vs_overlap.png, printed summary.

library(sf)
library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "effective_resolution")
dir.create(out_dir, showWarnings = FALSE)
nominal_fp <- 2.9                     # km2, OCO-2 nominal footprint (1.29 x 2.25 km)
edge_km <- 1.5

grid <- st_read(file.path(proj_dir, "sub_road_grid.shp"), quiet = TRUE)
grid$cell_id <- seq_len(nrow(grid))
grid$area <- as.numeric(st_area(grid)) / 1e6
edge <- st_cast(st_union(st_geometry(grid)), "MULTILINESTRING")
grid$interior <- as.numeric(st_distance(st_centroid(st_geometry(grid)), edge)) > edge_km * 1000
cells <- as.data.table(st_drop_geometry(st_read(file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg"),
  query = "SELECT cell_id, date, avg_kernel, n_footprints, res_area_km2, res_area_pr_km2, ak_row_sum FROM cells",
  quiet = TRUE)))
cells[, area := grid$area[cell_id]]
cells[, data_driven := ak_row_sum > 0.5 & grid$interior[cell_id]]

by_target <- cells[, .(n_constrained = sum(avg_kernel > 0.5),
                       res_constrained = median(res_area_pr_km2[avg_kernel > 0.5]),
                       res_fwhm_constrained = median(res_area_km2[avg_kernel > 0.5]),
                       area_data = sum(area[data_driven]),
                       area_finer = sum(area[data_driven & res_area_pr_km2 < nominal_fp]),
                       area_finer_1km2 = sum(area[data_driven & res_area_pr_km2 < 1])), by = date]
by_target[, frac_finer := area_finer / sum(grid$area)]
fwrite(by_target, file.path(out_dir, "resolution_by_target.csv"))

by_cell <- cells[, .(res_med = median(res_area_pr_km2[ak_row_sum > 0.5]), n_fp_med = as.numeric(median(n_footprints)),
                     n_fp_med_dd = as.numeric(median(n_footprints[ak_row_sum > 0.5])),   # same targets as res_med
                     ak_med = median(avg_kernel), share_data = mean(data_driven)), by = cell_id]
fwrite(by_cell, file.path(out_dir, "resolution_by_cell.csv"))

q <- function(x) sprintf("%.2f (IQR %.2f-%.2f)", median(x, na.rm = TRUE),
                         quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
cat(sprintf("\n%d targets, grid %.0f km2\n", nrow(by_target), sum(grid$area)))
cat("resolution of constrained cells (participation area, km2):", q(by_target$res_constrained), "\n")
cat("  as a length (sqrt, km):", q(sqrt(by_target$res_constrained)), "\n")
cat("  FWHM area of the same cells (km2):", q(by_target$res_fwhm_constrained), "\n")
cat("grid area resolved finer than 2.9 km2 (km2):", q(by_target$area_finer), "\n")
cat("  share of the grid:", q(by_target$frac_finer), "\n")
cat("grid area resolved finer than 1 km2 (km2):", q(by_target$area_finer_1km2), "\n")
cat(sprintf("cell-level: log resolution vs log overlapping soundings r = %.2f (old: all-target soundings, all cells)\n",
            by_cell[n_fp_med > 0 & is.finite(res_med), cor(log(res_med), log(n_fp_med))]))
# resolution against overlap on the same targets, interior cells only (edge rows are truncated)
sc <- merge(by_cell, data.table(cell_id = grid$cell_id, interior = grid$interior), by = "cell_id")[
  interior == TRUE & is.finite(res_med) & n_fp_med_dd > 0]
cat(sprintf("interior cells, matched targets: n = %d, Spearman rho = %.2f, Pearson r (log-log) = %.2f\n", nrow(sc),
            sc[, cor(res_med, n_fp_med_dd, method = "spearman")], sc[, cor(log(res_med), log(n_fp_med_dd))]))
print(sc[, .(cells = .N, res_median = median(res_med)),
         keyby = .(soundings = cut(n_fp_med_dd, c(0, 10, 20, 40, 60, 80, 120, Inf)))])
print(cells[data_driven == TRUE, .(cells = .N, res_med = median(res_area_pr_km2)),
            keyby = .(n_soundings = cut(n_footprints, c(-1, 0, 5, 10, 20, 40, 80, Inf)))])

# ---- figure: median over targets --------------------------------------------
g <- merge(grid[, c("cell_id", "interior", "geometry")], by_cell, by = "cell_id")
fine <- st_union(g[!is.na(g$res_med) & g$res_med < nominal_fp & g$interior, ])
map_theme <- theme_void(base_size = 9) +
  theme(panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.4),
        plot.title = element_text(size = 9, face = "bold"))
pa <- ggplot(g) + geom_sf(aes(fill = n_fp_med), color = NA) +
  scale_fill_viridis_c(option = "rocket", direction = -1, trans = "sqrt", name = "Soundings\nper cell") +
  labs(title = "Median overlapping soundings per cell") + map_theme
pb <- ggplot(g) + geom_sf(aes(fill = res_med), color = NA) +
  geom_sf(data = st_sf(geometry = fine), fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "mako", trans = "log10", limits = c(0.64, 12), oob = scales::squish,
                       breaks = c(0.64, 1, 2.9, 10), labels = c("0.64 (cell)", "1", "2.9 (OCO-2)", "10"),
                       na.value = "grey85", name = expression(km^2)) +
  labs(title = expression(bold("Median effective resolution (outline: finer than 2.9" ~ km^2 * ")"))) + map_theme
pc <- ggplot(by_target, aes(frac_finer * 100)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = "#2c6e8f", color = "white") +
  labs(x = "Grid area resolved finer than a nominal OCO-2 footprint (%)", y = "Targets",
       title = sprintf("Per target (median %.0f%%)", 100 * median(by_target$frac_finer))) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(out_dir, "plot_effective_resolution.png"), pa + pb + pc, width = 13, height = 4.4, dpi = 300, bg = "white")

# ---- figure: resolution against overlap (per cell, same targets) ------------
source(file.path(proj_dir, "fig_fonts.R"))
bins <- sc[, .(n = median(n_fp_med_dd), res = median(res_med), lo = quantile(res_med, 0.25), hi = quantile(res_med, 0.75)),
           by = .(b = cut(n_fp_med_dd, c(0, 10, 20, 40, 60, 80, 120, Inf)))]
pd <- ggplot(sc, aes(n_fp_med_dd, res_med)) +
  geom_hline(yintercept = nominal_fp, linetype = "22", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 1, y = nominal_fp, label = "nominal OCO-2 footprint (2.9 km²)", hjust = 0, vjust = -0.5,
           size = (fs_tick - 1) / .pt, family = fig_font, color = "grey40") +
  geom_point(size = 0.7, alpha = 0.45, color = "#2c6e8f", stroke = 0) +
  geom_errorbar(data = bins, aes(n, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0, linewidth = 0.45) +
  geom_point(data = bins, aes(n, res), inherit.aes = FALSE, size = 1.6) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  labs(x = "Overlapping soundings (median per cell)", y = expression("Effective resolution (km"^2*")")) +
  theme_bw(base_size = fs_title, base_family = fig_font) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = fs_tick - 0.5))
ggsave(file.path(out_dir, "plot_resolution_vs_overlap.png"), pd, width = 4.5, height = 3, dpi = 600, bg = "white",
       device = ragg::agg_png)

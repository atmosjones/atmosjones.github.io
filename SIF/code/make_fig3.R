# Figure 3 for the manuscript: what the two inversions recover for each crop.
#
#   a  Seasonal cycles by crop, 2014-2025. Crop-class inversion
#      (downscale_sif_by_crop_cdl30.R): every target (points, averaging kernel
#      >= 0.9) and its seasonal curve. ~800 m grid inversion
#      (downscale_sif_inversion.R): seasonal curve of the per-target median of
#      the cells whose majority crop is that class, using cells with averaging
#      kernel >= min_ak_grid (mostly set by the data) and >= 3 cells per target.
#      Grey bands: median active season from USDA NASS Crop Progress (active_season_bands.R);
#      grassland/pasture has no NASS stages. Triangles: median 50% date of each crop's
#      reproductive stage (wheat heading, corn silking, sorghum heading, soybean blooming,
#      cotton setting bolls), bar = range over 2014-2025. Green ribbon: 95% bootstrap
#      range of the crop-class curve (bootstrap_climatology.R, 1,000 resamples of targets).
#   b  Within-day crop contrasts (crop minus that target's mean over the crops
#      present in both): grid medians against crop-class values, all targets.
#
# Both curves use the seasonal fit of drought_analysis.R (loess on day of year,
# span 0.4 year, data replicated a year either side); crop-class points are
# weighted by 1/SD^2, grid medians by their number of cells.
# Output: fig3_manuscript.png / .pdf

library(sf)
library(data.table)
library(ggplot2)
library(patchwork)
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "fig_fonts.R"))   # Barlow, sizes, crop_short, save_fig()
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "strip_colors.R"))   # color_strips()
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "active_season_bands.R"))

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
             "Sorghum", "Cotton", "Grassland/Pasture")
crop_class_colors <- c(   # as in make_plots.R
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
  "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e")
min_ak_class <- 0.9
min_ak_grid  <- 0.5
min_cells    <- 3
loess_span   <- 0.4
sif_label <- expression(SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1))

seasonal_curve <- function(doy, y, w, new_doy = 1:365) {
  d <- data.table(doy, y, w)
  wrap <- rbind(d, copy(d)[, doy := doy - 365], copy(d)[, doy := doy + 365])
  predict(loess(y ~ doy, data = wrap, weights = w, span = loess_span / 3, degree = 2),
          newdata = data.frame(doy = new_doy))
}

# ---- crop-class inversion -----------------------------------------------------
cc <- fread(file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv"))
cc <- cc[crop_class %in% classes & avg_kernel >= min_ak_class][, date := as.Date(date)]
cc[, doy := as.integer(format(date, "%j"))]

# ---- grid inversion, by majority crop -------------------------------------------
g <- as.data.table(st_drop_geometry(st_read(file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg"),
                                            query = "SELECT date, sif757, avg_kernel, crop_name FROM cells",
                                            quiet = TRUE)))
g <- g[crop_name %in% classes & avg_kernel >= min_ak_grid][, date := as.Date(date)]
gm <- g[, .(grid_med = median(sif757), n_cells = .N), by = .(date, crop_class = crop_name)][n_cells >= min_cells]
gm[, doy := as.integer(format(date, "%j"))]

curves <- rbind(
  cc[, .(doy = 1:365, sif = seasonal_curve(doy, sif757, 1 / sif757_sd_scaled^2),
         method = "Crop-class inversion"), by = crop_class],
  gm[, .(doy = 1:365, sif = seasonal_curve(doy, grid_med, n_cells),
         method = "Grid inversion\n(~800-m cells,\nmajority crop)"), by = crop_class])
for (d in list(cc, gm, curves)) d[, crop_class := factor(crop_class, levels = classes)]
meth_cols <- c("Crop-class inversion" = "#1a9850", "Grid inversion\n(~800-m cells,\nmajority crop)" = "grey25")

bands <- active_season_bands()[crop_class %in% classes][, crop_class := factor(crop_class, levels = classes)]
handoff <- bands[part == 2]                     # double crop: wheat harvest, soybeans follow
band_layers <- list(
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                              fill = "Active season\n(USDA, median)"),
            color = NA, inherit.aes = FALSE),
  scale_fill_manual(values = c("Active season\n(USDA, median)" = "grey88"), name = NULL),
  geom_vline(data = handoff, aes(xintercept = xmin), linetype = 3, color = "grey55", linewidth = 0.4,
             show.legend = FALSE))
# reproductive stage of each crop: median 50% date over 2014-2025 (USDA NASS Crop Progress)
# (double crop: wheat heading, and soybean blooming from the all-soybean NASS series, which
# double-cropped soybeans, planted after the wheat harvest, probably reach later)
stage_of <- data.table(crop_class = c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Dbl Crop WinWht/Soybeans",
                                      "Corn", "Sorghum", "Soybeans", "Cotton"),
                       commodity = c("WHEAT", "WHEAT", "SOYBEANS", "CORN", "SORGHUM", "SOYBEANS", "COTTON"),
                       stage = c("HEADED", "HEADED", "BLOOMING", "SILKING", "HEADED", "BLOOMING", "SETTING BOLLS"),
                       label = c("headed", "headed", "soy bloom", "silking", "headed", "blooming", "setting bolls"))
x50 <- fread(file.path(proj_dir, "crop_progress", "stage_dates_50pct.csv"))[, doy := as.numeric(format(as.Date(d50), "%j"))]
stage_marks <- merge(stage_of, x50[, .(doy = median(doy), doy_min = min(doy), doy_max = max(doy)), by = .(commodity, stage)],
                     by = c("commodity", "stage"))[, crop_class := factor(crop_class, levels = classes)]
# the double crop's two labels sit close: push them apart
stage_marks[, hj := fifelse(crop_class == "Dbl Crop WinWht/Soybeans", fifelse(stage == "HEADED", 0.85, 0.15), 0.5)]
env <- fread(file.path(proj_dir, "bootstrap_climatology", "envelope.csv"))[crop_class %in% classes]
env[, crop_class := factor(crop_class, levels = classes)]
pa <- ggplot() +
  band_layers +
  geom_errorbarh(data = stage_marks, aes(xmin = doy_min, xmax = doy_max, y = 1.52), height = 0, color = "grey35",
                 linewidth = 0.4, inherit.aes = FALSE) +
  geom_point(data = stage_marks, aes(doy, 1.52), shape = 25, size = 1.2, fill = "grey20", color = "grey20",
             inherit.aes = FALSE) +
  geom_text(data = stage_marks, aes(doy, 1.52, label = label, hjust = hj), vjust = -0.7, size = (fs_tick - 1) / .pt, color = "grey20",
            inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3, show.legend = FALSE) +
  geom_point(data = cc, aes(doy, sif757), color = "#1a9850", size = 0.35, alpha = 0.45) +
  # 95% bootstrap range of the crop-class curve (bootstrap_climatology.R; same fit as here)
  geom_ribbon(data = env, aes(doy, ymin = lo, ymax = hi), fill = "#1a9850", alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = curves, aes(doy, sif, color = method, linetype = method), linewidth = 0.55) +
  facet_wrap(~crop_class, ncol = 4, labeller = as_labeller(crop_short)) +
  scale_color_manual(values = meth_cols) +
  scale_linetype_manual(values = c("Crop-class inversion" = 1,
                                   "Grid inversion\n(~800-m cells,\nmajority crop)" = 2)) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("Jan", "Apr", "Jul", "Oct")) +
  coord_cartesian(ylim = c(-0.25, 1.7)) +
  labs(x = NULL, y = sif_label, color = NULL, linetype = NULL,
       title = sprintf("Seasonal cycles by crop, %s-%s (points: crop-class estimates on each target)",
                       format(min(cc$date), "%Y"), format(max(cc$date), "%Y"))) +
  theme_bw(base_size = fs_title, base_family = fig_font) +
  theme(strip.text = element_text(size = fs_tick, face = "bold", margin = margin(2, 0, 2, 0)), legend.position = "inside",
        axis.text = element_text(size = fs_tick - 0.5), legend.text = element_text(size = fs_tick - 0.5),
        legend.position.inside = c(0.775, 0.24), legend.justification.inside = c(0, 0.5),
        legend.key.spacing.y = unit(3, "pt"), legend.key.size = unit(0.4, "cm"), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.background = element_rect(fill = NA))

# ---- b. within-day contrasts --------------------------------------------------
m <- merge(gm, cc[, .(date, crop_class, cls = sif757)], by = c("date", "crop_class"))
m[, n_crops := .N, by = date]
m <- m[n_crops >= 3]
m[, `:=`(grid_anom = grid_med - mean(grid_med), cls_anom = cls - mean(cls)), by = date]
r <- cor(m$cls_anom, m$grid_anom)
slope <- coef(lm(grid_anom ~ cls_anom, m))[[2]]
m[, crop_class := factor(crop_class, levels = classes)]
lim <- range(c(m$cls_anom, m$grid_anom))
pb <- ggplot(m, aes(cls_anom, grid_anom, color = crop_class)) +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 1.1, alpha = 0.75) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.5) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 3,
           label = sprintf("r = %.2f\nslope = %.2f\nn = %d crop-targets", r, slope, nrow(m))) +
  scale_color_manual(values = crop_class_colors) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Crop-class inversion: crop minus target mean",
       y = "Grid inversion: crop median minus target mean", color = NULL,
       title = "Within-day crop contrasts") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"))

# Figure 3 is panel a only; the within-day contrasts (pb) are a facet of Figure 4b
# (make_fig4.R reads grid_vs_class_contrasts.csv). pb is still drawn here for checking.
fwrite(m[, .(date, crop_class, cls_anom, grid_anom)], file.path(proj_dir, "grid_vs_class_contrasts.csv"))
ggsave(file.path(proj_dir, "plot_grid_vs_class_contrasts.png"), pb, width = 6, height = 4.5, dpi = 300, bg = "white")
fig <- wrap_elements(full = color_strips(pa + labs(title = NULL), crop_class_colors, classes, ncol = 4))
save_fig(fig, "fig3_manuscript", height = 4.2)

# numbers for the text
pk <- curves[, .(peak_doy = doy[which.max(sif)], peak = max(sif)), by = .(crop_class, method)]
# SIF peak timing against the NASS stages (crop-class curve), for the text / a table
cmp <- merge(pk[method == "Crop-class inversion", .(crop_class, sif_peak_doy = peak_doy)], stage_marks[, .(crop_class, stage = label,
             stage_doy = doy, stage_range = paste(doy_min, doy_max, sep = "-"))], by = "crop_class")
cmp[, peak_minus_stage_days := sif_peak_doy - stage_doy]
fwrite(cmp, file.path(proj_dir, "crop_progress", "sif_peak_vs_stage.csv"))
print(cmp)
print(dcast(pk, crop_class ~ method, value.var = c("peak_doy", "peak")), digits = 3)
cat(sprintf("within-day contrasts: r = %.2f, slope = %.2f, n = %d; grid medians from %d cell-targets\n",
            r, slope, nrow(m), nrow(g)))

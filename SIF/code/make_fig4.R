# Figure 4 for the manuscript: crop-level comparison of the target-mode inversion
# with the grid inversion, SIFnet and HLS NDVI.
#   Figure 4a = panel "c" below (seasonal cycles), Figure 4b = panel "b" (contrasts).
# The footprint-level panel is Figure S8. (GOSIF was dropped from the paper.)
#
#   S8 Footprint level: per-target correlation of each predictor with the OCO-2
#      footprint SIF, beside the cross-validated crop-class model on the same
#      targets and the noise ceiling (the best correlation possible given the
#      retrieval noise). Each product is compared on its own matched targets.
#   b  Crop level: within-day crop contrasts (class minus that day's class mean)
#      of the grid inversion (median of majority-crop cells; make_fig3.R,
#      grid_vs_class_contrasts.csv), SIFnet and HLS NDVI against the OCO-2
#      crop-class inversion, over the same seven classes (crops +
#      grassland/pasture) and class-days with an averaging kernel >= 0.9. Products on
#      y, inversion on x, so each fitted slope is the share of the contrast recovered;
#      1:1 line where both axes are SIF (grid, SIFnet).
#   c  Seasonal cycles by crop at 740 nm over the SIFnet period only (Apr 2018 -
#      Apr 2021): crop-class inversion on each target (points; avg. kernel >= 0.9)
#      and SIFnet's class mean in every 16-day composite (points), each with
#      the seasonal fit of make_fig3.R. Same band and years, no rescaling.
#      HLS NDVI class means from every clear scene in the period (dashed, right
#      axis; fixed mapping NDVI = SIF / ndvi_scale), for timing and shape only.
#      Grey bands: median active season from USDA NASS Crop Progress (active_season_bands.R).
#
# Needs: compare_sifnet.R, compare_ndvi.R and make_fig3.R outputs and
#   downscale_sif_by_crop_cdl30.R run with band 740.
# Output: fig4_manuscript.png / .pdf; figS8_footprint_skill.png / .pdf

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
min_avg_kernel <- 0.9

# ---- a. footprint level -------------------------------------------------------
sn <- fread(file.path(proj_dir, "sifnet_comparison", "footprint_stats.csv"))
nd <- fread(file.path(proj_dir, "ndvi_comparison", "footprint_stats.csv"))
lab <- function(name, band, n) sprintf("%s targets\n(%s nm, n = %d)", name, band, n)
set_s <- lab("SIFnet", 740, nrow(sn)); set_n <- lab("HLS NDVI", 757, nrow(nd))
fa <- rbind(
  sn[, .(set = set_s, what = "Crop-class model (CV)", r = r_crop_cv)],
  sn[, .(set = set_s, what = "Product", r = r_sifnet)],
  sn[, .(set = set_s, what = "Noise ceiling", r = noise_ceiling)],
  nd[, .(set = set_n, what = "Crop-class model (CV)", r = r_crop_cv)],
  nd[, .(set = set_n, what = "Product", r = r_ndvi_cv)],
  nd[, .(set = set_n, what = "Crop classes + NDVI (CV)", r = r_crop_ndvi_cv)],
  nd[, .(set = set_n, what = "Noise ceiling", r = noise_ceiling)])
fa[, set := factor(set, levels = c(set_s, set_n))]
fa[, what := factor(what, levels = c("Crop-class model (CV)", "Product", "Crop classes + NDVI (CV)", "Noise ceiling"))]
med <- fa[!is.na(r), .(m = median(r)), by = .(set, what)]
pa <- ggplot(fa[!is.na(r)], aes(x = set, y = r, fill = what)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.6, width = 0.75, position = position_dodge(width = 0.85)) +
  scale_fill_manual(values = c("Crop-class model (CV)" = "#1a9850", "Product" = "grey60",
                               "Crop classes + NDVI (CV)" = "#91cf60", "Noise ceiling" = "white"),
                    labels = c("Crop-class model (cross-validated)", "SIFnet / NDVI",
                               "Crop classes + NDVI (cross-validated)", "Noise ceiling")) +
  labs(x = NULL, y = "Per-target correlation with OCO-2 footprint SIF", fill = NULL,
       title = "Footprint level") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 2))

# ---- b. within-day crop contrasts --------------------------------------------
contr <- function(d, y, x, name, xlab) {
  d <- d[avg_kernel >= min_avg_kernel & crop_class %in% classes & !is.na(get(x)) & !is.na(get(y))]
  d[, `:=`(oco = get(y) - mean(get(y)), prod = get(x) - mean(get(x))), by = date]
  d[, .(date, crop_class, oco, prod, product = name, xlab = xlab)]
}
cb <- rbind(
  contr(fread(file.path(proj_dir, "sifnet_comparison", "class_comparison.csv")), "sif740", "sifnet_sif740",
        "SIFnet (740 nm)", "SIFnet: class minus day mean"),
  contr(fread(file.path(proj_dir, "ndvi_comparison", "class_comparison.csv")), "sif757", "ndvi",
        "HLS NDVI", "NDVI: class minus day mean"))
grid_c <- fread(file.path(proj_dir, "grid_vs_class_contrasts.csv"))      # make_fig3.R
cb <- rbind(grid_c[, .(date = as.Date(date), crop_class, oco = cls_anom, prod = grid_anom,
                       product = "Grid inversion (~800 m)", xlab = "Grid: class minus day mean")],
            cb[, date := as.Date(date)])
cb[, product := factor(product, levels = c("Grid inversion (~800 m)", "SIFnet (740 nm)", "HLS NDVI"))]
cb[, crop_class := factor(crop_class, levels = classes)]
# product regressed on the inversion: the slope is the share of the inversion's contrast
# the product recovers (as in contrast_amplitude.R and the text)
stats <- cb[, { f <- coef(lm(prod ~ oco)); .(r = cor(oco, prod), n = .N, m = f[[2]], b = f[[1]], sd_ratio = sd(prod) / sd(oco)) },
            by = product]
# both axes are centred on each day's mean, so b is zero by construction and not shown
stats[, label := sprintf("y = %.2fx\nr = %.2f, n = %d", m, r, n)]
print(stats)

pb <- ggplot(cb, aes(x = oco, y = prod, color = crop_class)) +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  geom_abline(data = data.frame(product = factor(c("Grid inversion (~800 m)", "SIFnet (740 nm)"), levels = levels(cb$product))),
              aes(slope = 1, intercept = 0), linetype = 2, color = "grey50") +   # 1:1 where both axes are SIF
  geom_point(size = 0.45, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.4) +
  geom_text(data = stats, aes(label = label), x = min(cb$oco) - 0.02, y = max(cb$prod) + 0.02, hjust = 0, vjust = 1,
            inherit.aes = FALSE, size = (fs_tick - 1) / .pt, lineheight = 0.95) +
  facet_wrap(~product, nrow = 1) +
  # same scale on x and y (shared by all facets), so the 1:1 lines run at 45 degrees
  coord_fixed(ratio = 1, xlim = range(cb$oco), ylim = range(cb$prod), expand = TRUE) +
  scale_color_manual(values = crop_class_colors, labels = crop_short) +
  labs(x = expression("Crop-class inversion: class minus day mean"~(mW~m^-2~sr^-1~nm^-1)),
       y = "Product: class minus\nday mean",
       color = NULL, title = "Within-day crop contrasts") +
  theme_bw(base_size = fs_title, base_family = fig_font) +   # printed sizes, as in Figures 1-3
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        strip.text = element_text(size = fs_tick, margin = margin(2, 0, 2, 0)), axis.text = element_text(size = fs_tick - 0.5),
        legend.text = element_text(size = fs_tick - 0.5), legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 1.5, alpha = 1)))

# ---- c. seasonal cycles at 740 nm, SIFnet period --------------------------------
loess_span <- 0.4
seasonal_curve <- function(doy, y, w, new_doy = 1:365) {      # as in make_fig3.R
  d <- data.table(doy, y, w)
  wrap <- rbind(d, copy(d)[, doy := doy - 365], copy(d)[, doy := doy + 365])
  predict(loess(y ~ doy, data = wrap, weights = w, span = loess_span / 3, degree = 2),
          newdata = data.frame(doy = new_doy))
}
snm <- fread(file.path(proj_dir, "sifnet_comparison", "sifnet_class_means.csv"))[crop_class %in% classes]
snm[, `:=`(start = as.Date(sifnet_start), mid = as.Date(mid_date))]
period <- c(min(snm$start), max(snm$start) + 15)
inv <- fread(file.path(proj_dir, "downscaled_sif757", "sif740_by_crop_cdl30.csv"))[
  crop_class %in% classes & avg_kernel >= min_avg_kernel][, date := as.Date(date)][date >= period[1] & date <= period[2]]
inv[, doy := as.integer(format(date, "%j"))]
snm[, doy := as.integer(format(mid, "%j"))]
ndvi_scale <- 2                                        # SIF axis units per NDVI unit
hls <- fread(file.path(proj_dir, "ndvi_comparison", "ndvi_residuals.csv"))[
  crop_class %in% classes & !is.na(ndvi)][, date := as.Date(date)][date >= period[1] & date <= period[2]]
hls[, doy := as.integer(format(date, "%j"))]
src_inv <- "Crop-class inversion"; src_sn <- "SIFnet"
src_nd <- "NDVI (right axis)"
pts <- rbind(inv[, .(crop_class, doy, sif = sif740, source = src_inv)],
             snm[, .(crop_class, doy, sif = sifnet_sif740, source = src_sn)])
cur <- rbind(inv[, .(doy = 1:365, sif = seasonal_curve(doy, sif740, 1 / sif740_sd_scaled^2), source = src_inv), by = crop_class],
             snm[, .(doy = 1:365, sif = seasonal_curve(doy, sifnet_sif740, rep(1, .N)), source = src_sn), by = crop_class],
             hls[, .(doy = 1:365, sif = ndvi_scale * seasonal_curve(doy, ndvi, rep(1, .N)), source = src_nd), by = crop_class])
for (d in list(pts, cur)) d[, crop_class := factor(crop_class, levels = classes)]
src_cols <- setNames(c("#1a9850", "#5e3c99", "grey35"), c(src_inv, src_sn, src_nd))
src_lty <- setNames(c(1, 1, 2), c(src_inv, src_sn, src_nd))
bands <- active_season_bands()[crop_class %in% classes][, crop_class := factor(crop_class, levels = classes)]
handoff <- bands[part == 2]                     # double crop: wheat harvest, soybeans follow
band_layers <- list(
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                              fill = "Active season\n(USDA, median)"),
            color = NA, inherit.aes = FALSE),
  scale_fill_manual(values = c("Active season\n(USDA, median)" = "grey88"), name = NULL),
  geom_vline(data = handoff, aes(xintercept = xmin), linetype = 3, color = "grey55", linewidth = 0.4,
             show.legend = FALSE))
pc <- ggplot() +
  band_layers +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3, show.legend = FALSE) +
  geom_point(data = pts, aes(doy, sif, color = source), size = 0.35, alpha = 0.5) +
  geom_line(data = cur[source == src_nd], aes(doy, sif, color = source, linetype = source), linewidth = 0.45) +
  geom_line(data = cur[source != src_nd], aes(doy, sif, color = source, linetype = source), linewidth = 0.55) +
  facet_wrap(~crop_class, ncol = 4, labeller = as_labeller(crop_short)) +
  scale_color_manual(values = src_cols, breaks = c(src_inv, src_sn, src_nd)) +
  scale_linetype_manual(values = src_lty, breaks = c(src_inv, src_sn, src_nd)) +
  scale_y_continuous(sec.axis = sec_axis(~ . / ndvi_scale, name = "HLS NDVI")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = NULL, y = expression(SIF[740] ~~ (mW~m^-2~sr^-1~nm^-1)), color = NULL, linetype = NULL,
       title = sprintf("Seasonal cycles at 740 nm, %s to %s",   # counts are in the caption
                       format(period[1], "%b %Y"), format(period[2], "%b %Y"))) +
  theme_bw(base_size = fs_title, base_family = fig_font) +   # printed sizes, as in Figure 3
  theme(strip.text = element_text(size = fs_tick, face = "bold", margin = margin(2, 0, 2, 0)), legend.position = "inside",
        axis.text = element_text(size = fs_tick - 0.5), legend.text = element_text(size = fs_tick - 0.5),
        legend.key.spacing.y = unit(3, "pt"), legend.key.size = unit(0.4, "cm"), panel.grid.minor = element_blank(),
        legend.position.inside = c(0.815, 0.24), legend.justification.inside = c(0, 0.5), plot.title = element_text(face = "bold"),
        legend.background = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA), size = 1.5, alpha = 1)))
pk <- cur[source != src_nd, .(peak_doy = doy[which.max(sif)], peak = max(sif), range = max(sif) - min(sif)), by = .(crop_class, source)]
print(dcast(pk, crop_class ~ source, value.var = c("peak_doy", "peak", "range")), digits = 3)
# spread among the seven class curves on each day: how far apart the crops are
spread <- cur[source != src_nd, .(sd_between = sd(sif)), by = .(doy, source)][, .(mean_sd_between = mean(sd_between)), by = source]
print(spread)
nd_pk <- cur[source == src_nd, .(ndvi_peak_doy = doy[which.max(sif)], ndvi_apr = max(sif[doy %in% 100:125]) / ndvi_scale,
                                 ndvi_summer = max(sif[doy %in% 170:260]) / ndvi_scale), by = crop_class]
print(nd_pk, digits = 2)

ggsave(file.path(proj_dir, "figS8_footprint_skill.png"), pa, width = 9, height = 5, dpi = 300, bg = "white")
ggsave(file.path(proj_dir, "figS8_footprint_skill.pdf"), pa, width = 9, height = 5, bg = "white")

# seasonal cycles on top (crop-coloured strips), contrasts below (product strips stay grey)
fig <- wrap_elements(full = color_strips(pc, crop_class_colors, classes, ncol = 4)) / pb +
  plot_layout(heights = c(1.25, 0.52)) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(plot.tag = element_text(size = fs_title + 0.5, face = "bold", family = fig_font))
save_fig(fig, "fig4_manuscript", height = 5.7)

# Drought response by growth stage, NDVI and SIF in the same format:
# active minus inactive crops (residual from each crop's climatology) by DSCI
# class. Active = the crop's climatology above the midpoint of its annual range.
#   a  HLS NDVI, monthly means of per-scene contrasts (drought_by_stage_ndvi.R)
#   b  HLS NDVI on OCO-2 target days only: the nearest clear scene to each
#      target, with the target's DSCI (drought_stage_ndvi_on_targets.R)
#   c  OCO-2 crop-class SIF, one value per target (drought_by_stage.R)
# All panels use the same growth-stage labels, from the HLS NDVI calendar
# (ndvi_stage_calendar.csv). A supplementary version repeats c with stages
# defined from the SIF climatology (stage_contrast_by_target_sifstage.csv).
# a -> b separates NDVI's sampling advantage from its spatial (30-m) advantage.
# Output: drought_by_stage/plot_stage_contrast_3panel.png / .pdf

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
dir_s <- file.path(proj_dir, "drought_by_stage")
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
nd <- fread(file.path(dir_s, "ndvi_stage_contrast_by_month.csv"))
sf <- fread(file.path(dir_s, "stage_contrast_by_target.csv"))
nt <- fread(file.path(dir_s, "ndvi_stage_contrast_on_targets.csv"))[!is.na(d_abs)]
for (d in list(nd, sf, nt)) d[, `:=`(drought = factor(drought, levels = lev),
                                 period = ifelse(year %in% 2022:2023, "2022-2023 drought", "Other years"))]
pal <- setNames(c("#ffffb2", "#fecc5c", "#e31a1c"), lev)

stat_label <- function(y, dsci, unit) {
  t <- suppressWarnings(cor.test(y, dsci, method = "spearman"))
  # 95% interval for r from the Fisher transform (approximate for Spearman)
  z <- atanh(t$estimate); se <- 1 / sqrt(length(y) - 3)
  sprintf("Spearman r = %.2f (95%% CI %.2f to %.2f)\np %s, n = %d %s",
          t$estimate, tanh(z - 1.96 * se), tanh(z + 1.96 * se),
          ifelse(t$p.value < 0.001, "< 0.001", sprintf("= %.2f", t$p.value)), length(y), unit)
}
panel <- function(d, y, ylab, title, unit) {
  n_by <- d[, .N, by = drought]
  ggplot(d, aes(drought, .data[[y]])) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.6) +
    geom_point(aes(shape = period), position = position_jitter(width = 0.12, seed = 1),
               size = 1.5, alpha = 0.75, fill = "grey15", color = "grey15") +
    annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.3, size = 3, label = stat_label(d[[y]], d$dsci, unit)) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c("2022-2023 drought" = 24, "Other years" = 21), name = NULL) +
    scale_x_discrete(labels = function(x) {
      i <- match(x, n_by$drought); paste0(sub(" \\(", "\n(", x), "\nn = ", n_by$N[i]) }) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
}
pa <- panel(nd, "d_abs", "Active minus inactive crops (NDVI residual)",
            "HLS NDVI, all clear scenes (monthly means)", "months")
pt <- panel(nt, "d_abs", "Active minus inactive crops (NDVI residual)",
            "HLS NDVI on OCO-2 target days (nearest scene)", "targets")
pb <- panel(sf, "d_abs", expression(atop("Active minus inactive crops", (SIF[757]~residual*","~mW~m^-2~sr^-1~nm^-1))),
            "OCO-2 crop-class SIF", "targets")
sf_alt <- fread(file.path(dir_s, "stage_contrast_by_target_sifstage.csv"))
sf_alt[, `:=`(drought = factor(drought, levels = lev), period = ifelse(year %in% 2022:2023, "2022-2023 drought", "Other years"))]
pb_alt <- panel(sf_alt, "d_abs", expression(atop("Active minus inactive crops", (SIF[757]~residual*","~mW~m^-2~sr^-1~nm^-1))),
                "OCO-2 crop-class SIF, stages from the SIF climatology", "targets")
# a and b are both NDVI: same y scale so the boxes compare directly
ylim_nd <- range(c(nd$d_abs, nt$d_abs))
ylim_nd <- ylim_nd + c(-0.04, 0.16) * diff(ylim_nd)
pa <- pa + coord_cartesian(ylim = ylim_nd); pt <- pt + coord_cartesian(ylim = ylim_nd)
g <- pa + pt + pb + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 12, face = "bold"))
ggsave(file.path(dir_s, "plot_stage_contrast_3panel.png"), g, width = 14.5, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(dir_s, "plot_stage_contrast_3panel.pdf"), g, width = 14.5, height = 5.4, bg = "white")
g_alt <- pb + labs(title = "Stages from the NDVI calendar (main)") + pb_alt + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 12, face = "bold"))
ggsave(file.path(dir_s, "plot_stage_contrast_sif_stage_definitions.png"), g_alt, width = 10, height = 5.4, dpi = 300, bg = "white")

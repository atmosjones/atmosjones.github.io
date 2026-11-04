# Figure 5 for the manuscript: drought, greenness and growth stage.
#   a  target-mean SIF residual by DSCI class, all targets        (drought_analysis.R)
#   b  the same by season                                          (drought_analysis.R)
#   c  SIF residual against HLS NDVI residual on OCO-2 days        (compare_ndvi.R output)
#   d  NDVI, active minus inactive crops by DSCI class, all scenes (drought_by_stage_usda.R)
#   e  SIF, the same on OCO-2 targets                               (drought_by_stage_usda.R)
# Growth stages in d-e from USDA NASS Crop Progress (crop_progress_calendar.R).
# Run drought_analysis.R and drought_by_stage_usda.R first.
# Output: fig5_manuscript.png / .pdf

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
p1 <- readRDS(file.path(proj_dir, "fig5_panels_drought.rds"))
p2 <- readRDS(file.path(proj_dir, "fig5_panels_stage.rds"))
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
pal <- setNames(c("#ffffb2", "#fecc5c", "#e31a1c"), lev)

d <- fread(file.path(proj_dir, "ndvi_comparison", "sif_vs_ndvi_residuals_by_day.csv"))
d[, drought := factor(fcase(dsci < cls$low_upper, lev[1], dsci < cls$high_lower, lev[2], default = lev[3]), levels = lev)]
t <- suppressWarnings(cor.test(d$sif_residual, d$ndvi_residual, method = "spearman"))
pc <- ggplot(d, aes(ndvi_residual, sif_residual)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey75", linewidth = 0.6) +
  geom_point(aes(fill = drought), shape = 21, size = 2.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 3,
           label = sprintf("Spearman r = %.2f (p %s)\nn = %d targets", t$estimate,
                           ifelse(t$p.value < 0.001, "< 0.001", sprintf("= %.3f", t$p.value)), nrow(d))) +
  scale_fill_manual(values = pal, name = "DSCI class") +
  labs(x = "NDVI residual (mean over crops)", y = expression(SIF[757]~residual~(mean~over~crops)),
       title = "SIF and NDVI departures, OCO-2 days") +
  theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))

small <- theme(plot.title = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 9))
y_top <- coord_cartesian(ylim = c(-0.27, 0.53))    # same scale in a and b, room for the statistics
top <- (p1$a + y_top + small | p1$b + y_top + small) + plot_layout(widths = c(1, 2.2))
noshape <- guides(shape = "none")
bottom <- (pc | p2$ndvi_all + labs(title = "NDVI: active minus inactive crops, all scenes") + small + noshape |
           p2$sif + labs(title = "SIF: active minus inactive crops, OCO-2 targets") + small + noshape) +
  plot_layout(widths = c(1, 1, 1))
fig <- top / bottom + plot_layout(heights = c(1, 1.05), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 12, face = "bold"))
ggsave(file.path(proj_dir, "fig5_manuscript.png"), fig, width = 14, height = 10.5, dpi = 300, bg = "white")
ggsave(file.path(proj_dir, "fig5_manuscript.pdf"), fig, width = 14, height = 10.5, bg = "white")
cat(sprintf("c: r = %.2f, p = %.2g, n = %d\n", t$estimate, t$p.value, nrow(d)))

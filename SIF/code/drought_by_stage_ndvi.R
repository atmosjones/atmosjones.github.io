# Drought response by crop and growth stage in HLS NDVI (compare_ndvi.R),
# the same tests as drought_by_stage.R for SIF:
#   stage    active where the crop's NDVI climatology (ndvi_fit) is above the
#            midpoint between its annual minimum and maximum
#   d        mean residual of active crops minus inactive crops in each scene
#   crops    each crop's residual relative to the scene mean
# Scenes a few days apart share the same drought state, so every test is also
# run on monthly means (year-month), a more conservative sample size.
# Output (drought_by_stage/): ndvi_stats.csv, ndvi_stage_contrast_by_month.csv,
#   ndvi_stage_calendar.csv (active or not, per crop and day of year),
#   plot_ndvi_drought_by_crop.png

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "drought_by_stage")
dir.create(out_dir, showWarnings = FALSE)
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
n <- fread(file.path(proj_dir, "ndvi_comparison", "ndvi_residuals.csv"))[, date := as.Date(date)]
n[, drought := factor(drought, levels = lev)]
n <- merge(n, n[, .(fit_min = min(ndvi_fit), fit_max = max(ndvi_fit)), by = crop_class], by = "crop_class")
n[, stage := fifelse(ndvi_fit >= (fit_min + fit_max) / 2, "Active", "Inactive")]
# the same rule as a calendar (crop x day of year) from the climatology curve, for
# labelling other data (drought_by_stage.R uses it for the SIF test)
calendar <- n[, {
  u <- unique(.SD[, .(doy, ndvi_fit)])[order(doy)]
  u <- rbind(copy(u)[, doy := doy - 365], u, copy(u)[, doy := doy + 365])
  curve <- approx(u$doy, u$ndvi_fit, xout = 1:366, ties = mean, rule = 2)$y
  .(doy = 1:366, active = curve >= (fit_min[1] + fit_max[1]) / 2)
}, by = crop_class]
fwrite(calendar, file.path(out_dir, "ndvi_stage_calendar.csv"))
n[, rel := residual / ndvi_fit]
n[, rel_to_scene := residual - mean(residual), by = date]
n[, ym := format(date, "%Y-%m")]
cat(sprintf("%d scenes, %d class-scenes, %d months\n", uniqueN(n$date), nrow(n), uniqueN(n$ym)))
print(dcast(n[, .(active = round(mean(stage == "Active"), 2)), by = .(crop_class, season)], crop_class ~ season,
            value.var = "active"))

sp <- function(x, y) { t <- suppressWarnings(cor.test(x, y, method = "spearman")); c(r = unname(t$estimate), p = t$p.value) }
row <- function(label, x, dsci, grp) {
  a <- sp(x, dsci)
  data.table(test = label, n = length(x), spearman_r = a[["r"]], spearman_p = a[["p"]],
             kruskal_p = if (nlevels(droplevels(grp)) > 1) kruskal.test(x ~ droplevels(grp))$p.value else NA_real_,
             med_none = median(x[grp == lev[1]]), med_dry = median(x[grp == lev[2]]), med_drought = median(x[grp == lev[3]]))
}
# stage contrast per scene, and per month
d <- n[, .(d_abs = mean(residual[stage == "Active"]) - mean(residual[stage == "Inactive"]),
           d_rel = mean(rel[stage == "Active"]) - mean(rel[stage == "Inactive"]),
           ok = any(stage == "Active") & any(stage == "Inactive"),
           dsci = dsci[1], drought = drought[1], season = season[1], ym = ym[1], year = year[1]), by = date][ok == TRUE]
dm <- d[, .(d_abs = mean(d_abs), d_rel = mean(d_rel), dsci = mean(dsci), season = season[1], year = year[1]), by = ym]
dm[, drought := factor(fcase(dsci < cls$low_upper, lev[1], dsci < cls$high_lower, lev[2], default = lev[3]), levels = lev)]
fwrite(dm, file.path(out_dir, "ndvi_stage_contrast_by_month.csv"))
res <- rbind(
  d[, row("Stage: active minus inactive, absolute (scenes)", d_abs, dsci, drought)],
  d[, row("Stage: active minus inactive, relative (scenes)", d_rel, dsci, drought)],
  dm[, row("Stage: active minus inactive, absolute (months)", d_abs, dsci, drought)],
  dm[, row("Stage: active minus inactive, relative (months)", d_rel, dsci, drought)],
  dm[!year %in% 2022:2023, row("Stage: absolute (months), without 2022-2023", d_abs, dsci, drought)],
  d[, row(paste("Stage, absolute (scenes):", season[1]), d_abs, dsci, drought), by = season][, season := NULL])
# each crop relative to the other crops in the scene, per scene and per month
cm <- n[, .(rel_to_scene = mean(rel_to_scene), dsci = mean(dsci)), by = .(crop_class, ym)]
cm[, drought := factor(fcase(dsci < cls$low_upper, lev[1], dsci < cls$high_lower, lev[2], default = lev[3]), levels = lev)]
res <- rbind(res,
  n[, row(paste("Relative to other crops (scenes):", crop_class[1]), rel_to_scene, dsci, drought), by = crop_class][, crop_class := NULL],
  cm[, row(paste("Relative to other crops (months):", crop_class[1]), rel_to_scene, dsci, drought), by = crop_class][, crop_class := NULL])
fwrite(res, file.path(out_dir, "ndvi_stats.csv"))
print(res, digits = 2)

# ---- figure: each crop relative to the scene mean, by DSCI class (monthly means)
crop_class_colors <- c("Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
  "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e")
cm[, crop_class := factor(crop_class, levels = names(crop_class_colors))]
pal <- setNames(c("#ffffb2", "#fecc5c", "#e31a1c"), lev)
pa <- ggplot(cm, aes(crop_class, rel_to_scene, fill = drought)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_boxplot(outlier.size = 0.5, width = 0.75) +
  scale_fill_manual(values = pal, name = "DSCI class") +
  scale_x_discrete(labels = function(x) sub("Dbl Crop WinWht/Soybeans", "Dbl crop\nWW/Soy", sub("/", "/\n", x))) +
  labs(x = NULL, y = "NDVI residual minus scene mean (monthly means)",
       title = "Each crop relative to the other crops, by drought class") +
  theme_bw() + theme(legend.position = "bottom", plot.title = element_text(size = 10, face = "bold"))
pb <- ggplot(dm, aes(drought, d_abs)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.6) +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = c("None", "Dry", "Drought")) +
  labs(x = NULL, y = "Active minus inactive crops (NDVI residual)", title = "Growth stage (monthly means)") +
  theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
ggsave(file.path(out_dir, "plot_ndvi_drought_by_crop.png"), pa + pb + plot_layout(widths = c(2.3, 1)),
       width = 12, height = 5, dpi = 300, bg = "white")

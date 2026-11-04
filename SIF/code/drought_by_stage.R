# Do crops respond to drought differently depending on growth stage?
#
# Growth stage (main): from the HLS NDVI calendar (drought_by_stage_ndvi.R,
# ndvi_stage_calendar.csv): a crop is active on days when its NDVI climatology is
# above the midpoint of its annual range, i.e. when its fields are actually green.
# NDVI is independent of the SIF being tested, and gives all panels of
# plot_stage_contrast.R the same stage labels. The SIF climatology is not used for
# this because spring wheat mixing lifts the summer crops' SIF curves and marks
# soybeans active in April-May and corn in May, when their fields are bare.
# Alternative (reported, *_sifstage): the same rule on each crop's SIF climatology.
#   inactive otherwise (dormant, not yet planted, or harvested)
# Crops on the same target share a landscape-wide anomaly, so stage effects are
# measured within each target: d = mean residual of active crops minus mean
# residual of inactive crops, on targets with both present. One value per
# target, tested against DSCI (Spearman) and across the DSCI classes
# (Kruskal-Wallis). Repeated with residuals relative to the seasonal curve
# (residual / sif_fit), since active crops have more SIF to lose.
# Also: each crop's residual relative to the target mean, by DSCI class.
# Output (drought_by_stage/): stage_contrast_by_target.csv (NDVI stages),
#   stage_contrast_by_target_sifstage.csv, stats.csv (both),
#   plot_drought_by_stage.png

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "drought_by_stage")
dir.create(out_dir, showWarnings = FALSE)
s <- fread(file.path(proj_dir, "drought_residuals.csv"))[, date := as.Date(date)]
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
s[, drought := factor(drought, levels = unlist(cls[, .(label_low, label_mid, label_high)]))]

# climatology range per crop from the fitted values over the year
rng <- s[, .(fit_min = min(sif_fit), fit_max = max(sif_fit)), by = crop_class]
s <- merge(s, rng, by = "crop_class")
s[, stage_sif := fifelse(sif_fit >= (fit_min + fit_max) / 2, "Active", "Inactive")]
cal <- fread(file.path(out_dir, "ndvi_stage_calendar.csv"))
s <- merge(s, cal, by = c("crop_class", "doy"), all.x = TRUE)
if (s[, anyNA(active)]) stop("crop or day missing from the NDVI calendar")
s[, stage_ndvi := fifelse(active, "Active", "Inactive")]
s[, stage := stage_ndvi]
cat("stage agreement, NDVI vs SIF definition (class-days):\n"); print(s[, .N, by = .(stage_ndvi, stage_sif)])
s[, rel := residual / pmax(sif_fit, 0.1)]                      # fraction of the expected SIF
cat("class-days by season and stage:\n"); print(dcast(s, season ~ stage, fun.aggregate = length, value.var = "rel"))
cat("\nactive share by crop and season:\n")
print(dcast(s[, .(active = round(mean(stage == "Active"), 2)), by = .(crop_class, season)], crop_class ~ season))

contrast <- function(st) {
  s[, .(n_active = sum(get(st) == "Active"), n_inactive = sum(get(st) == "Inactive"),
        d_abs = mean(residual[get(st) == "Active"]) - mean(residual[get(st) == "Inactive"]),
        d_rel = mean(rel[get(st) == "Active"]) - mean(rel[get(st) == "Inactive"]),
        dsci = dsci[1], drought = drought[1], season = season[1], year = year[1]), by = date][
    n_active > 0 & n_inactive > 0]
}
d <- contrast("stage_ndvi")
d_sifstage <- contrast("stage_sif")
fwrite(d, file.path(out_dir, "stage_contrast_by_target.csv"))
fwrite(d_sifstage, file.path(out_dir, "stage_contrast_by_target_sifstage.csv"))

test <- function(x, dsci, grp, label) {
  t <- suppressWarnings(cor.test(x, dsci, method = "spearman"))
  data.table(test = label, n = length(x), spearman_r = unname(t$estimate), spearman_p = t$p.value,
             kruskal_p = kruskal.test(x ~ grp)$p.value,
             median_none = median(x[as.integer(grp) == 1]), median_dry = median(x[as.integer(grp) == 2]),
             median_drought = median(x[as.integer(grp) == 3]))
}
st <- rbind(
  d[, test(d_abs, dsci, drought, "Active minus inactive, absolute residual")],
  d[, test(d_rel, dsci, drought, "Active minus inactive, relative residual")],
  d[!year %in% 2022:2023, test(d_abs, dsci, drought, "Absolute, without 2022-2023")],
  d[!year %in% 2022:2023, test(d_rel, dsci, drought, "Relative, without 2022-2023")],
  d_sifstage[, test(d_abs, dsci, drought, "SIF-defined stages: active minus inactive, absolute")],
  d_sifstage[, test(d_rel, dsci, drought, "SIF-defined stages: active minus inactive, relative")])
# each crop relative to the target mean, by DSCI class
s[, rel_to_day := residual - mean(residual), by = date]
st_crop <- s[, test(rel_to_day, dsci, drought, paste("Relative to other crops:", crop_class[1])), by = crop_class][, crop_class := NULL]
st_season <- d[, test(d_abs, dsci, drought, paste("Active minus inactive, absolute:", season[1])), by = season][, season := NULL]
res <- rbind(st, st_season, st_crop, fill = TRUE)
fwrite(res, file.path(out_dir, "stats.csv"))
cat("\n"); print(res, digits = 2)

# ---- figure -------------------------------------------------------------------
pal <- setNames(c("#ffffb2", "#fecc5c", "#e31a1c"), levels(s$drought))
box <- function(dd, y, ylab, title, lab) {
  ggplot(dd, aes(drought, .data[[y]])) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.3, alpha = 0.7) +
    annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.4, size = 3, label = lab) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
    labs(x = NULL, y = ylab, title = title) + theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold"))
}
lab_of <- function(r) sprintf("Spearman r = %.2f (p = %.2f)\nKruskal-Wallis p = %.2f\nn = %d targets",
                              r$spearman_r, r$spearman_p, r$kruskal_p, r$n)
pa <- box(d, "d_abs", "Active minus inactive crops (residual SIF)", "Absolute residuals", lab_of(st[1]))
pb <- box(d, "d_rel", "Active minus inactive crops (residual / expected SIF)", "Relative residuals", lab_of(st[2]))
ggsave(file.path(out_dir, "plot_drought_by_stage.png"), pa + pb, width = 10, height = 4.8, dpi = 300, bg = "white")

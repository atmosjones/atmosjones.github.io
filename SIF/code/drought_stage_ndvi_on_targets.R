# Is NDVI's growth-stage drought signal (drought_by_stage_ndvi.R) a matter of
# spatial resolution (pure 30-m crop pixels) or of sampling (427 scenes vs 104
# OCO-2 targets)? Repeat the NDVI test using only OCO-2 target days: for each
# target, the nearest clear HLS scene (any gap; gap reported, and a <= 5-day
# subset tested). DSCI is the target's own value, as in the SIF test.
#   set "all targets"   every target whose matched scene has active and inactive crops
#   set "SIF targets"   the targets of the SIF test (drought_by_stage.R)
# Stage and residuals as in drought_by_stage_ndvi.R.
# Output (drought_by_stage/): ndvi_stage_contrast_on_targets.csv, printed tests.

library(data.table)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
dir_s <- file.path(proj_dir, "drought_by_stage")
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
n <- fread(file.path(proj_dir, "ndvi_comparison", "ndvi_residuals.csv"))[, date := as.Date(date)]
n <- merge(n, n[, .(fit_min = min(ndvi_fit), fit_max = max(ndvi_fit)), by = crop_class], by = "crop_class")
n[, stage := fifelse(ndvi_fit >= (fit_min + fit_max) / 2, "Active", "Inactive")]
scene <- n[, .(d_abs = if (any(stage == "Active") && any(stage == "Inactive"))
                 mean(residual[stage == "Active"]) - mean(residual[stage == "Inactive"]) else NA_real_), by = date]

tg <- unique(fread(file.path(proj_dir, "drought_residuals.csv"))[, .(date = as.Date(date), dsci, season, year)])
sif_t <- fread(file.path(dir_s, "stage_contrast_by_target.csv"))[, date := as.Date(date)]
sc <- sort(unique(scene$date))
tg[, hls_date := sc[vapply(date, function(d) which.min(abs(as.numeric(sc - d))), 1L)]]
tg[, gap_days := as.integer(hls_date - date)]
tg <- merge(tg, scene, by.x = "hls_date", by.y = "date")
tg[, drought := factor(fcase(dsci < cls$low_upper, lev[1], dsci < cls$high_lower, lev[2], default = lev[3]), levels = lev)]
tg[, in_sif_set := date %in% sif_t$date]
fwrite(tg, file.path(dir_s, "ndvi_stage_contrast_on_targets.csv"))

cat(sprintf("%d targets; nearest HLS scene: median gap %d days, %d within 5 days, %d within 16 days, max %d\n",
            nrow(tg), as.integer(median(abs(tg$gap_days))), tg[abs(gap_days) <= 5, .N], tg[abs(gap_days) <= 16, .N],
            max(abs(tg$gap_days))))
cat(sprintf("distinct scenes used: %d\n", uniqueN(tg$hls_date)))
test <- function(d, label) {
  d <- d[!is.na(d_abs)]
  t <- suppressWarnings(cor.test(d$d_abs, d$dsci, method = "spearman"))
  z <- atanh(t$estimate); se <- 1 / sqrt(nrow(d) - 3)
  data.table(set = label, n = nrow(d), r = unname(t$estimate), lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
             p = t$p.value, kw_p = kruskal.test(d_abs ~ drought, d)$p.value,
             med_none = median(d[drought == lev[1], d_abs]), med_dry = median(d[drought == lev[2], d_abs]),
             med_drought = median(d[drought == lev[3], d_abs]))
}
res <- rbind(test(tg, "NDVI, all OCO-2 targets (nearest scene)"),
             test(tg[abs(gap_days) <= 5], "NDVI, all targets, scene within 5 days"),
             test(tg[in_sif_set == TRUE], "NDVI, the SIF-test targets"),
             test(tg[in_sif_set == TRUE & abs(gap_days) <= 5], "NDVI, the SIF-test targets, within 5 days"),
             test(tg[!year %in% 2022:2023], "NDVI, all targets, without 2022-2023"),
             test(sif_t[, .(d_abs, dsci, drought = factor(drought, levels = lev))], "SIF, NDVI-defined stages (reference)"))
fwrite(res, file.path(dir_s, "ndvi_on_targets_stats.csv"))
print(res, digits = 2)

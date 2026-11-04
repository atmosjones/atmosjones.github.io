# Does the drought result survive uncertainty in the seasonal climatology?
# In each bootstrap replicate, targets are resampled with replacement and every
# crop's seasonal curve is refitted exactly as in drought_analysis.R (leave-year-out:
# the curve used for a year is fitted without that year). The anomalies of the
# ORIGINAL 104 targets are then recomputed against the replicate's curves and the
# main drought tests are rerun on them:
#   target-mean anomaly vs DSCI (Spearman), across DSCI classes (Kruskal-Wallis),
#   class medians, autumn and winter correlations, without 2022-2023, and the part
#   of the SIF anomaly not explained by the NDVI anomaly vs DSCI.
# Replicate 0 is the unresampled fit and must reproduce drought_analysis.R.
# Output (bootstrap_climatology/): drought_tests_by_replicate.csv, drought_tests_summary.csv,
#   plot_drought_tests_bootstrap.png

library(data.table)
library(ggplot2)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "bootstrap_climatology")
n_boot <- 1000; loess_span <- 0.4; set.seed(7)
s <- fread(file.path(proj_dir, "drought_residuals.csv"))[, date := as.Date(date)]
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
nd <- fread(file.path(proj_dir, "ndvi_comparison", "sif_vs_ndvi_residuals_by_day.csv"))[, date := as.Date(date)]

seasonal_fit <- function(fit_data, new_doy) {                 # as in drought_analysis.R
  wrap <- rbind(fit_data, copy(fit_data)[, doy := doy - 365], copy(fit_data)[, doy := doy + 365])
  m <- loess(sif757 ~ doy, data = wrap, weights = 1 / sif757_sd_scaled^2, span = loess_span / 3, degree = 2)
  predict(m, newdata = data.frame(doy = new_doy))
}
anomalies <- function(train_all) {                            # curves from train_all, applied to s
  out <- s[, .(date, crop_class, year, doy, sif757)]
  out[, fit := {
    v <- numeric(.N)
    tr <- train_all[crop_class == .BY$crop_class]
    for (y in unique(year)) {
      rows <- year == y
      v[rows] <- tryCatch(seasonal_fit(tr[year != y], doy[rows]), error = function(e) NA_real_)
    }
    v
  }, by = crop_class]
  out[, residual := sif757 - fit]
}
tests <- function(a) {
  t <- a[, .(residual = mean(residual, na.rm = TRUE)), by = date]
  t <- merge(t, unique(s[, .(date, dsci, drought, season, year)]), by = "date")
  t[, drought := factor(drought, levels = lev)]
  sp <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman"))
  r_all <- sp(t$residual, t$dsci)
  r_no <- t[!year %in% 2022:2023, sp(residual, dsci)]
  m <- merge(t, nd[, .(date, ndvi_residual)], by = "date")
  m[, beyond := resid(lm(residual ~ ndvi_residual))]
  data.table(r = unname(r_all$estimate), p = r_all$p.value,
             kw_p = kruskal.test(residual ~ drought, t)$p.value,
             med_none = t[drought == lev[1], median(residual)], med_dry = t[drought == lev[2], median(residual)],
             med_drought = t[drought == lev[3], median(residual)],
             r_autumn = t[season == "Autumn", unname(sp(residual, dsci)$estimate)],
             p_autumn = t[season == "Autumn", sp(residual, dsci)$p.value],
             r_winter = t[season == "Winter", unname(sp(residual, dsci)$estimate)],
             p_winter = t[season == "Winter", sp(residual, dsci)$p.value],
             r_no2223 = unname(r_no$estimate), p_no2223 = r_no$p.value,
             r_sif_ndvi = cor(m$residual, m$ndvi_residual, method = "spearman"),
             r_beyond_ndvi = unname(sp(m$beyond, m$dsci)$estimate), p_beyond_ndvi = sp(m$beyond, m$dsci)$p.value)
}
train0 <- s[, .(crop_class, year, doy, sif757, sif757_sd_scaled, date)]
base <- tests(anomalies(train0))[, rep := 0L]
cat("replicate 0 (should match drought_analysis.R: r = -0.26, p = 0.007, KW p = 0.018):\n"); print(base, digits = 2)
if (abs(base$r - (-0.26)) > 0.01 || abs(base$kw_p - 0.018) > 0.002) stop("replicate 0 does not reproduce drought_analysis.R")

dates <- unique(s$date)
reps <- rbindlist(lapply(seq_len(n_boot), function(b) {
  d <- train0[data.table(date = sample(dates, replace = TRUE)), on = "date", allow.cartesian = TRUE]
  tests(anomalies(d))[, rep := b]
}))
fwrite(rbind(base, reps), file.path(out_dir, "drought_tests_by_replicate.csv"))

q <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, .025), quantile(x, .975))
summ <- data.table(
  test = c("Spearman r, target-mean anomaly vs DSCI", "  p", "  share of replicates with p < 0.05",
           "Kruskal-Wallis p across DSCI classes", "  share with p < 0.05",
           "median anomaly: none", "median anomaly: dry", "median anomaly: drought", "drought minus none",
           "autumn r", "  share with p < 0.05", "winter r", "  share with p < 0.05",
           "r without 2022-2023", "  share with p < 0.05",
           "SIF vs NDVI anomaly r", "SIF beyond NDVI vs DSCI r", "  share with p < 0.05"),
  published = c(base$r, base$p, NA, base$kw_p, NA, base$med_none, base$med_dry, base$med_drought,
                base$med_drought - base$med_none, base$r_autumn, NA, base$r_winter, NA, base$r_no2223, NA,
                base$r_sif_ndvi, base$r_beyond_ndvi, NA),
  bootstrap_median_95 = c(q(reps$r), q(reps$p), sprintf("%.2f", mean(reps$p < 0.05)), q(reps$kw_p),
                          sprintf("%.2f", mean(reps$kw_p < 0.05)), q(reps$med_none), q(reps$med_dry),
                          q(reps$med_drought), q(reps$med_drought - reps$med_none), q(reps$r_autumn),
                          sprintf("%.2f", mean(reps$p_autumn < 0.05)), q(reps$r_winter),
                          sprintf("%.2f", mean(reps$p_winter < 0.05)), q(reps$r_no2223),
                          sprintf("%.2f", mean(reps$p_no2223 < 0.05)), q(reps$r_sif_ndvi),
                          q(reps$r_beyond_ndvi), sprintf("%.2f", mean(reps$p_beyond_ndvi < 0.05))))
fwrite(summ, file.path(out_dir, "drought_tests_summary.csv"))
print(summ, digits = 3)

long <- melt(reps[, .(`Spearman r (all)` = r, `Autumn r` = r_autumn, `Winter r` = r_winter,
                      `Without 2022-2023 r` = r_no2223, `Drought minus none (median)` = med_drought - med_none)],
             measure.vars = 1:5)
base_long <- data.table(variable = levels(long$variable),
                        value = c(base$r, base$r_autumn, base$r_winter, base$r_no2223, base$med_drought - base$med_none))
p <- ggplot(long, aes(value)) +
  geom_histogram(bins = 40, fill = "grey60", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_vline(data = base_long, aes(xintercept = value), color = "#d95f02", linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  labs(x = NULL, y = "Bootstrap replicates",
       title = sprintf("Drought tests with the seasonal climatology refitted in %d bootstrap replicates (orange: published value)", n_boot)) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9))
ggsave(file.path(out_dir, "plot_drought_tests_bootstrap.png"), p, width = 13, height = 3.2, dpi = 300, bg = "white")

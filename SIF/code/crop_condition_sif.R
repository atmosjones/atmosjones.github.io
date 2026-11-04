# Does crop-class SIF track USDA crop condition? (design fixed before running)
#   condition  weekly % of the Oklahoma crop rated good or excellent (NASS Crop
#              Progress, fetch_crop_progress.R), minus that crop's 2014-2025 mean
#   matching   each crop-class SIF anomaly (drought_residuals.csv) in its USDA
#              active season (crop_progress_calendar.R) is paired with the nearest
#              condition report for that crop within 7 days; double crop uses the
#              wheat rating before the wheat harvest and the soybean rating after
#   test 1     one value per target: mean SIF anomaly of the rated active crops vs
#              mean condition anomaly of the same crops (Spearman)
#   test 2     crop-specific, targets with >= 2 rated crops: each crop's SIF anomaly
#              minus the others' mean vs its condition anomaly minus the others' mean
#              (Spearman; 95% CI from 2,000 bootstrap resamples of targets)
# Context: correlation of the condition anomaly with the DSCI.
# Output (crop_progress/): condition_matched.csv, condition_stats.csv,
#   plot_sif_vs_condition.png

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
cp_dir <- file.path(proj_dir, "crop_progress")
source(file.path(proj_dir, "crop_progress_calendar.R"))      # windows, active_usda()
set.seed(1)

cond <- fread(file.path(cp_dir, "ok_crop_condition_raw.csv"))
cond[, `:=`(date = as.Date(week_ending), pct = as.numeric(gsub(",", "", Value)),
            rating = sub(".*MEASURED IN PCT ", "", short_desc))]
ge <- cond[rating %in% c("GOOD", "EXCELLENT") & !is.na(pct), .(ge = sum(pct)), by = .(commodity = commodity_desc, date)]
ge[, ge_anom := ge - mean(ge), by = commodity]

s <- fread(file.path(proj_dir, "drought_residuals.csv"))[, date := as.Date(date)][year %between% range(windows$year)]
s[, active := active_usda(crop_class, date)]
s <- s[active == TRUE]
wheat_harvest <- windows[crop_class == "Winter Wheat", .(year, wh = end)]
s <- merge(s, wheat_harvest, by = "year")
s[, commodity := fcase(crop_class == "Winter Wheat", "WHEAT", crop_class == "Corn", "CORN",
                       crop_class == "Sorghum", "SORGHUM", crop_class == "Soybeans", "SOYBEANS",
                       crop_class == "Cotton", "COTTON",
                       crop_class == "Dbl Crop WinWht/Soybeans", ifelse(date <= wh, "WHEAT", "SOYBEANS"))]
match_cond <- function(com, d) {
  g <- ge[commodity == com]
  if (!nrow(g)) return(c(NA_real_, NA_real_))
  i <- which.min(abs(as.numeric(g$date - d)))
  if (abs(as.numeric(g$date[i] - d)) > 7) c(NA_real_, NA_real_) else c(g$ge[i], g$ge_anom[i])
}
m <- t(mapply(match_cond, s$commodity, s$date))
s[, `:=`(ge = m[, 1], ge_anom = m[, 2])]
s <- s[!is.na(ge_anom)]
fwrite(s[, .(date, crop_class, commodity, residual, ge, ge_anom, dsci, season, year)], file.path(cp_dir, "condition_matched.csv"))
cat(sprintf("matched crop-targets: %d on %d targets\n", nrow(s), uniqueN(s$date)))
print(s[, .N, by = crop_class])

sp <- function(x, y) { t <- suppressWarnings(cor.test(x, y, method = "spearman")); c(unname(t$estimate), t$p.value) }
fz <- function(r, n) { z <- atanh(r); se <- 1 / sqrt(n - 3); tanh(z + c(-1.96, 1.96) * se) }
# test 1
t1 <- s[, .(sif = mean(residual), cond = mean(ge_anom), dsci = dsci[1], year = year[1]), by = date]
a1 <- sp(t1$sif, t1$cond); ci1 <- fz(a1[1], nrow(t1))
a1b <- t1[!year %in% 2022:2023, sp(sif, cond)]
# test 2
s[, n_day := .N, by = date]
t2 <- s[n_day >= 2][, `:=`(sif_rel = residual - (sum(residual) - residual) / (.N - 1),
                           cond_rel = ge_anom - (sum(ge_anom) - ge_anom) / (.N - 1)), by = date]
a2 <- sp(t2$sif_rel, t2$cond_rel)
days <- unique(t2$date)
boot <- replicate(2000, {
  b <- data.table(date = sample(days, replace = TRUE))
  bb <- t2[b, on = "date", allow.cartesian = TRUE]
  suppressWarnings(cor(bb$sif_rel, bb$cond_rel, method = "spearman"))
})
st <- data.table(
  test = c("1. target mean: SIF anomaly vs condition anomaly", "1. without 2022-2023",
           "2. crop-specific (relative to other crops)", "context: condition anomaly vs DSCI (targets)"),
  n = c(nrow(t1), t1[!year %in% 2022:2023, .N], nrow(t2), nrow(t1)),
  spearman_r = c(a1[1], a1b[1], a2[1], sp(t1$cond, t1$dsci)[1]),
  ci_lo = c(ci1[1], NA, quantile(boot, .025, na.rm = TRUE), NA),
  ci_hi = c(ci1[2], NA, quantile(boot, .975, na.rm = TRUE), NA),
  p = c(a1[2], a1b[2], a2[2], sp(t1$cond, t1$dsci)[2]))
fwrite(st, file.path(cp_dir, "condition_stats.csv"))
cat(sprintf("test 2: %d crop-targets on %d targets\n", nrow(t2), length(days)))
print(st, digits = 2)

# ---- figure -----------------------------------------------------------------------
crop_class_colors <- c("Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
                       "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3")
lab <- function(r, lo, hi, p, n, unit) sprintf("Spearman r = %.2f (95%% CI %.2f to %.2f)\np = %.2f, n = %d %s", r, lo, hi, p, n, unit)
pa <- ggplot(t1, aes(cond, sif)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", fill = "grey75", linewidth = 0.6) +
  geom_point(aes(shape = ifelse(year %in% 2022:2023, "2022-2023", "Other years")), size = 2, fill = "grey30") +
  scale_shape_manual(values = c("2022-2023" = 24, "Other years" = 21), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3,
           label = lab(st$spearman_r[1], st$ci_lo[1], st$ci_hi[1], st$p[1], st$n[1], "targets")) +
  labs(x = "Condition anomaly (% good + excellent, mean of rated crops)", y = "SIF anomaly (mean of rated crops)",
       title = "Target means") + theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
pb <- ggplot(t2, aes(cond_rel, sif_rel, color = crop_class)) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, color = "black", fill = "grey75", linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = crop_class_colors, name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3,
           label = lab(st$spearman_r[3], st$ci_lo[3], st$ci_hi[3], st$p[3], st$n[3], "crop-targets")) +
  labs(x = "Condition anomaly relative to the other crops", y = "SIF anomaly relative to the other crops",
       title = "Crop-specific (relative to the other rated crops that day)") +
  theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
ggsave(file.path(cp_dir, "plot_sif_vs_condition.png"), pa + pb + plot_layout(widths = c(1, 1.25)),
       width = 12, height = 5, dpi = 300, bg = "white")

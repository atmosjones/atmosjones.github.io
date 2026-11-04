# Bootstrap uncertainty of the crop-class seasonal climatology (Figure 3 curves).
# Targets (dates) are resampled with replacement, keeping all crops observed on a
# target together, and each crop's LOESS curve is refitted exactly as in
# make_fig3.R (span 0.4 year, quadratic, weights 1/SD^2, data wrapped a year either
# side). The envelope is the 2.5-97.5% range of the bootstrap curves on each day.
# Also: the share of each crop's estimates lying outside the envelope, and the
# envelope width by season, to judge whether anomalies should be read against it.
# Output (bootstrap_climatology/): envelope.csv, summary.csv, plot_bootstrap_envelopes.png

library(data.table)
library(ggplot2)
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "strip_colors.R"))
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "active_season_bands.R"))

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "bootstrap_climatology")
dir.create(out_dir, showWarnings = FALSE)
classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
             "Sorghum", "Cotton", "Grassland/Pasture")
crop_class_colors <- c("Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
  "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e")
n_boot <- 1000; loess_span <- 0.4; set.seed(42)

seasonal_curve <- function(doy, y, w, new_doy = 1:365) {      # as in make_fig3.R
  d <- data.table(doy, y, w)
  wrap <- rbind(d, copy(d)[, doy := doy - 365], copy(d)[, doy := doy + 365])
  predict(loess(y ~ doy, data = wrap, weights = w, span = loess_span / 3, degree = 2),
          newdata = data.frame(doy = new_doy))
}
cc <- fread(file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv"))[
  crop_class %in% classes & avg_kernel >= 0.9][, date := as.Date(date)]
cc[, `:=`(doy = as.integer(format(date, "%j")), w = 1 / sif757_sd_scaled^2)]
dates <- unique(cc$date)

fit <- cc[, .(doy = 1:365, sif = seasonal_curve(doy, sif757, w)), by = crop_class]
boot <- rbindlist(lapply(seq_len(n_boot), function(b) {
  s <- data.table(date = sample(dates, replace = TRUE))
  d <- cc[s, on = "date", allow.cartesian = TRUE]
  d[, .(doy = 1:365, sif = tryCatch(seasonal_curve(doy, sif757, w), error = function(e) rep(NA_real_, 365))),
    by = crop_class][, rep := b]
}))
env <- boot[, .(lo = quantile(sif, 0.025, na.rm = TRUE), hi = quantile(sif, 0.975, na.rm = TRUE)), by = .(crop_class, doy)]
env <- merge(env, fit, by = c("crop_class", "doy"))
fwrite(env, file.path(out_dir, "envelope.csv"))

# how wide, and how many estimates fall outside (the curve's uncertainty, not the scatter)
season_of <- function(doy) c("Winter", "Spring", "Summer", "Autumn")[c(1,1,2,2,2,3,3,3,4,4,4,1)[as.integer(format(as.Date(doy - 1, origin = "2021-01-01"), "%m"))]]
env[, season := season_of(doy)]
pts <- merge(cc[, .(crop_class, doy, sif757)], env, by = c("crop_class", "doy"))
summ <- merge(env[, .(median_halfwidth = median((hi - lo) / 2), max_halfwidth = max((hi - lo) / 2)), by = crop_class],
              pts[, .(n = .N, share_outside = mean(sif757 < lo | sif757 > hi)), by = crop_class], by = "crop_class")
by_season <- dcast(env[, .(hw = median((hi - lo) / 2)), by = .(crop_class, season)], crop_class ~ season, value.var = "hw")
fwrite(merge(summ, by_season, by = "crop_class"), file.path(out_dir, "summary.csv"))
print(merge(summ, by_season, by = "crop_class"), digits = 2)
# typical residual scatter, for scale
cat(sprintf("median |estimate - curve|: %.2f\n", pts[, median(abs(sif757 - sif))]))

# ---- figure (Figure 3 layout) -------------------------------------------------------
for (d in list(env, cc, pts)) d[, crop_class := factor(crop_class, levels = classes)]
bands <- active_season_bands()[crop_class %in% classes][, crop_class := factor(crop_class, levels = classes)]
p <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey90", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(data = cc, aes(doy, sif757), color = "#1a9850", size = 0.7, alpha = 0.4) +
  geom_ribbon(data = env, aes(doy, ymin = lo, ymax = hi), fill = "#1a9850", alpha = 0.3) +
  geom_line(data = env, aes(doy, sif), color = "#1a9850", linewidth = 0.8) +
  facet_wrap(~crop_class, ncol = 4) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("Jan", "Apr", "Jul", "Oct")) +
  coord_cartesian(ylim = c(-0.25, 1.6)) +
  labs(x = NULL, y = expression(SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1)),
       title = sprintf("Crop-class seasonal cycles with 95%% bootstrap envelopes (%d resamples of %d targets)", n_boot, length(dates))) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8, face = "bold"), plot.title = element_text(face = "bold"))
g <- color_strips(p, crop_class_colors, classes, ncol = 4)
png(file.path(out_dir, "plot_bootstrap_envelopes.png"), width = 10, height = 5.8, units = "in", res = 300)
grid::grid.draw(g); invisible(dev.off())

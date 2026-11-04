# Figure S1: viewing-zenith-angle (VZA) dependence of target-mode SIF757 and the
# correction curve used to normalise it to nadir (vza_correction.R).
#
# Data: the soundings the correction was fitted to (vza_check/soundings_vza.csv; target
# mode, Quality_Flag 0/1, days with error-weighted day-mean SIF > 0.15). Each sounding is
# expressed relative to its day's mean, so season and surface drop out, then scaled by
# the fit's value at nadir, so the y-axis is SIF relative to a nadir view and the
# dashed curve is the correction factor f(VZA) itself.
# Grey: error-weighted mean of each day in each 5-degree bin (>= 20 soundings; single soundings
# scatter too widely to show); black: mean over all soundings in each bin (+/- 2 SE).
# Output: fig_s1_vza.png (600 dpi)

suppressMessages({library(data.table); library(ggplot2)})
set.seed(1)
proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
source(file.path(proj_dir, "fig_fonts.R"))
apply_vza_correction <- TRUE; source(file.path(proj_dir, "vza_correction.R"))

o <- fread(file.path(proj_dir, "vza_check", "soundings_vza.csv"))
o <- o[is.finite(sif) & is.finite(unc) & unc > 0 & is.finite(VZA)]
o[, w := 1 / unc^2]
o[, dm := sum(w * sif) / sum(w), by = date]
o <- o[dm > 0.15]
o[, `:=`(ratio = sif / dm, ww = w * dm^2)]

# the correction's fit: ratio = c * (1 + a1 v + a2 v^2); check it reproduces vza_coef
m <- lm(ratio ~ VZA + I(VZA^2), o, weights = ww)
c0 <- coef(m)[[1]]
a_fit <- coef(m)[2:3] / c0
cat(sprintf("refit a1 = %.4g, a2 = %.4g (vza_correction.R: %.4g, %.4g)\n",
            a_fit[1], a_fit[2], vza_coef[["757"]][["a1"]], vza_coef[["757"]][["a2"]]))
o[, rel := ratio / c0]                                       # SIF relative to nadir

bins <- o[, .(rel = sum(ww * rel) / sum(ww), se = sqrt(1 / sum(ww)) / c0, n = .N),
          by = .(vza = pmin(floor(VZA / 5) * 5, 65) + 2.5)][order(vza)]
daybins <- o[, .(rel = sum(ww * rel) / sum(ww), n = .N),
             by = .(date, vza = pmin(floor(VZA / 5) * 5, 65) + 2.5)][n >= 20]
daybins[, vza_j := vza + runif(.N, -1.5, 1.5)]               # spread within the bin
curve <- data.table(vza = seq(0, 70, 0.5))[, f := vza_factor(vza, "757")]
ylim <- c(0.5, 2.2)

p <- ggplot() +
  geom_point(data = daybins, aes(vza_j, rel), color = "grey60", size = 0.5, alpha = 0.6, stroke = 0) +
  geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_errorbar(data = bins, aes(vza, ymin = rel - 2 * se, ymax = rel + 2 * se), width = 0, linewidth = 0.4) +
  geom_point(data = bins, aes(vza, rel), size = 1.4) +
  geom_line(data = curve, aes(vza, f), linetype = "22", linewidth = 0.6, color = "#b2182b") +
  coord_cartesian(xlim = c(0, 71), ylim = ylim, expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  labs(x = "Viewing zenith angle (\u00b0)", y = expression(SIF[757] ~ "relative to nadir view")) +
  theme_bw(base_size = fs_title, base_family = fig_font) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = fs_tick - 0.5),
        plot.margin = margin(6, 8, 4, 4))
save_fig(p, "fig_s1_vza", height = 3, width = 4)
cat(sprintf("%d soundings, %d days; %.1f%% of day-bins outside the y-range shown\n", nrow(o), uniqueN(o$date),
            100 * daybins[, mean(rel < ylim[1] | rel > ylim[2])]))
cat(sprintf("%d day-bins shown\n", nrow(daybins)))
print(bins, digits = 3)

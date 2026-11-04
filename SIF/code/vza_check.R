# Does the viewing zenith angle (VZA) affect retrieved SIF?
#
# Target mode views the same ~25 km area from a wide range of VZAs within one
# overpass (here ~4-70 deg), so each target day is a natural experiment: same
# day, same surface, same sun, different view angle. Comparing across days
# instead would mostly show season.
#
# Target-mode soundings with Quality_Flag 0/1. Tests:
#   a. All soundings: SIF757 vs VZA (2-D histogram) - context only, confounded by season
#   b. Within-day: SIF minus that day's (error-weighted) mean, binned by VZA, and the
#      residual of the crop-class model (downscale_sif_by_crop_cdl30.R design, fitted
#      per day), which removes the part of SIF explained by the crops a footprint covers
#   c. Per-day weighted regression slope of SIF on VZA (per 10 deg), as % of the day mean
#   d. Side of the pass: VAz is the azimuth from the footprint to the sensor (early
#      soundings ~155 deg: the satellite is south of the site on its northbound pass).
#      With the sun at ~205 deg, VAz 90-270 = sensor on the sun's side (backscatter),
#      otherwise forward scatter. An effect symmetric between the two sides follows VZA
#      itself rather than the sun-view geometry (e.g. the hotspot).
# Also: same day and grid cell (SIF and VZA demeaned within day x section-road cell of
# the footprint centre), which removes almost all differences in the surface seen;
# VZA and phase angle fitted jointly; relative SIF (SIF / continuum radiance).
#
# Output (vza_check/): soundings_vza.csv, daily_vza_slopes.csv, plot_sif_vs_vza.png,
# and a printed summary.

library(sf)
library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
oco_file <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
design_file <- file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30_design.rds")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
out_dir <- file.path(proj_dir, "vza_check")
quality_flags <- c(0, 1)
prior_sd <- 1.0          # as in downscale_sif_by_crop_cdl30.R
dir.create(out_dir, showWarnings = FALSE)

# ---- Data -------------------------------------------------------------------
os <- st_read(oco_file, quiet = TRUE, query = paste(
  "SELECT Metadata_SoundingId, time_utc, Quality_Flag, VZA, VAz, SZA, SAz,",
  "Science_SIF_757nm AS sif, Science_SIF_Uncertainty_757nm AS unc,",
  "Science_SIF_Relative_757nm AS sif_rel, geom",
  "FROM oco2_sif WHERE Metadata_MeasurementMode = 2"))
os <- os[os$Quality_Flag %in% quality_flags & !is.na(os$sif) & !is.na(os$unc), ]
grid <- st_read(grid_file, quiet = TRUE)
grid$cell <- seq_len(nrow(grid))
cent <- st_centroid(st_convex_hull(st_transform(st_geometry(os), st_crs(grid))))
os$cell <- st_join(st_sf(geometry = cent), grid[, "cell"])$cell
o <- as.data.table(st_drop_geometry(os))
o[, `:=`(sid = as.character(Metadata_SoundingId), date = as.Date(substr(time_utc, 1, 10)))]

# phase angle: angle between the directions to the sun and to the sensor, seen
# from the footprint (0 = sensor looking straight down-sun = hotspot)
rad <- pi / 180
o[, phase := acos(pmin(1, cos(SZA * rad) * cos(VZA * rad) +
                          sin(SZA * rad) * sin(VZA * rad) * cos((SAz - VAz) * rad))) / rad]
o[, w := 1 / unc^2]
o[, day_mean := sum(w * sif) / sum(w), by = date]
o[, anom := sif - day_mean]
o[, side := fifelse(VAz >= 90 & VAz < 270, "Sensor on the sun's side", "Sensor opposite the sun")]

# crop-class model residuals (per day, all that day's footprints in the fit)
design <- readRDS(design_file)
fit_res <- rbindlist(lapply(names(design), function(d) {
  ds <- design[[d]]
  Kc <- ds$Kc[, colSums(ds$Kc) > 0, drop = FALSE]
  xa <- sum(ds$y / ds$se2) / sum(1 / ds$se2)
  KtSe <- t(Kc / ds$se2)
  S_hat <- solve(KtSe %*% Kc + diag(1 / prior_sd^2, ncol(Kc)))
  x <- as.vector(xa + S_hat %*% (KtSe %*% (ds$y - Kc %*% rep(xa, ncol(Kc)))))
  data.table(sid = as.character(ds$sounding_id), crop_resid = as.vector(ds$y - Kc %*% x))
}))
o <- merge(o, fit_res, by = "sid", all.x = TRUE)
fwrite(o[, .(sid, date, cell, VZA, VAz, SZA, SAz, phase, side, sif, unc, sif_rel, anom, crop_resid)],
       file.path(out_dir, "soundings_vza.csv"))

# ---- Per-day slopes ---------------------------------------------------------
day_slopes <- o[, {
  cf <- summary(lm(sif ~ VZA, weights = w))$coefficients
  cg <- summary(lm(sif ~ VZA + phase, weights = w))$coefficients
  rr <- if (sum(!is.na(crop_resid)) > 50) summary(lm(crop_resid ~ VZA, weights = w))$coefficients else NULL
  .(n = .N, vza_min = min(VZA), vza_max = max(VZA), day_mean = day_mean[1],
    slope10 = 10 * cf["VZA", 1], slope10_se = 10 * cf["VZA", 2],
    slope10_pct = 100 * 10 * cf["VZA", 1] / day_mean[1],
    slope10_given_phase = 10 * cg["VZA", 1], slope10_given_phase_se = 10 * cg["VZA", 2],
    phase_slope10_given_vza = 10 * cg["phase", 1], phase_slope10_given_vza_se = 10 * cg["phase", 2],
    resid_slope10 = if (is.null(rr)) NA_real_ else 10 * rr["VZA", 1],
    resid_slope10_se = if (is.null(rr)) NA_real_ else 10 * rr["VZA", 2],
    r_vza = cor(sif, VZA), r_phase = cor(sif, phase), r_vza_phase = cor(VZA, phase))
}, by = date]
setorder(day_slopes, date)
fwrite(day_slopes, file.path(out_dir, "daily_vza_slopes.csv"))

# ---- Pooled within-day fits (day fixed effects), weighted --------------------
pool <- summary(lm(sif ~ VZA + factor(date), data = o, weights = w))$coefficients["VZA", ]
pool_both <- summary(lm(sif ~ VZA + phase + factor(date), data = o, weights = w))$coefficients
pool_res <- summary(lm(crop_resid ~ VZA, data = o[!is.na(crop_resid)], weights = w))$coefficients["VZA", ]
pool_rel <- summary(lm(sif_rel ~ VZA + factor(date), data = o[!is.na(sif_rel)], weights = w))$coefficients["VZA", ]
pool_side <- o[, {
  cc <- summary(lm(sif ~ VZA + factor(date), weights = w))$coefficients["VZA", ]
  .(slope10 = 10 * cc[1], se10 = 10 * cc[2], n = .N)
}, by = side]
# same day, same grid cell: weighted demeaning within (date, cell)
cf <- o[!is.na(cell)]
cf[, grp := .N, by = .(date, cell)]
cf <- cf[grp >= 5]
cf[, `:=`(sif_d = sif - sum(w * sif) / sum(w), vza_d = VZA - sum(w * VZA) / sum(w)), by = .(date, cell)]
pool_cell <- summary(lm(sif_d ~ 0 + vza_d, data = cf, weights = w))$coefficients["vza_d", ]

# ---- Figure -----------------------------------------------------------------
deg_lab <- function(s) parse(text = s)
bins <- function(dt, var, xvar, width = 5) {
  dt[!is.na(get(var)), .(x = (floor(get(xvar) / width) + 0.5) * width,
                         m = sum(w * get(var)) / sum(w), se = sqrt(1 / sum(w)), n = .N),
     by = .(xb = floor(get(xvar) / width))][n >= 200][order(x)]
}
dsif <- expression(Delta*SIF[757]~(mW~m^-2~sr^-1~nm^-1))
vza_bins_lab <- expression("Viewing zenith angle ("*degree*"), 5"*degree~"bins")

pa <- ggplot(o, aes(VZA, sif)) +
  geom_bin2d(bins = c(60, 80)) +
  scale_fill_viridis_c(trans = "log10", name = "Soundings") +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.3) +
  coord_cartesian(ylim = quantile(o$sif, c(0.005, 0.995))) +
  labs(x = expression("Viewing zenith angle ("*degree*")"),
       y = expression(SIF[757]~(mW~m^-2~sr^-1~nm^-1)), title = "All target-mode soundings") +
  theme_bw(base_size = 9)

ba <- rbind(bins(o, "anom", "VZA")[, what := "SIF minus day mean"],
            bins(o, "crop_resid", "VZA")[, what := "Crop-class model residual"])
pb <- ggplot(ba, aes(x, m, color = what)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = m - 2 * se, ymax = m + 2 * se), size = 0.25,
                  position = position_dodge(width = 1.5)) +
  geom_line(position = position_dodge(width = 1.5), linewidth = 0.4) +
  scale_color_manual(values = c("SIF minus day mean" = "grey30", "Crop-class model residual" = "#1a9850")) +
  labs(x = vza_bins_lab, y = dsif, color = NULL,
       title = "Within each target day (weighted mean, 2 SE bars)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")

lo <- 100 * (day_slopes$slope10 - 2 * day_slopes$slope10_se) / day_slopes$day_mean
hi <- 100 * (day_slopes$slope10 + 2 * day_slopes$slope10_se) / day_slopes$day_mean
pc <- ggplot(day_slopes, aes(x = date, y = slope10_pct)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = 100 * (slope10 - 2 * slope10_se) / day_mean,
                      ymax = 100 * (slope10 + 2 * slope10_se) / day_mean), size = 0.2) +
  coord_cartesian(ylim = c(min(0, quantile(lo, 0.01)), quantile(hi, 0.99))) +
  labs(x = NULL, y = expression("Change in SIF per 10"*degree~"VZA (% of day mean)"),
       title = "Per-day slope of SIF on VZA (2 SE bars)") +
  theme_bw(base_size = 9)

bs <- rbindlist(lapply(split(o, o$side), function(d) bins(d, "anom", "VZA")[, side := d$side[1]]))
pd <- ggplot(bs, aes(x, m, color = side)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = m - 2 * se, ymax = m + 2 * se), size = 0.25,
                  position = position_dodge(width = 1.5)) +
  geom_line(position = position_dodge(width = 1.5), linewidth = 0.4) +
  scale_color_manual(values = c("Sensor on the sun's side" = "#d95f02", "Sensor opposite the sun" = "#1f78b4")) +
  labs(x = vza_bins_lab, y = dsif, color = NULL, title = "Within-day anomaly by side of the pass") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")

fig <- (pa + pb) / (pc + pd) + plot_annotation(tag_levels = "a", tag_suffix = ".")
ggsave(file.path(out_dir, "plot_sif_vs_vza.png"), fig, width = 10, height = 8, dpi = 300)

# ---- Summary ----------------------------------------------------------------
mean_sif <- weighted.mean(o$sif, o$w)
cat(sprintf("\n%s soundings on %d target days; VZA %.0f-%.0f deg (median per-day range %.0f-%.0f)\n",
            format(nrow(o), big.mark = ","), uniqueN(o$date), min(o$VZA), max(o$VZA),
            median(day_slopes$vza_min), median(day_slopes$vza_max)))
cat(sprintf("Pooled within-day slope: %.4f +/- %.4f per 10 deg (mean SIF %.2f -> %.1f%% per 10 deg)\n",
            10 * pool[1], 10 * pool[2], mean_sif, 100 * 10 * pool[1] / mean_sif))
cat(sprintf("  with phase angle in the model: VZA %.4f +/- %.4f, phase %.4f +/- %.4f per 10 deg\n",
            10 * pool_both["VZA", 1], 10 * pool_both["VZA", 2], 10 * pool_both["phase", 1], 10 * pool_both["phase", 2]))
cat(sprintf("  crop-class model residuals: %.4f +/- %.4f per 10 deg\n", 10 * pool_res[1], 10 * pool_res[2]))
cat(sprintf("  same day and grid cell (%s soundings): %.4f +/- %.4f per 10 deg\n",
            format(nrow(cf), big.mark = ","), 10 * pool_cell[1], 10 * pool_cell[2]))
for (i in seq_len(nrow(pool_side))) cat(sprintf("  %s: %.4f +/- %.4f per 10 deg (n = %s)\n", pool_side$side[i],
            pool_side$slope10[i], pool_side$se10[i], format(pool_side$n[i], big.mark = ",")))
rel_mean <- weighted.mean(o$sif_rel, o$w, na.rm = TRUE)
cat(sprintf("  relative SIF (SIF / continuum): %.1f%% of its mean per 10 deg\n", 100 * 10 * pool_rel[1] / rel_mean))
cat(sprintf("Per-day slopes: median %.1f%% per 10 deg (IQR %.1f to %.1f); %d of %d days significantly < 0, %d > 0 (|t| > 2)\n",
            median(day_slopes$slope10_pct), quantile(day_slopes$slope10_pct, 0.25), quantile(day_slopes$slope10_pct, 0.75),
            sum(day_slopes$slope10 / day_slopes$slope10_se < -2), nrow(day_slopes),
            sum(day_slopes$slope10 / day_slopes$slope10_se > 2)))
cat(sprintf("Median within-day correlation of VZA with phase angle: %.2f\n", median(day_slopes$r_vza_phase)))
cat(sprintf("Median per-day |r(SIF, VZA)|: %.3f\n", median(abs(day_slopes$r_vza))))

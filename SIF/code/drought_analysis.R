# This script pulls in drought status for the area around Lamont during the time
# of the study do see if SIF correlates with drought.
#
#
# 1. Average seasonality: for each crop class, a smooth curve of SIF against
#    day of year from all years (2014-2024) of OCO-2 target days.
# 2. Residuals: every crop-class value minus that curve.
# 3. Drought: residuals compared with US Drought Monitor status on the same
#    date, by season.
#
# SIF comes from the latest inversion, downscale_sif_by_crop_cdl30.R (one value
# per crop class per target day, from 30 m CDL fractions). The previous
# version of this analysis, using downscaled_JFD/, is drought_analysis_JFD.R.
#
# Drought: the USDM export (categorical percentages per county) for the four
# counties the study grid spans, averaged with weights = each county's share
# of the grid (Grant 31%, Garfield 27%, Kay 22%, Noble 19%). Drought index
# Drought index: the USDM Drought Severity and Coverage Index on its standard
# 0-500 scale, DSCI = D0 + 2 D1 + 3 D2 + 4 D3 + 5 D4 from the categorical (non-
# cumulative) percentages (= the sum of the cumulative percentages in D0+ ... D4).
# For box plots, target days are split into three DSCI classes with physical
# meaning: none (DSCI < 1, no part of the area even abnormally dry), dry (1-199),
# and drought (>= 200, equivalent to the whole area in at least moderate drought,
# D1). Cut-offs and labels are written to drought_dsci_classes.csv for
# compare_ndvi.R and make_plots.R. (Terciles of the target days were used before;
# their upper cut-off, 172, split a cluster of targets and had no meaning.)
#
# Outputs (in proj_dir):
#   drought_index_study_area.csv    weekly USDM categories and index for the study grid
#   drought_residuals.csv           every class-day: SIF, seasonal fit, residual, drought
#   drought_residual_stats.csv      residual vs drought statistics per crop (and season)
#   plot_drought_usdm.png           drought categories over the study area by year
#   plot_drought_seasonality.png    average seasonality per crop with all years
#   plot_drought_resid_boxplots.png residuals by season and drought status
#   plot_drought_resid_vs_dsci.png  residuals against the drought index, by season
#   plot_drought_dsci_boxes.png     target-mean residual by DSCI class: all targets and by season

library(ggplot2)
library(data.table)
library(sf)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")

# ---- Settings ---------------------------------------------------------------
sif_file    <- file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv")
usdm_file   <- file.path(proj_dir, "dm_export_20140101_20261001.csv")
county_file <- file.path(proj_dir, "lamont_counties.gpkg")   # Census cb_2023 county outlines
grid_file   <- file.path(proj_dir, "sub_road_grid.shp")

analysis_classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
                      "Grassland/Pasture", "Sorghum", "Cotton")
min_avg_kernel <- 0.9   # keep class-days determined by the data, not the prior
loess_span     <- 0.4   # smoothness of the seasonal curve (fraction of a year)
leave_year_out <- TRUE  # fit each year's seasonal curve without that year, so a
                        # drought year doesn't pull the average toward itself

sif_label   <- expression(SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1))
resid_label <- expression(Residual~SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1))
season_of <- function(m) fcase(m %in% c(12, 1, 2), "Winter", m %in% 3:5, "Spring",
                               m %in% 6:8, "Summer", default = "Autumn")

# ---- Drought over the study grid --------------------------------------------
grid <- st_read(grid_file, quiet = TRUE)
counties <- st_transform(st_read(county_file, quiet = TRUE), st_crs(grid))
sf_use_s2(FALSE)
overlap <- suppressWarnings(st_intersection(counties, st_union(st_geometry(grid))))
weights <- data.table(County = paste(overlap$NAME, "County"),
                      w = as.numeric(st_area(overlap)) / sum(as.numeric(st_area(overlap))))

usdm <- fread(usdm_file)[State == "OK"]
usdm <- merge(usdm, weights, by = "County")
cats <- c("None", "D0", "D1", "D2", "D3", "D4")
drought <- usdm[, lapply(.SD, function(x) sum(x * w) / sum(w)), by = MapDate, .SDcols = cats]
drought[, date := as.Date(as.character(MapDate), "%Y%m%d")]
drought[, dsci := D0 + 2 * D1 + 3 * D2 + 4 * D3 + 5 * D4]   # USDM DSCI, 0-500
setorder(drought, date)
fwrite(drought[, .(date, None, D0, D1, D2, D3, D4, dsci)],
       file.path(proj_dir, "drought_index_study_area.csv"))   # also used by compare_ndvi.R

# ---- Crop-class SIF and its average seasonality -----------------------------
sif <- fread(sif_file)[crop_class %in% analysis_classes & avg_kernel >= min_avg_kernel]
sif[, date := as.Date(date)]
sif[, `:=`(year = as.integer(format(date, "%Y")),
           doy = as.integer(format(date, "%j")),
           month = as.integer(format(date, "%m")))]
sif[, season := season_of(month)]
sif[, crop_class := factor(crop_class, levels = analysis_classes)]

# Seasonal curve: loess of SIF on day of year, weighted by 1 / SD^2, with the
# data repeated a year either side so December and January join smoothly
seasonal_fit <- function(fit_data, new_doy) {
  wrap <- rbind(fit_data, copy(fit_data)[, doy := doy - 365],
                copy(fit_data)[, doy := doy + 365])
  m <- loess(sif757 ~ doy, data = wrap, weights = 1 / sif757_sd_scaled^2,
             span = loess_span / 3, degree = 2)
  predict(m, newdata = data.frame(doy = new_doy))
}

sif[, sif_fit := {
  out <- numeric(.N)
  for (y in unique(year)) {
    rows <- year == y
    train <- if (leave_year_out) .SD[!rows] else .SD
    out[rows] <- seasonal_fit(train, doy[rows])
  }
  out
}, by = crop_class, .SDcols = c("doy", "sif757", "sif757_sd_scaled", "year")]
sif[, residual := sif757 - sif_fit]
sif[, residual_z := residual / sif757_sd_scaled]   # in units of the estimate's SD

# Smooth curve for plotting (all years)
curves <- sif[, .(doy = 1:365, sif_fit = seasonal_fit(.SD, 1:365)),
              by = crop_class, .SDcols = c("doy", "sif757", "sif757_sd_scaled")]

# ---- Attach drought status --------------------------------------------------
# USDM maps are weekly; interpolate linearly between maps to each SIF date
sif[, dsci := approx(as.numeric(drought$date), drought$dsci, as.numeric(date), rule = 2)$y]
# Fixed DSCI classes (see header); same columns as the earlier terciles file
terciles <- data.table(low_upper = 1, high_lower = 200,
                       label_low  = "None (DSCI 0)",
                       label_mid  = "Dry (DSCI 1-199)",
                       label_high = "Drought (DSCI >= 200)")
fwrite(terciles, file.path(proj_dir, "drought_dsci_classes.csv"))
dsci_class <- function(x) factor(fcase(x < terciles$low_upper, terciles$label_low,
                                       x < terciles$high_lower, terciles$label_mid,
                                       default = terciles$label_high),
                                 levels = unlist(terciles[, .(label_low, label_mid, label_high)]))
sif[, drought := dsci_class(dsci)]
fwrite(sif, file.path(proj_dir, "drought_residuals.csv"))

# ---- Statistics -------------------------------------------------------------
# Spearman correlation of residual with the drought index, and a
# Kruskal-Wallis test of residuals across drought groups
test_groups <- function(r, g) {
  g <- droplevels(g)
  if (nlevels(g) < 2 || any(table(g) < 3)) return(NA_real_)
  kruskal.test(r, g)$p.value
}
stats_crop <- sif[, .(n = .N,
                      spearman_r = cor(residual, dsci, method = "spearman"),
                      spearman_p = suppressWarnings(cor.test(residual, dsci, method = "spearman")$p.value),
                      kruskal_p = test_groups(residual, drought)),
                  by = crop_class][, season := "All"]
stats_season <- sif[, .(n = .N,
                        spearman_r = if (.N > 4) cor(residual, dsci, method = "spearman") else NA_real_,
                        spearman_p = if (.N > 4) suppressWarnings(cor.test(residual, dsci, method = "spearman")$p.value) else NA_real_,
                        kruskal_p = test_groups(residual, drought)),
                    by = .(crop_class, season)]
stats <- rbind(stats_crop, stats_season)
setorder(stats, crop_class, season)
fwrite(stats, file.path(proj_dir, "drought_residual_stats.csv"))
print(stats_crop, digits = 2)

# Residuals of different crops on the same day are correlated (shared
# day-to-day effects: viewing geometry, clouds, regional weather), so the
# per-crop tests above are not independent and overstate the evidence. Test
# at the level of days instead: (a) the mean residual across crops, and
# (b) each crop's residual relative to the other crops that day, which asks
# whether some crops respond to drought more than others.
spearman <- function(x, y) {
  t <- suppressWarnings(cor.test(x, y, method = "spearman"))
  data.table(n = length(x), spearman_r = unname(t$estimate), spearman_p = t$p.value)
}
day_resid <- sif[, .(residual = mean(residual), dsci = dsci[1], year = year[1]), by = date]
same_day_r <- cor(as.matrix(dcast(sif, date ~ crop_class, value.var = "residual")[, -1]),
                  use = "pairwise.complete.obs")
sif[, residual_vs_crops := residual - mean(residual), by = date]
stats_day <- rbind(
  cbind(test = "Day-mean residual vs drought index", spearman(day_resid$residual, day_resid$dsci)),
  cbind(test = "Day-mean residual vs drought index, without 2022-2023",
        day_resid[!year %in% 2022:2023, spearman(residual, dsci)]),
  sif[, cbind(test = paste("Relative to other crops:", crop_class[1]),
              spearman(residual_vs_crops, dsci)), by = crop_class][, crop_class := NULL])
fwrite(stats_day, file.path(proj_dir, "drought_residual_stats_daylevel.csv"))
cat(sprintf("\nMean correlation between crops' residuals on the same day: %.2f\n",
            mean(same_day_r[upper.tri(same_day_r)], na.rm = TRUE)))
print(stats_day, digits = 2)

# ---- Plots ------------------------------------------------------------------
# Drought categories over the study area, by year
ddf <- melt(drought[date >= as.Date("2014-01-01")], id.vars = "date",
            measure.vars = c("D0", "D1", "D2", "D3", "D4"), variable.name = "category")
ddf[, `:=`(year = format(date, "%Y"), doy = as.integer(format(date, "%j")))]
my_discrete_palette <- c('#fef0d9', '#fdcc8a', '#fc8d59', '#e34a33', '#b30000')
g <- ggplot(ddf, aes(x = doy, y = value, fill = category)) +
  geom_area() +
  scale_fill_manual(values = my_discrete_palette) +
  facet_wrap(~year, nrow = 3) +
  geom_vline(xintercept = c(2, 5, 8, 11) / 12 * 365.25, linetype = 2) +
  labs(x = "Day of Year", y = "Percent of study area", fill = "") +
  theme_bw()
ggsave(file.path(proj_dir, "plot_drought_usdm.png"), plot = g, width = 10, height = 6, dpi = 300)

# Average seasonality per crop with every year's values
g <- ggplot(sif, aes(x = doy, y = sif757, color = factor(year))) +
  geom_point() +
  geom_line(data = curves, aes(x = doy, y = sif_fit), color = "red", linewidth = 1.5,
            inherit.aes = FALSE) +
  labs(x = "Day of Year", y = sif_label, color = "") +
  facet_wrap(~crop_class) +
  guides(color = guide_legend(ncol = 2)) +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.15))
ggsave(file.path(proj_dir, "plot_drought_seasonality.png"), plot = g,
       width = 10, height = 7, dpi = 300)

# Residuals by season and drought status
three_pal <- c('#ffffb2', '#fecc5c', '#e31a1c')
g <- ggplot(sif, aes(x = season, y = residual, fill = drought)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 0.6, alpha = 0.5) +
  facet_wrap(~crop_class, ncol = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = three_pal, drop = FALSE) +
  scale_x_discrete(limits = c("Spring", "Summer", "Autumn", "Winter")) +
  labs(x = "Season", y = resid_label, fill = "DSCI class") +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(.85, .12))
ggsave(file.path(proj_dir, "plot_drought_resid_boxplots.png"), plot = g,
       width = 10, height = 8, dpi = 300)

# Residuals against the drought index, by season with all crops together.
# Crop residuals on the same day are correlated, so the fit and statistics use
# the target-mean residual (one point per target); crop points are background.
crop_class_colors <- c(   # as in make_plots.R
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d", "Soybeans" = "#1b9e77",
  "Corn" = "#d95f02", "Sorghum" = "#e7298a", "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e")
seasons <- c("Spring", "Summer", "Autumn", "Winter")
day_season <- sif[, .(residual = mean(residual), dsci = dsci[1], drought = drought[1], season = season[1], year = year[1]), by = date]
day_season[, `:=`(season = factor(season, levels = seasons),
                  period = ifelse(year %in% 2022:2023, "2022-2023 drought", "Other years"))]
season_stats <- day_season[, {
  t <- suppressWarnings(cor.test(residual, dsci, method = "spearman"))
  t2 <- suppressWarnings(cor.test(residual[period == "Other years"], dsci[period == "Other years"], method = "spearman"))
  .(n = .N, spearman_r = unname(t$estimate), spearman_p = t$p.value,
    spearman_r_no2223 = unname(t2$estimate), spearman_p_no2223 = t2$p.value)
}, keyby = season]
fwrite(season_stats, file.path(proj_dir, "drought_residual_stats_season.csv"))
print(season_stats, digits = 2)
season_stats[, label := sprintf("r[s] == '%.2f' ~~ (p == '%.2f'*','~~n == %d)", spearman_r, spearman_p, n)]
y_zoom <- c(-0.5, 0.55)   # target means all lie well inside; a few crop points are clipped
cat(sprintf("crop points outside the plotted range: %d of %d\n",
            sif[residual < y_zoom[1] | residual > y_zoom[2], .N], nrow(sif)))
crop_pts <- copy(sif)[, `:=`(season = factor(season, levels = seasons),
                             crop_class = factor(crop_class, levels = names(crop_class_colors)))]
g <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_point(data = crop_pts, aes(dsci, residual, color = crop_class), size = 0.9, alpha = 0.35) +
  geom_smooth(data = day_season, aes(dsci, residual), method = "lm", formula = y ~ x,
              color = "black", fill = "grey70", linewidth = 0.7) +
  geom_point(data = day_season, aes(dsci, residual, shape = period), size = 2.2, fill = "black",
             color = "black", stroke = 0.6) +
  geom_text(data = season_stats, aes(label = label), x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
            parse = TRUE, size = 3) +
  facet_wrap(~season, ncol = 2) +
  coord_cartesian(ylim = y_zoom) +
  scale_color_manual(values = crop_class_colors, name = "Crop class (each target)") +
  scale_shape_manual(values = c("2022-2023 drought" = 24, "Other years" = 21), name = "Target mean") +
  labs(x = "Drought Severity and Coverage Index (DSCI; 0 = none, 500 = all D4)", y = resid_label) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), order = 2),
         shape = guide_legend(order = 1))
ggsave(file.path(proj_dir, "plot_drought_resid_vs_dsci.png"), plot = g,
       width = 10, height = 7.5, dpi = 300)

# Target-mean residual by DSCI class, all crops together: (a) all targets,
# (b) by season. One point per target (crop residuals on a day are correlated);
# Kruskal-Wallis across classes in (a), Spearman with DSCI per season in (b).
kw_all <- kruskal.test(residual ~ drought, data = day_season)
rs_all <- suppressWarnings(cor.test(day_season$residual, day_season$dsci, method = "spearman"))
terc_n <- day_season[, .N, by = .(season, drought)]
box_common <- list(
  geom_hline(yintercept = 0, linetype = 2, color = "grey50"),
  scale_fill_manual(values = three_pal, drop = FALSE, name = "DSCI class"),
  scale_shape_manual(values = c("2022-2023 drought" = 24, "Other years" = 21), name = NULL),
  theme_bw(),
  theme(plot.title = element_text(size = 11, face = "bold"), strip.text = element_text(face = "bold")))
# axis labels computed here (not in a function) so the saved panel is self-contained
box_m <- day_season[, .(med = median(residual), n = .N), keyby = drought]
box_labels <- setNames(paste0(sub(" \\(", "\n(", as.character(box_m$drought)), "\nn = ", box_m$n,
                              ", median ", sprintf("%+.2f", box_m$med)), as.character(box_m$drought))
pa <- ggplot(day_season, aes(drought, residual)) +
  box_common +
  geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.6, alpha = 0.9) +
  geom_point(aes(shape = period), position = position_jitter(width = 0.12, seed = 1),
             fill = "grey15", color = "grey15", size = 1.6, alpha = 0.8) +
  annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.5, size = 3.2,
           label = sprintf("Kruskal-Wallis p = %.3f\nSpearman with DSCI r = %.2f (p = %.3f)\nn = %d targets",
                           kw_all$p.value, rs_all$estimate, rs_all$p.value, nrow(day_season))) +
  scale_x_discrete(labels = box_labels) +
  labs(x = NULL, y = resid_label, title = "All targets (mean over crops)") +
  guides(fill = "none", shape = guide_legend(override.aes = list(size = 2.2)))
pb <- ggplot(day_season, aes(drought, residual)) +
  box_common +
  geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.65, alpha = 0.9) +
  geom_point(aes(shape = period), position = position_jitter(width = 0.12, seed = 1),
             fill = "grey15", color = "grey15", size = 1.3, alpha = 0.8) +
  geom_text(data = terc_n, aes(drought, -Inf, label = paste0("n=", N)), vjust = -0.5, size = 2.6, color = "grey30") +
  geom_text(data = season_stats, aes(x = 0.5, y = Inf, label = label), hjust = 0, vjust = 1.4,
            parse = TRUE, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~season, nrow = 1) +
  scale_x_discrete(labels = c("None", "Dry", "Drought")) +
  labs(x = "DSCI class", y = NULL, title = "By season") +
  guides(fill = "none", shape = "none")
yr <- range(day_season$residual) + c(-0.06, 0.22)
g <- (pa + coord_cartesian(ylim = yr)) + (pb + coord_cartesian(ylim = yr)) +
  plot_layout(widths = c(1, 2.2), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 12, face = "bold"))
saveRDS(list(a = pa, b = pb), file.path(proj_dir, "fig5_panels_drought.rds"))   # for make_fig5.R
ggsave(file.path(proj_dir, "plot_drought_dsci_boxes.png"), plot = g, width = 11, height = 5, dpi = 300)
cat(sprintf("DSCI class boxes: KW p = %.3f; medians %s\n", kw_all$p.value,
            paste(day_season[, .(m = sprintf("%+.3f", median(residual))), keyby = drought]$m, collapse = " / ")))

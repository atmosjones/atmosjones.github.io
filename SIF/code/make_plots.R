# Plots for the cropland SIF downscaling paper.
#
# Reads the downscaled SIF (downscale_sif_inversion.R) with crop types added
# (add_crop_to_sif.R) and writes figures to the project folder.

library(ggplot2)
library(sf)
library(data.table)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
sif_file <- file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg")

sif_label <- expression(SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1))  # = W m-2 sr-1 um-1

# Crop class colors, shared by all plots
crop_class_colors <- c(
  "Winter Wheat" = "#e6ab02", "Dbl Crop WinWht/Soybeans" = "#a6761d",
  "Soybeans" = "#1b9e77", "Corn" = "#d95f02", "Sorghum" = "#e7298a",
  "Cotton" = "#7570b3", "Grassland/Pasture" = "#66a61e",
  "Fallow/Idle Cropland" = "#bdb76b", "Trees/Wetlands" = "#1f4e3d",
  "Developed" = "grey55", "Other" = "grey80")

# ---- Data -------------------------------------------------------------------
sif <- as.data.table(st_drop_geometry(st_read(sif_file, quiet = TRUE)))
sif[, date := as.Date(date)]
sif[, year := as.integer(format(date, "%Y"))]

# ---- Time series of downscaled SIF by crop, one year -----------------------
# Each point is one grid cell on one day; the line and band are the median
# and interquartile range across cells of that crop.
#   min_cells       crops with fewer cells that year are left out
#   min_avg_kernel  drop cells the inversion left close to the prior
#                   (0 keeps all; 0.5 keeps cells mostly set by the data)
plot_sif_timeseries_by_crop <- function(yr, min_cells = 15, min_avg_kernel = 0) {
  d <- sif[year == yr & avg_kernel >= min_avg_kernel]
  crop_cells <- d[, .(n_cells = uniqueN(cell_id)), by = crop_name]
  keep <- crop_cells[n_cells >= min_cells][order(-n_cells), crop_name]
  d <- d[crop_name %in% keep]
  d[, crop_name := factor(crop_name, levels = keep)]

  stats <- d[, .(median = median(sif757),
                 q25 = quantile(sif757, 0.25),
                 q75 = quantile(sif757, 0.75)),
             by = .(crop_name, date)]
  facet_labels <- setNames(sprintf("%s (n = %d)", keep, crop_cells[match(keep, crop_name), n_cells]),
                           keep)

  ggplot(d, aes(x = date, y = sif757)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_jitter(width = 2, height = 0, size = 0.4, alpha = 0.25, color = "grey40") +
    geom_ribbon(data = stats, aes(x = date, ymin = q25, ymax = q75),
                inherit.aes = FALSE, fill = "#1a9850", alpha = 0.3) +
    geom_line(data = stats, aes(x = date, y = median), color = "#1a9850",
              linewidth = 0.8) +
    geom_point(data = stats, aes(x = date, y = median), color = "#1a9850",
               size = 1.2) +
    facet_wrap(~crop_name, labeller = as_labeller(facet_labels)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b",
                 limits = as.Date(paste0(yr, c("-01-01", "-12-31")))) +
    labs(x = NULL, y = sif_label,
         title = sprintf("Downscaled OCO-2 SIF over Lamont, OK, %d, by majority crop", yr),
         subtitle = sprintf("%d OCO-2 target days; points = grid cells, line = median, band = IQR%s",
                            uniqueN(d$date),
                            if (min_avg_kernel > 0) sprintf("; cells with averaging kernel >= %.1f", min_avg_kernel) else "")) +
    theme_bw() +
    theme(strip.text = element_text(size = 8),
          plot.subtitle = element_text(size = 9))
}

g <- plot_sif_timeseries_by_crop(2019)
ggsave(file.path(proj_dir, "plot_sif_timeseries_2019_by_crop.png"), plot = g,
       width = 10, height = 7, dpi = 300)

# Same, keeping only cells mostly set by the data (fewer cells per crop, so a
# lower minimum); which cells qualify changes from day to day with the scan
g <- plot_sif_timeseries_by_crop(2019, min_cells = 5, min_avg_kernel = 0.5)
ggsave(file.path(proj_dir, "plot_sif_timeseries_2019_by_crop_ak05.png"), plot = g,
       width = 10, height = 7, dpi = 300)

# ---- Time series of SIF estimated per crop class, one year ------------------
# From downscale_sif_by_crop.R: one SIF value per crop class per day, with
# +/- 2 SD error bars (posterior SD scaled by the fit's chi2).
crop_est <- fread(file.path(proj_dir, "downscaled_sif757", "sif757_by_crop.csv"))
crop_est[, date := as.Date(date)]

plot_crop_sif_timeseries <- function(yr, min_avg_kernel = 0.9,
                                     drop_classes = c("Developed", "Other")) {
  d <- crop_est[format(date, "%Y") == yr & avg_kernel >= min_avg_kernel &
                !crop_class %in% drop_classes]
  order_by_area <- d[, .(area = mean(area_km2)), by = crop_class][order(-area), crop_class]
  d[, crop_class := factor(crop_class, levels = order_by_area)]

  ggplot(d, aes(x = date, y = sif757)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(color = "#1a9850", linewidth = 0.8) +
    geom_errorbar(aes(ymin = sif757 - 2 * sif757_sd_scaled,
                      ymax = sif757 + 2 * sif757_sd_scaled),
                  width = 6, color = "#1a9850") +
    geom_point(color = "#1a9850", size = 1.5) +
    facet_wrap(~crop_class) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b",
                 limits = as.Date(paste0(yr, c("-01-01", "-12-31")))) +
    labs(x = NULL, y = sif_label,
         title = sprintf("OCO-2 SIF over Lamont, OK, %d, estimated per crop class", yr),
         subtitle = sprintf("%d target days; one inversion per day with a SIF value per crop class; error bars = 2 SD",
                            uniqueN(d$date))) +
    theme_bw() +
    theme(strip.text = element_text(size = 8),
          plot.subtitle = element_text(size = 9))
}

g <- plot_crop_sif_timeseries(2019)
ggsave(file.path(proj_dir, "plot_sif_by_crop_class_2019.png"), plot = g,
       width = 10, height = 7, dpi = 300)

# ---- Crop-class SIF: 30 m CDL fractions vs grid-cell majority crop ---------
# downscale_sif_by_crop_cdl30.R (exact 30 m crop fractions per footprint)
# against downscale_sif_by_crop.R (each grid cell = its majority crop)
crop_est_30 <- fread(file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv"))
crop_est_30[, date := as.Date(date)]

plot_crop_sif_methods <- function(yr, min_avg_kernel = 0.9,
                                  drop_classes = c("Other")) {
  d <- rbind(crop_est[, .(date, crop_class, sif757, sif757_sd_scaled, avg_kernel,
                          method = "Grid cells (majority crop)")],
             crop_est_30[, .(date, crop_class, sif757, sif757_sd_scaled, avg_kernel,
                             method = "30 m CDL fractions")])
  d <- d[format(date, "%Y") == yr & avg_kernel >= min_avg_kernel &
         !crop_class %in% drop_classes]
  order_by_area <- crop_est_30[format(date, "%Y") == yr, .(a = mean(area_km2)),
                              by = crop_class][order(-a), crop_class]
  d[, crop_class := factor(crop_class, levels = intersect(order_by_area, unique(crop_class)))]
  d[, method := factor(method, levels = c("30 m CDL fractions", "Grid cells (majority crop)"))]

  ggplot(d, aes(x = date, y = sif757, color = method)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = sif757 - 2 * sif757_sd_scaled,
                      ymax = sif757 + 2 * sif757_sd_scaled), width = 6) +
    geom_point(size = 1.3) +
    facet_wrap(~crop_class) +
    scale_color_manual(values = c("30 m CDL fractions" = "#1a9850",
                                  "Grid cells (majority crop)" = "grey55")) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b",
                 limits = as.Date(paste0(yr, c("-01-01", "-12-31")))) +
    labs(x = NULL, y = sif_label, color = NULL,
         title = sprintf("OCO-2 SIF over Lamont, OK, %d, estimated per crop class", yr),
         subtitle = "Error bars = 2 SD") +
    theme_bw() +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
}

g <- plot_crop_sif_methods(2019)
ggsave(file.path(proj_dir, "plot_sif_by_crop_class_2019_cdl30.png"), plot = g,
       width = 10, height = 8, dpi = 300)

# ---- The two downscalings: crop-class inversion vs ~800 m grid --------------
# a: one year - crop-class estimates (downscale_sif_by_crop_cdl30.R, +/- 2 SD)
#    against the grid cells (downscale_sif_inversion.R) whose majority crop is
#    that class: median and IQR across cells with averaging kernel >= min_ak.
# b: all days - within-day crop contrasts (crop minus that day's mean over the
#    crops present in both) of the grid medians against the crop-class values.
#    Grid cells are mixed (the majority crop covers ~30% of a typical cell) and
#    pulled toward the day mean by the prior, so their contrasts should be weaker.
cmp_classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
                 "Sorghum", "Cotton", "Grassland/Pasture")
grid_by_crop <- function(min_ak, min_cells = 3) {
  sif[crop_name %in% cmp_classes & avg_kernel >= min_ak,
      .(grid_med = median(sif757), grid_q25 = quantile(sif757, 0.25),
        grid_q75 = quantile(sif757, 0.75), n_cells = .N),
      by = .(date, crop_class = crop_name)][n_cells >= min_cells]
}
paired_contrasts <- function(min_ak) {
  g <- grid_by_crop(min_ak)
  cc <- crop_est_30[crop_class %in% cmp_classes & avg_kernel >= 0.9,
                    .(date, crop_class, cls = sif757)]
  m <- merge(g, cc, by = c("date", "crop_class"))
  m[, n_crops := .N, by = date]
  m <- m[n_crops >= 3]
  m[, `:=`(grid_anom = grid_med - mean(grid_med), cls_anom = cls - mean(cls)), by = date]
  m
}

plot_grid_vs_crop_class <- function(yr = 2019, min_ak = 0.5) {
  g <- grid_by_crop(min_ak)[format(date, "%Y") == yr]
  cc <- crop_est_30[format(date, "%Y") == yr & crop_class %in% cmp_classes & avg_kernel >= 0.9]
  lv <- cmp_classes[cmp_classes %in% c(g$crop_class, cc$crop_class)]
  g[, crop_class := factor(crop_class, levels = lv)]
  cc[, crop_class := factor(crop_class, levels = lv)]
  cols <- c("Crop-class inversion (30 m CDL fractions)" = "#1a9850",
            "Grid cells of that majority crop (median, IQR)" = "grey35")
  pa <- ggplot() +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_ribbon(data = g, aes(x = date, ymin = grid_q25, ymax = grid_q75), fill = "grey60", alpha = 0.35) +
    geom_line(data = g, aes(x = date, y = grid_med, color = names(cols)[2]), linewidth = 0.6) +
    geom_point(data = g, aes(x = date, y = grid_med, color = names(cols)[2]), size = 1) +
    geom_errorbar(data = cc, aes(x = date, ymin = sif757 - 2 * sif757_sd_scaled,
                                 ymax = sif757 + 2 * sif757_sd_scaled, color = names(cols)[1]), width = 6) +
    geom_line(data = cc, aes(x = date, y = sif757, color = names(cols)[1]), linewidth = 0.6) +
    geom_point(data = cc, aes(x = date, y = sif757, color = names(cols)[1]), size = 1.3) +
    facet_wrap(~crop_class, ncol = 4) +
    scale_color_manual(values = cols) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b",
                 limits = as.Date(paste0(yr, c("-01-01", "-12-31")))) +
    labs(x = NULL, y = sif_label, color = NULL,
         title = sprintf("Crop-class inversion vs ~800 m grid inversion, %d", yr),
         subtitle = sprintf("Grid: cells with averaging kernel >= %.1f and at least 3 cells of the crop that day", min_ak)) +
    theme_bw() + theme(strip.text = element_text(size = 8), legend.position = "bottom",
                       plot.subtitle = element_text(size = 9))

  m <- paired_contrasts(min_ak)
  m_all <- paired_contrasts(0)
  r <- cor(m$cls_anom, m$grid_anom); r_all <- cor(m_all$cls_anom, m_all$grid_anom)
  slope <- coef(lm(grid_anom ~ cls_anom, m))[2]
  m[, crop_class := factor(crop_class, levels = cmp_classes)]
  pb <- ggplot(m, aes(x = cls_anom, y = grid_anom, color = crop_class)) +
    geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 1.3, alpha = 0.75) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.5) +
    scale_color_manual(values = crop_class_colors) +
    coord_equal() +
    labs(x = "Crop-class inversion: crop minus day mean",
         y = "Grid inversion: crop median minus day mean", color = NULL,
         title = "Within-day crop contrasts, all target days",
         subtitle = sprintf("r = %.2f, slope = %.2f (%d crop-days; cells with averaging kernel >= %.1f)\nall cells: r = %.2f",
                            r, slope, nrow(m), min_ak, r_all)) +
    theme_bw() + theme(plot.subtitle = element_text(size = 9))
  cat(sprintf("grid vs crop-class contrasts: r = %.2f, slope = %.2f, n = %d (AK >= %.1f); all cells r = %.2f, slope = %.2f, n = %d\n",
              r, slope, nrow(m), min_ak, r_all, coef(lm(grid_anom ~ cls_anom, m_all))[2], nrow(m_all)))
  pa / pb + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "a", tag_suffix = ".")
}
ggsave(file.path(proj_dir, "plot_grid_vs_crop_class.png"), plot_grid_vs_crop_class(),
       width = 10, height = 12, dpi = 300)

# ---- SIFNet comparison (compare_sifnet.R) -----------------------------------
cmp_dir <- file.path(proj_dir, "sifnet_comparison")
sif740_label <- expression(SIF[740] ~~ (mW~m^-2~sr^-1~nm^-1))
main_classes <- c("Winter Wheat", "Grassland/Pasture", "Dbl Crop WinWht/Soybeans",
                  "Soybeans", "Corn", "Sorghum", "Cotton", "Trees/Wetlands")

# Crop-class SIF through time: SIFNet 16-day composites (line) and OCO-2
# crop-class estimates on target days (points, +/- 2 SD)
plot_sifnet_class_timeseries <- function(min_avg_kernel = 0.9) {
  sn <- fread(file.path(cmp_dir, "sifnet_class_means.csv"))[crop_class %in% main_classes]
  oc <- fread(file.path(cmp_dir, "class_comparison.csv"))[crop_class %in% main_classes &
                                                          avg_kernel >= min_avg_kernel]
  sn[, crop_class := factor(crop_class, levels = main_classes)]
  oc[, crop_class := factor(crop_class, levels = main_classes)]
  ggplot() +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(data = sn, aes(x = as.Date(mid_date), y = sifnet_sif740, color = "SIFNet (16-day)"),
              linewidth = 0.6) +
    geom_errorbar(data = oc, aes(x = as.Date(date), ymin = sif740 - 2 * sif740_sd_scaled,
                                 ymax = sif740 + 2 * sif740_sd_scaled,
                                 color = "OCO-2 crop-class inversion"), width = 15) +
    geom_point(data = oc, aes(x = as.Date(date), y = sif740,
                              color = "OCO-2 crop-class inversion"), size = 1.3) +
    facet_wrap(~crop_class, ncol = 2) +
    scale_color_manual(values = c("SIFNet (16-day)" = "grey40",
                                  "OCO-2 crop-class inversion" = "#1a9850")) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
    labs(x = NULL, y = sif740_label, color = NULL,
         title = "Crop-class SIF over Lamont, OK: OCO-2 inversion vs SIFNet") +
    theme_bw() +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
}

# Footprint level: how well each product explains the OCO-2 footprints on each
# day (correlation), against the most any predictor could reach given the
# OCO-2 retrieval noise
plot_sifnet_footprint_skill <- function() {
  fs <- fread(file.path(cmp_dir, "footprint_stats.csv"))
  d <- melt(fs, id.vars = "date", measure.vars = c("r_sifnet", "r_crop_cv", "noise_ceiling"))
  d[, variable := factor(variable, levels = c("r_crop_cv", "r_sifnet", "noise_ceiling"),
                         labels = c("OCO-2 crop-class model (cross-validated)",
                                    "SIFNet", "Noise ceiling"))]
  ggplot(d, aes(x = as.Date(date), y = value, color = variable, shape = variable)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#1a9850", "grey40", "black")) +
    scale_shape_manual(values = c(16, 17, 4)) +
    labs(x = NULL, y = "Correlation with OCO-2 footprint SIF (740 nm)",
         color = NULL, shape = NULL,
         title = "Footprint-level agreement with OCO-2 target observations") +
    theme_bw() + theme(legend.position = "bottom")
}

# Crop contrasts within a day: each class's departure from that day's mean,
# OCO-2 against SIFNet
plot_sifnet_class_contrasts <- function(min_avg_kernel = 0.9) {
  cc <- fread(file.path(cmp_dir, "class_comparison.csv"))[crop_class %in% main_classes &
                                                          avg_kernel >= min_avg_kernel]
  cc[, `:=`(oco_anom = sif740 - mean(sif740), sifnet_anom = sifnet_sif740 - mean(sifnet_sif740)),
     by = date]
  r <- cor(cc$oco_anom, cc$sifnet_anom)
  ggplot(cc, aes(x = sifnet_anom, y = oco_anom, color = crop_class)) +
    geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "SIFNet: class minus day mean", y = "OCO-2: class minus day mean",
         color = NULL, title = "Crop contrasts within each day (SIF 740 nm)",
         subtitle = sprintf("r = %.2f; %d class-days", r, nrow(cc))) +
    theme_bw()
}

ggsave(file.path(proj_dir, "plot_sifnet_class_timeseries.png"),
       plot_sifnet_class_timeseries(), width = 10, height = 9, dpi = 300)
ggsave(file.path(proj_dir, "plot_sifnet_footprint_skill.png"),
       plot_sifnet_footprint_skill(), width = 9, height = 4.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_sifnet_class_contrasts.png"),
       plot_sifnet_class_contrasts(), width = 7, height = 5.5, dpi = 300)

# ---- GOSIF comparison (compare_gosif.R) -------------------------------------
gosif_dir <- file.path(proj_dir, "gosif_comparison")
sif757_daily_label <- expression(Daily~mean~SIF[757] ~~ (mW~m^-2~sr^-1~nm^-1))

# Crop-class SIF through time in daily-mean units (GOSIF's own): GOSIF 8-day
# composites (line) and OCO-2 crop-class estimates x that day's OCO-2 daily
# correction factor (points, +/- 2 SD)
plot_gosif_class_timeseries <- function(min_avg_kernel = 0.9) {
  gs <- fread(file.path(gosif_dir, "gosif_class_means.csv"))[crop_class %in% main_classes]
  oc <- fread(file.path(gosif_dir, "class_comparison.csv"))[crop_class %in% main_classes &
                                                            avg_kernel >= min_avg_kernel]
  oc[, `:=`(oco_daily = sif757 * daily_factor, oco_daily_sd = sif757_sd_scaled * daily_factor)]
  gs[, crop_class := factor(crop_class, levels = main_classes)]
  oc[, crop_class := factor(crop_class, levels = main_classes)]
  # break the line where a class is missing (e.g. no cotton in that year's CDL)
  setorder(gs, crop_class, mid_date)
  gs[, run := cumsum(c(1, diff(as.Date(mid_date)) > 16)), by = crop_class]
  ggplot() +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(data = gs, aes(x = as.Date(mid_date), y = gosif_daily_sif757, color = "GOSIF (8-day)",
                             group = interaction(crop_class, run)),
              linewidth = 0.5) +
    geom_errorbar(data = oc, aes(x = as.Date(date), ymin = oco_daily - 2 * oco_daily_sd,
                                 ymax = oco_daily + 2 * oco_daily_sd,
                                 color = "OCO-2 crop-class inversion"), width = 30) +
    geom_point(data = oc, aes(x = as.Date(date), y = oco_daily,
                              color = "OCO-2 crop-class inversion"), size = 1.3) +
    facet_wrap(~crop_class, ncol = 2) +
    scale_color_manual(values = c("GOSIF (8-day)" = "grey40",
                                  "OCO-2 crop-class inversion" = "#1a9850")) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = sif757_daily_label, color = NULL,
         title = "Crop-class SIF over Lamont, OK: OCO-2 inversion vs GOSIF") +
    theme_bw() +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
}

plot_gosif_footprint_skill <- function() {
  fs <- fread(file.path(gosif_dir, "footprint_stats.csv"))
  d <- melt(fs, id.vars = "date", measure.vars = c("r_gosif", "r_crop_cv", "noise_ceiling"))
  d[, variable := factor(variable, levels = c("r_crop_cv", "r_gosif", "noise_ceiling"),
                         labels = c("OCO-2 crop-class model (cross-validated)",
                                    "GOSIF", "Noise ceiling"))]
  ggplot(d, aes(x = as.Date(date), y = value, color = variable, shape = variable)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#1a9850", "grey40", "black")) +
    scale_shape_manual(values = c(16, 17, 4)) +
    labs(x = NULL, y = "Correlation with OCO-2 footprint SIF (757 nm)",
         color = NULL, shape = NULL,
         title = "Footprint-level agreement with OCO-2 target observations: GOSIF") +
    theme_bw() + theme(legend.position = "bottom")
}

plot_gosif_class_contrasts <- function(min_avg_kernel = 0.9) {
  cc <- fread(file.path(gosif_dir, "class_comparison.csv"))[crop_class %in% main_classes &
                                                            avg_kernel >= min_avg_kernel]
  cc[, `:=`(oco_anom = sif757 - mean(sif757), gosif_anom = gosif_sif757 - mean(gosif_sif757)),
     by = date]
  r <- cor(cc$oco_anom, cc$gosif_anom)
  ggplot(cc, aes(x = gosif_anom, y = oco_anom, color = crop_class)) +
    geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "GOSIF: class minus day mean", y = "OCO-2: class minus day mean",
         color = NULL, title = "Crop contrasts within each day (SIF 757 nm, instantaneous)",
         subtitle = sprintf("r = %.2f; %d class-days", r, nrow(cc))) +
    theme_bw()
}

# Both products together: per-day footprint correlation with OCO-2, against
# the crop-class model at the matching wavelength and the noise ceiling
plot_product_footprint_skill <- function() {
  g <- fread(file.path(gosif_dir, "footprint_stats.csv"))
  s <- fread(file.path(cmp_dir, "footprint_stats.csv"))
  d <- rbind(
    g[, .(band = "757 nm (vs GOSIF)", what = "OCO-2 crop-class model (CV)", r = r_crop_cv)],
    g[, .(band = "757 nm (vs GOSIF)", what = "GOSIF", r = r_gosif)],
    g[, .(band = "757 nm (vs GOSIF)", what = "Noise ceiling", r = noise_ceiling)],
    s[, .(band = "740 nm (vs SIFNet)", what = "OCO-2 crop-class model (CV)", r = r_crop_cv)],
    s[, .(band = "740 nm (vs SIFNet)", what = "SIFNet", r = r_sifnet)],
    s[, .(band = "740 nm (vs SIFNet)", what = "Noise ceiling", r = noise_ceiling)])
  d[, what := factor(what, levels = c("OCO-2 crop-class model (CV)", "SIFNet", "GOSIF", "Noise ceiling"))]
  ggplot(d[!is.na(r)], aes(x = what, y = r, fill = what)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_boxplot(outlier.size = 0.8, width = 0.6) +
    facet_wrap(~band, scales = "free_x") +
    scale_fill_manual(values = c("OCO-2 crop-class model (CV)" = "#1a9850", "SIFNet" = "grey55",
                                 "GOSIF" = "grey80", "Noise ceiling" = "white"), guide = "none") +
    scale_x_discrete(labels = function(x) gsub("OCO-2 crop-class model \\(CV\\)",
                                               "OCO-2 crop-class\nmodel (CV)", x)) +
    labs(x = NULL, y = "Per-day correlation with OCO-2 footprint SIF",
         title = "Footprint-level agreement with OCO-2: crop-class model, SIFNet and GOSIF") +
    theme_bw()
}

ggsave(file.path(proj_dir, "plot_gosif_class_timeseries.png"),
       plot_gosif_class_timeseries(), width = 10, height = 9, dpi = 300)
ggsave(file.path(proj_dir, "plot_gosif_footprint_skill.png"),
       plot_gosif_footprint_skill(), width = 9, height = 4.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_gosif_class_contrasts.png"),
       plot_gosif_class_contrasts(), width = 7, height = 5.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_product_footprint_skill.png"),
       plot_product_footprint_skill(), width = 9, height = 5, dpi = 300)

# ---- Crop area and rotation (crop_area_rotation.R) --------------------------

# Share of the study area in each crop class, by year
plot_crop_area <- function() {
  a <- fread(file.path(proj_dir, "crop_area_by_year.csv"))
  a[, crop_class := factor(crop_class, levels = names(crop_class_colors))]
  ggplot(a, aes(x = factor(year), y = share * 100, fill = crop_class)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = crop_class_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Share of study area (%)", fill = NULL,
         title = "Crop classes in the study area (USDA CDL, 30 m)") +
    theme_bw()
}

# Fraction of cropland pixels whose crop changed from one year to the next
plot_crop_rotation <- function() {
  r <- fread(file.path(proj_dir, "crop_rotation_by_year.csv"))
  d <- melt(r, id.vars = "to_year", measure.vars = c("rotation_rate", "rotation_rate_interior"))
  d[, variable := factor(variable, levels = c("rotation_rate", "rotation_rate_interior"),
                         labels = c("All cropland pixels", "Field interiors (3x3 uniform)"))]
  ggplot(d, aes(x = to_year, y = value * 100, color = variable)) +
    geom_line() + geom_point(size = 2) +
    scale_color_manual(values = c("black", "grey55")) +
    scale_x_continuous(breaks = r$to_year, labels = paste0(r$from_year, "\n", r$to_year)) +
    scale_y_continuous(limits = c(0, 70)) +
    labs(x = NULL, y = "Cropland pixels changing crop (%)", color = NULL,
         title = "Crop rotation rate from one year to the next") +
    theme_bw() + theme(legend.position = "bottom")
}

ggsave(file.path(proj_dir, "plot_crop_area_by_year.png"), plot_crop_area(),
       width = 9, height = 5.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_crop_rotation.png"), plot_crop_rotation(),
       width = 8, height = 4.5, dpi = 300)

# ---- NDVI comparison (compare_ndvi.R) ---------------------------------------
ndvi_dir_cmp <- file.path(proj_dir, "ndvi_comparison")
terciles <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))     # drought_analysis.R
drought_pal <- setNames(c('#ffffb2', '#fecc5c', '#e31a1c'),
                        unlist(terciles[, .(label_low, label_mid, label_high)]))
dsci_class <- function(x) factor(fcase(x < terciles$low_upper, terciles$label_low,
                                       x < terciles$high_lower, terciles$label_mid,
                                       default = terciles$label_high), levels = names(drought_pal))

# NDVI residuals (mean across crops per HLS scene) by season and drought status
plot_ndvi_drought <- function() {
  nd <- fread(file.path(ndvi_dir_cmp, "ndvi_residuals.csv"))
  d <- nd[, .(residual = mean(residual), drought = drought[1], season = season[1]), by = date]
  d[, drought := factor(drought, levels = names(drought_pal))]
  ggplot(d, aes(x = season, y = residual, fill = drought)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 0.7, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = drought_pal, drop = FALSE) +
    scale_x_discrete(limits = c("Spring", "Summer", "Autumn", "Winter")) +
    labs(x = "Season", y = "NDVI residual (mean across crops)", fill = "DSCI class",
         title = sprintf("HLS NDVI departure from the average seasonal cycle (%d scenes)", nrow(d))) +
    theme_bw()
}

# SIF residual against NDVI residual on OCO-2 days (day means across crops)
plot_sif_vs_ndvi_residuals <- function() {
  d <- fread(file.path(ndvi_dir_cmp, "sif_vs_ndvi_residuals_by_day.csv"))
  d[, drought := dsci_class(dsci)]
  r <- cor(d$sif_residual, d$ndvi_residual, method = "spearman")
  ggplot(d, aes(x = ndvi_residual, y = sif_residual)) +
    geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.6) +
    geom_point(aes(fill = drought), shape = 21, size = 2.5) +
    scale_fill_manual(values = drought_pal, drop = FALSE) +
    labs(x = "NDVI residual", y = resid_label_sif, fill = "DSCI class",
         title = "SIF and NDVI departures from their seasonal cycles on OCO-2 days",
         subtitle = sprintf("Spearman r = %.2f; %d days", r, nrow(d))) +
    theme_bw()
}
resid_label_sif <- expression(SIF[757]~residual ~~ (mW~m^-2~sr^-1~nm^-1))

# Footprint level: crop classes vs NDVI vs both, cross-validated
plot_ndvi_footprint_skill <- function() {
  fs <- fread(file.path(ndvi_dir_cmp, "footprint_stats.csv"))
  d <- melt(fs, id.vars = "date", measure.vars = c("r_crop_cv", "r_ndvi_cv", "r_crop_ndvi_cv", "noise_ceiling"))
  d[, variable := factor(variable, levels = c("r_crop_cv", "r_ndvi_cv", "r_crop_ndvi_cv", "noise_ceiling"),
                         labels = c("Crop classes", "NDVI", "Crop classes + NDVI", "Noise ceiling"))]
  ggplot(d, aes(x = variable, y = value)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_boxplot(outlier.shape = NA, fill = "grey95") +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, color = "#1a9850") +
    labs(x = NULL, y = "Correlation with OCO-2 footprint SIF",
         title = sprintf("Footprint-level skill, cross-validated (%d OCO-2 days)", nrow(fs))) +
    theme_bw()
}

# Crop-class SIF against crop-class NDVI on OCO-2 days
plot_sif_vs_ndvi_classes <- function() {
  d <- fread(file.path(ndvi_dir_cmp, "class_comparison.csv"))[!is.na(ndvi)]
  ggplot(d, aes(x = ndvi, y = sif757, color = crop_class)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = crop_class_colors) +
    labs(x = "NDVI (HLS, interpolated to the OCO-2 day)", y = sif_label, color = NULL,
         title = "Crop-class SIF against crop-class NDVI",
         subtitle = sprintf("r = %.2f; %d class-days", cor(d$ndvi, d$sif757), nrow(d))) +
    theme_bw()
}

ggsave(file.path(proj_dir, "plot_ndvi_drought_boxplots.png"), plot_ndvi_drought(),
       width = 8, height = 4.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_sif_vs_ndvi_residuals.png"), plot_sif_vs_ndvi_residuals(),
       width = 7, height = 5, dpi = 300)
ggsave(file.path(proj_dir, "plot_ndvi_footprint_skill.png"), plot_ndvi_footprint_skill(),
       width = 7, height = 4.5, dpi = 300)
ggsave(file.path(proj_dir, "plot_sif_vs_ndvi_classes.png"), plot_sif_vs_ndvi_classes(),
       width = 8, height = 5.5, dpi = 300)

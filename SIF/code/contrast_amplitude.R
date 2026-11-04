# How much of the crop contrast does each product recover?
#
# One measure for every product: the within-day crop contrast (class minus that
# day's mean over the classes compared) of the product, regressed on the same
# contrast from the crop-class inversion (downscale_sif_by_crop_cdl30.R).
#   slope     share of the reference contrast the product recovers (1 = all)
#   sd_ratio  SD of the product's contrasts / SD of the reference contrasts
#   r         correlation of the two
# Products: ~800 m grid inversion (median of majority-crop cells, as in
# make_fig3.R) and SIFnet (740 nm; reference = crop-class inversion at 740 nm).
# NDVI is not in SIF units and is left out. (GOSIF was dropped from the paper.)
# Two sets: (1) each product on all its own class-days; (2) common class-days,
# where both products and the reference are present (SIFnet years, 2018-2021).
# 95% intervals from 2,000 bootstrap resamples of targets.
# Output: contrast_amplitude/contrast_amplitude.csv, plot_contrast_amplitude.png

library(sf)
library(data.table)
library(ggplot2)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "contrast_amplitude")
dir.create(out_dir, showWarnings = FALSE)
classes <- c("Winter Wheat", "Dbl Crop WinWht/Soybeans", "Soybeans", "Corn",
             "Sorghum", "Cotton", "Grassland/Pasture")
min_ak_class <- 0.9; min_ak_grid <- 0.5; min_cells <- 3; min_classes <- 3; n_boot <- 2000
set.seed(1)

# ---- class-day values: reference and products ---------------------------------
ref757 <- fread(file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv"))[
  crop_class %in% classes & avg_kernel >= min_ak_class, .(date = as.Date(date), crop_class, ref = sif757)]
g <- as.data.table(st_drop_geometry(st_read(file.path(proj_dir, "downscaled_sif757", "downscaled_sif757.gpkg"),
       query = "SELECT date, sif757, avg_kernel, crop_name FROM cells", quiet = TRUE)))
grid <- g[crop_name %in% classes & avg_kernel >= min_ak_grid,
          .(prod = median(sif757), n = .N), by = .(date = as.Date(date), crop_class = crop_name)][n >= min_cells, !"n"]
sn <- fread(file.path(proj_dir, "sifnet_comparison", "class_comparison.csv"))[
  crop_class %in% classes & avg_kernel >= min_ak_class & !is.na(sifnet_sif740),
  .(date = as.Date(date), crop_class, ref = sif740, prod = sifnet_sif740)]
sets <- list("Grid inversion (~800 m)" = merge(ref757, grid, by = c("date", "crop_class")),
             "SIFnet (0.005 deg)" = sn)

# ---- contrasts and the three measures -------------------------------------------
contrasts <- function(d) {
  d <- copy(d)[, k := .N, by = date][k >= min_classes]
  d[, `:=`(ref_c = ref - mean(ref), prod_c = prod - mean(prod)), by = date][]
}
measures <- function(d) d[, .(slope = cov(prod_c, ref_c) / var(ref_c), sd_ratio = sd(prod_c) / sd(ref_c),
                              r = cor(prod_c, ref_c))]
summarise <- function(d, product, set) {
  d <- contrasts(d); days <- unique(d$date)
  boot <- rbindlist(lapply(seq_len(n_boot), function(b) {
    s <- data.table(date = sample(days, replace = TRUE))[, rep := .I]
    measures(d[s, on = "date", allow.cartesian = TRUE])
  }))
  m <- measures(d)
  data.table(set, product, targets = length(days), class_days = nrow(d),
             slope = m$slope, slope_lo = quantile(boot$slope, .025), slope_hi = quantile(boot$slope, .975),
             sd_ratio = m$sd_ratio, sd_lo = quantile(boot$sd_ratio, .025), sd_hi = quantile(boot$sd_ratio, .975),
             r = m$r)
}
own <- rbindlist(lapply(names(sets), function(p) summarise(sets[[p]], p, "Each product, own targets")))

# common class-days: present in both the grid and SIFnet sets
key <- Reduce(function(a, b) merge(a, b, by = c("date", "crop_class")),
              lapply(sets, function(d) d[, .(date, crop_class)]))
common <- rbindlist(lapply(names(sets), function(p)
  summarise(merge(sets[[p]], key, by = c("date", "crop_class")), p, "Common class-days (2018-2021)")))
res <- rbind(own, common)
fwrite(res, file.path(out_dir, "contrast_amplitude.csv"))
print(res, digits = 2)

# ---- figure -----------------------------------------------------------------------
res[, product := factor(product, levels = rev(names(sets)))]
p <- ggplot(res, aes(y = product, x = slope, xmin = slope_lo, xmax = slope_hi, color = set)) +
  geom_vline(xintercept = c(0, 1), color = c("grey70", "grey40"), linetype = c(1, 2)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.35) +
  scale_color_manual(values = c("Each product, own targets" = "#2c6e8f", "Common class-days (2018-2021)" = "#d95f02")) +
  labs(x = "Share of the crop-class inversion's within-day crop contrast recovered (regression slope)",
       y = NULL, color = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "plot_contrast_amplitude.png"), p, width = 7, height = 3, dpi = 300, bg = "white")

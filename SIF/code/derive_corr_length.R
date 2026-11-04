# Derive the prior correlation length L of the grid inversion
# (downscale_sif_inversion.R: S_a,jk = sigma_a^2 exp(-d_jk / L), currently L = 1.5 km
# by assumption) from the crop mosaic itself.
#
# For each target day: a crop-driven SIF field on the grid,
#     s_j = sum_c F_jc x_c
# with F_jc the fraction of cell j in crop class c (30 m CDL, crop-year rule as in
# crop_classes.R) and x_c that day's crop-class SIF (downscale_sif_by_crop_cdl30.R).
# Its spatial autocorrelation rho(d) between cell centroids (binned by distance)
# is fitted with the inversion's own form, rho(d) = exp(-d / L). This weights the
# Figure 1 mosaic by the SIF contrasts actually present that day, so L can change
# with season (wheat vs. pasture in April, summer crops vs. stubble in August).
#
# Also reported: the SD of the crop-driven field, a cell-scale counterpart of the
# per-day sigma_a estimated from footprint variance in downscale_sif_inversion.R.
#
# Output (corr_length/): cell_class_fractions.csv, correlograms.csv,
#   daily_corr_length.csv, plot_corr_length.png, and a printed summary.

library(sf)
library(terra)
library(data.table)
library(exactextractr)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
cdl_dir  <- file.path(proj_dir, "cdl_30m")
grid_file <- file.path(proj_dir, "sub_road_grid.shp")
class_file <- file.path(proj_dir, "downscaled_sif757", "sif757_by_crop_cdl30.csv")
summary_file <- file.path(proj_dir, "downscaled_sif757", "daily_summary.csv")
out_dir <- file.path(proj_dir, "corr_length")
bin_km  <- 0.4          # distance bins for the correlogram
max_km  <- 8            # fit / plot range
current_L <- 1.5
dir.create(out_dir, showWarnings = FALSE)
source(file.path(proj_dir, "crop_classes.R"))   # class_levels, class_raster()

grid <- st_read(grid_file, quiet = TRUE)
grid$cell <- seq_len(nrow(grid))
xy <- st_coordinates(st_centroid(st_geometry(grid))) / 1000           # km (NAD83 / OK North)
D <- as.matrix(dist(xy))
ut <- upper.tri(D)
dvec <- D[ut]
bin <- floor(dvec / bin_km)
keep_pair <- dvec <= max_km

# ---- cell x class fractions for each crop map (crop year rule depends on date) ----
cls_dt <- fread(class_file)[, date := as.Date(date)]
days <- sort(unique(cls_dt$date))
map_key <- function(d) paste0(format(d, "%Y"), ifelse(format(d, "%m-%d") >= fall_start, "F", ""))
fractions <- list()
for (k in unique(vapply(days, map_key, ""))) {
  d0 <- days[vapply(days, map_key, "") == k][1]
  cr <- class_raster(d0)
  g <- st_transform(grid, crs(cr))
  fr <- exact_extract(cr, g, function(v, cov) {
    out <- numeric(length(class_levels))
    s <- tapply(cov, factor(v, levels = seq_along(class_levels)), sum)
    s[is.na(s)] <- 0
    s / sum(s)
  }, progress = FALSE)
  fractions[[k]] <- t(fr)                                            # cells x classes
}
fr_long <- rbindlist(lapply(names(fractions), function(k)
  data.table(map = k, cell = seq_len(nrow(fractions[[k]])), fractions[[k]])))
setnames(fr_long, c("map", "cell", class_levels))
fwrite(fr_long, file.path(out_dir, "cell_class_fractions.csv"))

# ---- per-day correlogram and fit ------------------------------------------------
fit_L <- function(cg) {
  cg <- cg[d <= max_km & is.finite(rho)]
  f <- tryCatch(nls(rho ~ exp(-d / L), data = cg, start = list(L = 1), weights = n,
                    lower = list(L = 0.05), upper = list(L = 50), algorithm = "port"),
                error = function(e) NULL)
  if (is.null(f)) NA_real_ else coef(f)[["L"]]
}
res <- list(); cgs <- list()
for (i in seq_along(days)) {
  day <- days[i]
  x <- cls_dt[date == day]
  xv <- setNames(rep(NA_real_, length(class_levels)), class_levels)
  xv[x$crop_class] <- x$sif757
  F <- fractions[[map_key(day)]]
  use <- !is.na(xv)
  s <- as.vector(F[, use, drop = FALSE] %*% xv[use] / rowSums(F[, use, drop = FALSE]))
  a <- s - mean(s, na.rm = TRUE)
  v <- mean(a^2, na.rm = TRUE)
  prod <- outer(a, a)[ut]
  cg <- data.table(bin = bin[keep_pair], p = prod[keep_pair])[, .(rho = mean(p, na.rm = TRUE) / v, n = .N), by = bin]
  cg[, d := (bin + 0.5) * bin_km]
  setorder(cg, d)
  cgs[[i]] <- cg[, date := day][]
  res[[i]] <- data.table(date = day, L_km = fit_L(cg), sd_crop_field = sqrt(v),
                         rho_neighbour = cg[d <= 1.2 & d > 0.4, weighted.mean(rho, n)])
}
cgs <- rbindlist(cgs); res <- rbindlist(res)
res[, doy := as.integer(format(date, "%j"))]
res[, season := c("Winter", "Spring", "Summer", "Autumn")[c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)[as.integer(format(date, "%m"))]]]
sa <- fread(summary_file)[, .(date = as.Date(date), sigma_a = prior_sd)]
res <- merge(res, sa, by = "date", all.x = TRUE)
fwrite(cgs, file.path(out_dir, "correlograms.csv"))
fwrite(res, file.path(out_dir, "daily_corr_length.csv"))

# pooled correlogram (pair-weighted over days) and its fit
pool <- cgs[, .(rho = weighted.mean(rho, n), n = sum(n)), by = d][order(d)]
L_pool <- fit_L(pool)

# ---- figure ---------------------------------------------------------------------
ex <- res[!is.na(L_km)][c(which.min(abs(doy - 100)), which.min(abs(doy - 220)))]   # an April and an August day
cg_ex <- cgs[date %in% ex$date][, lab := sprintf("%s (L = %.2f km)", format(date, "%d %b %Y"),
                                                  ex$L_km[match(date, ex$date)])]
dd <- seq(0, max_km, 0.05)
pa <- ggplot() +
  geom_point(data = pool, aes(d, rho), color = "black", size = 1.4) +
  geom_line(data = data.table(d = dd, rho = exp(-dd / L_pool)), aes(d, rho, linetype = "Fit to all days"), linewidth = 0.6) +
  geom_line(data = data.table(d = dd, rho = exp(-dd / current_L)), aes(d, rho, linetype = "Current prior (L = 1.5 km)"),
            color = "#d95f02", linewidth = 0.6) +
  geom_point(data = cg_ex, aes(d, rho, color = lab), size = 1, alpha = 0.8) +
  scale_linetype_manual(values = c("Fit to all days" = 1, "Current prior (L = 1.5 km)" = 2)) +
  scale_color_manual(values = c("#1b9e77", "#7570b3")) +
  geom_hline(yintercept = 0, color = "grey70") +
  labs(x = "Distance between cell centres (km)", y = "Autocorrelation of crop-driven SIF",
       linetype = NULL, color = NULL,
       title = sprintf("Correlogram of the crop-driven SIF field (all days: L = %.2f km)", L_pool)) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.box = "vertical")
pb <- ggplot(res[!is.na(L_km)], aes(doy, L_km, color = season)) +
  geom_hline(yintercept = current_L, linetype = 2, color = "#d95f02") +
  geom_hline(yintercept = L_pool, color = "black") +
  geom_point(size = 1.6) +
  scale_color_manual(values = c(Winter = "#7570b3", Spring = "#1b9e77", Summer = "#d95f02", Autumn = "#a6761d")) +
  labs(x = "Day of year", y = "Fitted correlation length L (km)", color = NULL,
       title = "Per-day L (dashed: current 1.5 km; solid: all-day fit)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "plot_corr_length.png"), pa + pb, width = 10, height = 4.8, dpi = 300)

# ---- summary --------------------------------------------------------------------
cat(sprintf("\n%d target days; pooled correlogram fit L = %.2f km (current prior %.1f km)\n", nrow(res), L_pool, current_L))
cat(sprintf("per-day L: median %.2f km, IQR %.2f-%.2f km\n", median(res$L_km, na.rm = TRUE),
            quantile(res$L_km, .25, na.rm = TRUE), quantile(res$L_km, .75, na.rm = TRUE)))
print(res[, .(days = .N, median_L = median(L_km, na.rm = TRUE), median_rho_neighbour = median(rho_neighbour),
              median_sd_crop_field = median(sd_crop_field), median_sigma_a = median(sigma_a, na.rm = TRUE)), by = season], digits = 3)
cat(sprintf("correlation between neighbouring cells (0.4-1.2 km), median: %.2f\n", median(res$rho_neighbour)))
cat(sprintf("crop-driven field SD vs sigma_a (per day): r = %.2f; median ratio sd/sigma_a = %.2f\n",
            cor(res$sd_crop_field, res$sigma_a, use = "complete"), median(res$sd_crop_field / res$sigma_a, na.rm = TRUE)))

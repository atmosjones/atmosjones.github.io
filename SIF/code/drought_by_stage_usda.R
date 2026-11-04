# Growth-stage drought test with stages from USDA NASS Crop Progress
# (crop_progress_calendar.R: year-specific active windows, independent of the
# satellite data). Grassland/pasture has no NASS stages and is left out of both
# groups. Same contrast and panels as plot_stage_contrast.R:
#   a  HLS NDVI, all clear scenes, active minus inactive per scene, monthly means
#   b  HLS NDVI on OCO-2 target days (nearest clear scene), target's DSCI
#   c  OCO-2 crop-class SIF, active minus inactive per target
# Output (drought_by_stage/): usda_stage_contrast_{ndvi_months,ndvi_targets,sif}.csv,
#   usda_stats.csv, plot_stage_contrast_usda.png / .pdf

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
dir_s <- file.path(proj_dir, "drought_by_stage")
source(file.path(proj_dir, "crop_progress_calendar.R"))   # windows, active_usda()
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
lev <- unlist(cls[, .(label_low, label_mid, label_high)])
dclass <- function(x) factor(fcase(x < cls$low_upper, lev[1], x < cls$high_lower, lev[2], default = lev[3]), levels = lev)
yrs <- range(windows$year)

contrast <- function(d, by) d[!is.na(active), .(
  n_active = sum(active), n_inactive = sum(!active),
  d_abs = if (any(active) && any(!active)) mean(residual[active]) - mean(residual[!active]) else NA_real_), by = by]

# ---- SIF -------------------------------------------------------------------------
s <- fread(file.path(proj_dir, "drought_residuals.csv"))[, date := as.Date(date)][year %between% yrs]
s[, active := active_usda(crop_class, date)]
sif <- merge(contrast(s, "date"), unique(s[, .(date, dsci, year, season)]), by = "date")[!is.na(d_abs)]
sif[, drought := dclass(dsci)]

# ---- NDVI, all scenes -> monthly means -------------------------------------------
n <- fread(file.path(proj_dir, "ndvi_comparison", "ndvi_residuals.csv"))[, date := as.Date(date)][year %between% yrs]
n[, active := active_usda(crop_class, date)]
nsc <- merge(contrast(n, "date"), unique(n[, .(date, dsci, year)]), by = "date")[!is.na(d_abs)]
nsc[, ym := format(date, "%Y-%m")]
nm <- nsc[, .(d_abs = mean(d_abs), dsci = mean(dsci), year = year[1]), by = ym]
nm[, drought := dclass(dsci)]

# ---- NDVI on target days (nearest clear scene) -------------------------------------
tg <- unique(s[, .(date, dsci, year)])
sc <- sort(unique(n$date))
tg[, hls_date := sc[vapply(date, function(d) which.min(abs(as.numeric(sc - d))), 1L)]]
nt <- merge(tg, nsc[, .(hls_date = date, d_abs)], by = "hls_date")
nt[, drought := dclass(dsci)]

fwrite(nm, file.path(dir_s, "usda_stage_contrast_ndvi_months.csv"))
fwrite(nt, file.path(dir_s, "usda_stage_contrast_ndvi_targets.csv"))
fwrite(sif, file.path(dir_s, "usda_stage_contrast_sif.csv"))

# ---- stats -----------------------------------------------------------------------
test <- function(d, label) {
  t <- suppressWarnings(cor.test(d$d_abs, d$dsci, method = "spearman"))
  z <- atanh(t$estimate); se <- 1 / sqrt(nrow(d) - 3)
  data.table(test = label, n = nrow(d), r = unname(t$estimate), lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
             p = t$p.value, kw_p = kruskal.test(d_abs ~ drought, d)$p.value,
             med_none = median(d[drought == lev[1], d_abs]), med_dry = median(d[drought == lev[2], d_abs]),
             med_drought = median(d[drought == lev[3], d_abs]))
}
common <- intersect(nt$date, sif$date)
st <- rbind(test(nm, "a. NDVI, all scenes (monthly means)"),
            test(nt, "b. NDVI on OCO-2 target days"),
            test(sif, "c. SIF"),
            test(nm[!year %in% 2022:2023], "a. without 2022-2023"),
            test(nt[!year %in% 2022:2023], "b. without 2022-2023"),
            test(sif[!year %in% 2022:2023], "c. without 2022-2023"),
            test(nt[date %in% common], "b. on targets common to b and c"),
            test(sif[date %in% common], "c. on targets common to b and c"))
fwrite(st, file.path(dir_s, "usda_stats.csv"))
cat("class-days classified (SIF):", s[!is.na(active), .N], "of", nrow(s), "; active:", s[active == TRUE, .N], "\n")
cat("targets with both stages: SIF", nrow(sif), "; NDVI targets", nrow(nt), "; common", length(common), "\n")
cat("by season (SIF targets):\n"); print(sif[, .N, by = season])
print(st, digits = 2)

# ---- figure ----------------------------------------------------------------------
for (d in list(nm, nt, sif)) d[, period := ifelse(year %in% 2022:2023, "2022-2023 drought", "Other years")]
pal <- setNames(c("#ffffb2", "#fecc5c", "#e31a1c"), lev)
stat_label <- function(r) sprintf("Spearman r = %.2f (95%% CI %.2f to %.2f)\np %s, n = %d", r$r, r$lo, r$hi,
                                  ifelse(r$p < 0.001, "< 0.001", sprintf("= %.2f", r$p)), r$n)
panel <- function(d, ylab, title, r) {
  n_by <- d[, .N, by = drought]
  ggplot(d, aes(drought, d_abs)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_boxplot(aes(fill = drought), outlier.shape = NA, width = 0.6) +
    geom_point(aes(shape = period), position = position_jitter(width = 0.12, seed = 1),
               size = 1.5, alpha = 0.75, fill = "grey15", color = "grey15") +
    annotate("text", x = 0.5, y = Inf, hjust = 0, vjust = 1.3, size = 3, label = stat_label(r)) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c("2022-2023 drought" = 24, "Other years" = 21), name = NULL) +
    scale_x_discrete(labels = function(x) { i <- match(x, n_by$drought); paste0(sub(" \\(", "\n(", x), "\nn = ", n_by$N[i]) }) +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
}
ylim_nd <- range(c(nm$d_abs, nt$d_abs)); ylim_nd <- ylim_nd + c(-0.04, 0.2) * diff(ylim_nd)
ylim_sf <- range(sif$d_abs); ylim_sf <- ylim_sf + c(-0.04, 0.2) * diff(ylim_sf)
pa <- panel(nm, "Active minus inactive crops (NDVI residual)", "HLS NDVI, all clear scenes (monthly means)", st[1]) +
  coord_cartesian(ylim = ylim_nd)
pb <- panel(nt, "Active minus inactive crops (NDVI residual)", "HLS NDVI on OCO-2 target days (nearest scene)", st[2]) +
  coord_cartesian(ylim = ylim_nd)
pc <- panel(sif, expression(atop("Active minus inactive crops", (SIF[757]~residual*","~mW~m^-2~sr^-1~nm^-1))),
            "OCO-2 crop-class SIF", st[3]) + coord_cartesian(ylim = ylim_sf)
saveRDS(list(ndvi_all = pa, ndvi_targets = pb, sif = pc), file.path(proj_dir, "fig5_panels_stage.rds"))   # for make_fig5.R
g <- pa + pb + pc + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ".",
                  caption = "Growth stage from USDA NASS Crop Progress (Oklahoma): active from 50% planted (wheat: jointing) to 50% mature; grassland not classified.") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 12, face = "bold"))
ggsave(file.path(dir_s, "plot_stage_contrast_usda.png"), g, width = 14.5, height = 5.6, dpi = 300, bg = "white")
ggsave(file.path(dir_s, "plot_stage_contrast_usda.pdf"), g, width = 14.5, height = 5.6, bg = "white")

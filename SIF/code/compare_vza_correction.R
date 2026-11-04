# Before/after the VZA normalisation (vza_correction.R): how the main results
# change. "Before" = results_before_vza_correction/ (copied just before the
# correction was added); "after" = the current outputs.
#
# Output: vza_correction_effect.csv (one row per quantity) and a printed summary.

library(data.table)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
old_dir  <- file.path(proj_dir, "results_before_vza_correction")
rd <- function(rel, which) fread(file.path(if (which == "old") old_dir else proj_dir, rel))
rows <- list()
add <- function(what, before, after) rows[[length(rows) + 1]] <<- data.table(quantity = what, before = before, after = after)

# ---- Per-cell inversion -----------------------------------------------------
for (w in c("old", "new")) assign(paste0("ds_", w), rd("downscaled_sif757/daily_summary.csv", w))
m <- merge(ds_old, ds_new, by = "date", suffixes = c("_old", "_new"))
add("days inverted", nrow(ds_old), nrow(ds_new))
add("per-cell: median day-mean observed SIF757", median(m$obs_mean_old), median(m$obs_mean_new))
add("per-cell: median DOFS", median(m$dofs_old), median(m$dofs_new))
add("per-cell: median cells with averaging kernel > 0.5", median(m$cells_constrained_old), median(m$cells_constrained_new))
add("per-cell: median reduced chi2", median(m$chi2_reduced_old), median(m$chi2_reduced_new))

# ---- Crop-class inversion (757) -----------------------------------------------
cc <- merge(rd("downscaled_sif757/sif757_by_crop_cdl30.csv", "old"),
            rd("downscaled_sif757/sif757_by_crop_cdl30.csv", "new"),
            by = c("date", "crop_class"), suffixes = c("_old", "_new"))
cc <- cc[avg_kernel_new >= 0.9 & avg_kernel_old >= 0.9 & !crop_class %in% c("Developed", "Other")]
cc[, `:=`(anom_old = sif757_old - mean(sif757_old), anom_new = sif757_new - mean(sif757_new)), by = date]
add("crop classes: median ratio after/before (class-days)", NA, median(cc$sif757_new / cc$sif757_old, na.rm = TRUE))
add("crop classes: r(before, after) of within-day contrasts", NA, cor(cc$anom_old, cc$anom_new))
add("crop classes: SD of within-day contrasts", sd(cc$anom_old), sd(cc$anom_new))
add("crop classes: slope of after-contrast on before-contrast", NA, coef(lm(anom_new ~ anom_old, cc))[2])
# ranking of classes by mean contrast
rk <- cc[, .(old = mean(anom_old), new = mean(anom_new)), by = crop_class]
add("crop classes: Spearman r of class ranking (mean contrast)", NA, cor(rk$old, rk$new, method = "spearman"))

# ---- SIFNet and GOSIF comparisons ---------------------------------------------
sn_o <- rd("sifnet_comparison/footprint_stats.csv", "old"); sn_n <- rd("sifnet_comparison/footprint_stats.csv", "new")
add("SIFNet: median footprint r (SIFNet)", median(sn_o$r_sifnet), median(sn_n$r_sifnet))
add("SIFNet: median footprint r (crop-class CV, 740)", median(sn_o$r_crop_cv), median(sn_n$r_crop_cv))
add("SIFNet: median SIFNet / OCO-2 ratio", median(sn_o$sifnet_ratio), median(sn_n$sifnet_ratio))
go_o <- rd("gosif_comparison/footprint_stats.csv", "old"); go_n <- rd("gosif_comparison/footprint_stats.csv", "new")
add("GOSIF: median footprint r (GOSIF)", median(go_o$r_gosif, na.rm = TRUE), median(go_n$r_gosif, na.rm = TRUE))
add("GOSIF: median footprint r (crop-class CV, 757)", median(go_o$r_crop_cv), median(go_n$r_crop_cv))
add("GOSIF: median GOSIF / OCO-2 ratio", median(go_o$gosif_ratio), median(go_n$gosif_ratio))
add("noise ceiling (757, median)", median(go_o$noise_ceiling), median(go_n$noise_ceiling))
gc_o <- rd("gosif_comparison/class_comparison.csv", "old"); gc_n <- rd("gosif_comparison/class_comparison.csv", "new")
f_cc <- function(d) { d <- d[avg_kernel >= 0.9 & !crop_class %in% c("Developed", "Other")]
  d[, `:=`(a = sif757 - mean(sif757), b = gosif_sif757 - mean(gosif_sif757)), by = date]; cor(d$a, d$b) }
add("GOSIF: r of within-day crop contrasts", f_cc(gc_o), f_cc(gc_n))

out <- rbindlist(rows)
fwrite(out, file.path(proj_dir, "vza_correction_effect.csv"))
print(out, digits = 3)

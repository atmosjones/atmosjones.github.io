# Table: OCO-2 target-mode observations over the study area, by year.
#
#   targets           days with target-mode soundings in raster_extent.shp
#                     (one overpass per target)
#   soundings         all target-mode soundings
#   quality_soundings Quality_Flag 0 (best) or 1 (good), as kept by the inversions
#   pct_quality       quality_soundings / soundings
#   targets_used      days inverted by downscale_sif_inversion.R (>= min_obs usable
#                     footprints after the quality and 50% grid-coverage filters)
#   soundings_used    footprints in those inversions
#   lite_version      OCO-2 Lite SIF product version(s) that year
#
# Output: table_target_counts.csv and table_target_counts.tex (booktabs) in the
# project folder.

library(sf)
library(data.table)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
oco_file <- file.path(proj_dir, "oco2_sif_raster_extent.gpkg")
summary_file <- file.path(proj_dir, "downscaled_sif757", "daily_summary.csv")
quality_flags <- c(0, 1)   # as in downscale_sif_inversion.R

o <- as.data.table(st_read(oco_file, quiet = TRUE, query = paste(
  "SELECT time_utc, Quality_Flag, lite_build FROM oco2_sif WHERE Metadata_MeasurementMode = 2")))
o[, `:=`(date = as.Date(substr(time_utc, 1, 10)), year = as.integer(substr(time_utc, 1, 4)))]
o[, version := fifelse(grepl("^B110", lite_build), "V11r", "V11.2r")]

by_year <- o[, .(targets = uniqueN(date),
                 soundings = .N,
                 quality_soundings = sum(Quality_Flag %in% quality_flags),
                 lite_version = paste(sort(unique(version), decreasing = TRUE), collapse = ", ")),
             by = year]
used <- fread(summary_file)[, .(year = as.integer(substr(date, 1, 4)), n_obs)]
used <- used[, .(targets_used = .N, soundings_used = sum(n_obs)), by = year]
tab <- merge(by_year, used, by = "year", all.x = TRUE)
tab[is.na(targets_used), `:=`(targets_used = 0L, soundings_used = 0L)]
tab[, pct_quality := round(100 * quality_soundings / soundings, 1)]
setcolorder(tab, c("year", "targets", "soundings", "quality_soundings", "pct_quality",
                   "targets_used", "soundings_used", "lite_version"))
setorder(tab, year)
tot <- tab[, .(year = NA_integer_, targets = sum(targets), soundings = sum(soundings),
               quality_soundings = sum(quality_soundings),
               pct_quality = round(100 * sum(quality_soundings) / sum(soundings), 1),
               targets_used = sum(targets_used), soundings_used = sum(soundings_used),
               lite_version = "")]
out <- rbind(tab, tot)
fwrite(out[, year := ifelse(is.na(year), "Total", as.character(year))],
       file.path(proj_dir, "table_target_counts.csv"))

# LaTeX (booktabs)
fmt <- function(x) formatC(x, format = "d", big.mark = ",")
first_partial <- format(min(o$date), "%d %b %Y"); last_partial <- format(max(o$date), "%d %b %Y")
rows <- out[, sprintf("%s & %s & %s & %s & %.1f & %s & %s \\\\",
                      year, fmt(targets), fmt(soundings), fmt(quality_soundings), pct_quality,
                      fmt(targets_used), fmt(soundings_used))]
rows[length(rows)] <- paste("\\midrule\n", rows[length(rows)])
tex <- c("\\begin{table}[ht]", "\\centering",
         sprintf("\\caption{OCO-2 target-mode observations over the study area by year (%s to %s). Quality soundings have \\texttt{Quality\\_Flag} 0 or 1; targets and soundings used are those in the downscaling inversion (days with at least 50 quality soundings, each at least half inside the grid). OCO-2 Lite SIF V11r through March 2024, V11.2r after.}",
                 first_partial, last_partial),
         "\\label{tab:targets}",
         "\\begin{tabular}{lrrrrrr}", "\\toprule",
         "Year & Targets & Soundings & Quality soundings & Quality (\\%) & Targets used & Soundings used \\\\",
         "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(proj_dir, "table_target_counts.tex"))
print(out)

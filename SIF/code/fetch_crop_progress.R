# Download USDA NASS Crop Progress for Oklahoma (state level, weekly), 2014-2025,
# for the crops in the stage test: winter wheat, corn, sorghum, soybeans, cotton.
# Uses the Quick Stats API; the key is read from ~/.nass_api_key and never
# written anywhere else or printed.
# Output: crop_progress/ok_crop_progress_raw.csv  (PROGRESS: stage percentages)
#         crop_progress/ok_crop_condition_raw.csv  (CONDITION: % very poor ... excellent)

library(data.table)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "crop_progress")
dir.create(out_dir, showWarnings = FALSE)
key <- trimws(readLines(path.expand("~/.nass_api_key"), n = 1, warn = FALSE))
years <- 2014:2025
commodities <- c("WHEAT", "CORN", "SORGHUM", "SOYBEANS", "COTTON")

get_one <- function(commodity, statistic = "PROGRESS") {
  q <- c(key = key, source_desc = "SURVEY", sector_desc = "CROPS", statisticcat_desc = statistic,
         commodity_desc = commodity, state_alpha = "OK", agg_level_desc = "STATE",
         freq_desc = "WEEKLY", year__GE = min(years), year__LE = max(years), format = "CSV")
  url <- paste0("https://quickstats.nass.usda.gov/api/api_GET/?",
                paste(names(q), vapply(q, URLencode, "", reserved = TRUE), sep = "=", collapse = "&"))
  f <- tempfile(fileext = ".csv")
  ok <- tryCatch({ download.file(url, f, quiet = TRUE); TRUE }, error = function(e) FALSE)
  if (!ok) stop("download failed for ", commodity, " (check the key or the service)")   # url not printed: it holds the key
  d <- fread(f, colClasses = "character")
  unlink(f)
  d
}
raw <- rbindlist(lapply(commodities, get_one), fill = TRUE)
fwrite(raw, file.path(out_dir, "ok_crop_progress_raw.csv"))
cat(nrow(raw), "progress rows\n")
cond <- rbindlist(lapply(commodities, get_one, statistic = "CONDITION"), fill = TRUE)
fwrite(cond, file.path(out_dir, "ok_crop_condition_raw.csv"))
cat(nrow(cond), "condition rows\n")
print(cond[, .(n = .N, first = min(week_ending), last = max(week_ending)), keyby = .(commodity_desc, short_desc)])

# Year-specific active windows for each crop class from USDA NASS Crop Progress
# (Oklahoma, state level; fetch_crop_progress.R). A crop is active from the date
# 50% of the state crop reached the start stage to the date 50% reached the end
# stage (linear between weekly reports):
#   Winter Wheat               jointing  -> harvested        (fall/winter growth = inactive)
#   Corn, Sorghum              planted   -> mature
#   Soybeans                   planted   -> dropping leaves
#   Cotton                     planted   -> bolls opening
#   Dbl Crop WinWht/Soybeans   wheat jointing -> wheat harvested, then
#                              wheat harvested -> soybean dropping leaves
#   Grassland/Pasture          no NASS stage data: not classified
# Where a stage is missing in a year, its median day of year over 2014-2025 is used.
# Output: crop_progress/active_windows.csv (crop_class, year, start, end, filled)
#         crop_progress/stage_dates_50pct.csv
# active_usda(crop_class, date) is defined for other scripts to source.

library(data.table)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
cp_dir <- file.path(proj_dir, "crop_progress")
years <- 2014:2025

r <- fread(file.path(cp_dir, "ok_crop_progress_raw.csv"))
r[, `:=`(stage = sub(".*MEASURED IN PCT ", "", short_desc), date = as.Date(week_ending),
         pct = as.numeric(gsub(",", "", Value)))]
r <- r[!is.na(pct) & !is.na(date)]
r[, cal_year := as.integer(format(date, "%Y"))]           # calendar year of the report week
x50 <- r[, { o <- order(date); d <- date[o]; p <- pct[o]; i <- which(p >= 50)[1]
  if (is.na(i)) .(d50 = as.Date(NA))
  else if (i == 1) .(d50 = d[1])
  else .(d50 = d[i - 1] + (50 - p[i - 1]) / (p[i] - p[i - 1]) * as.numeric(d[i] - d[i - 1])) },
  by = .(commodity = commodity_desc, stage, cal_year)]
x50 <- x50[!is.na(d50) & cal_year %in% years]
fwrite(x50, file.path(cp_dir, "stage_dates_50pct.csv"))

stage_date <- function(com, stg, y) {
  v <- x50[commodity == com & stage == stg & cal_year == y, d50]
  if (length(v)) return(list(d = v[1], filled = FALSE))
  med <- x50[commodity == com & stage == stg, median(as.numeric(format(d50, "%j")))]
  list(d = as.Date(sprintf("%d-01-01", y)) + round(med) - 1, filled = TRUE)
}
spec <- list(
  "Winter Wheat" = list(c("WHEAT", "JOINTING", "WHEAT", "HARVESTED")),
  "Corn"         = list(c("CORN", "PLANTED", "CORN", "MATURE")),
  "Sorghum"      = list(c("SORGHUM", "PLANTED", "SORGHUM", "MATURE")),
  "Soybeans"     = list(c("SOYBEANS", "PLANTED", "SOYBEANS", "DROPPING LEAVES")),
  "Cotton"       = list(c("COTTON", "PLANTED", "COTTON", "BOLLS OPENING")),
  "Dbl Crop WinWht/Soybeans" = list(c("WHEAT", "JOINTING", "WHEAT", "HARVESTED"),
                                    c("WHEAT", "HARVESTED", "SOYBEANS", "DROPPING LEAVES")))
windows <- rbindlist(lapply(names(spec), function(cc) rbindlist(lapply(years, function(y)
  rbindlist(lapply(spec[[cc]], function(s) {
    a <- stage_date(s[1], s[2], y); b <- stage_date(s[3], s[4], y)
    data.table(crop_class = cc, year = y, start = a$d, end = b$d, filled = a$filled | b$filled)
  }))))))
fwrite(windows, file.path(cp_dir, "active_windows.csv"))

# TRUE/FALSE for classified crops, NA for unclassified classes (grassland, others)
active_usda <- function(cc, date) {
  vapply(seq_along(date), function(i) {
    ww <- windows[crop_class == cc[i]]
    if (!nrow(ww)) return(NA)
    any(date[i] >= ww$start & date[i] <= ww$end)
  }, logical(1))
}

if (sys.nframe() == 0) {   # run directly: summary
  cat(sprintf("filled from the median: %d of %d windows\n", windows[filled == TRUE, .N], nrow(windows)))
  print(windows[filled == TRUE])
  f <- function(v) format(as.Date(round(median(as.numeric(format(v, "%j")))) - 1, origin = "2021-01-01"), "%d %b")
  print(windows[, .(median_start = f(start), median_end = f(end), median_days = median(as.numeric(end - start))),
                by = .(crop_class, part = ifelse(crop_class == "Dbl Crop WinWht/Soybeans" & format(start, "%m") > "05", "soy", "main"))])
}

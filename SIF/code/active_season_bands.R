# Median active season of each crop class over 2014-2025, as day-of-year bands
# for background shading in seasonal-cycle plots (make_fig3.R, make_fig4.R).
# From crop_progress/active_windows.csv (crop_progress_calendar.R: USDA NASS Crop
# Progress, Oklahoma; 50% planted, or jointing for winter wheat, to 50% mature).
# Double-cropped wheat/soybeans has two bands (wheat, then soybeans); classes
# without NASS stages (grassland/pasture) have none.

active_season_bands <- function(proj_dir = path.expand("~/Documents/cropland_downscaling_paper")) {
  w <- data.table::fread(file.path(proj_dir, "crop_progress", "active_windows.csv"))
  w[, `:=`(start = as.Date(start), end = as.Date(end))]
  w[, part := ifelse(crop_class == "Dbl Crop WinWht/Soybeans" & as.integer(format(start, "%m")) >= 6, 2L, 1L)]
  w[, .(xmin = median(as.numeric(format(start, "%j"))), xmax = median(as.numeric(format(end, "%j"))),
        n_years = .N), by = .(crop_class, part)]
}

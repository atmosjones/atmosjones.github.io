# Viewing-zenith-angle (VZA) normalisation of OCO-2 target-mode SIF to nadir.
#
# Within each target day, retrieved SIF rises with VZA (~10% of the day mean per
# 10 deg at 757 nm; vza_check.R): the same at both sides of the pass and within
# the same grid cell, so it follows VZA itself, not the surface seen or the
# sun-view geometry. Whether it is a canopy directional effect or a retrieval
# artefact can't be told from these data, so it is treated as a normalisation:
# every sounding is put on a nadir-view footing,
#
#   SIF_nadir = SIF / f(VZA),   f(VZA) = 1 + a1 VZA + a2 VZA^2  (VZA in degrees),  f(0) = 1
#
# and its uncertainty is divided by the same factor. Coefficients: weighted
# within-day fit of SIF / day-mean SIF on VZA (target mode, Quality_Flag 0/1,
# days with day-mean SIF > 0.15), all soundings 2014-2026. Leaving each year out,
# the 757 nm curve removes ~3/4 of the within-day VZA trend in the held-out year
# (remaining median 1.5% per 10 deg). The 740 nm band is fitted separately.
#
#   VZA      0     20     40     60     70
#   f_757  1.00   1.09   1.27   1.55   1.72
#   f_740  1.00   1.06   1.27   1.65   1.90
#
# Set apply_vza_correction <- FALSE before sourcing to switch it off everywhere.

if (!exists("apply_vza_correction")) apply_vza_correction <- TRUE

vza_coef <- list(                       # per degree and per degree squared
  "757" = c(a1 = 2.212e-3, a2 = 1.159e-4),
  "740" = c(a1 = -1.219e-3, a2 = 2.016e-4))

vza_factor <- function(vza, band = "757") {
  a <- vza_coef[[band]]
  1 + a[["a1"]] * vza + a[["a2"]] * vza^2
}

# Divide the SIF and uncertainty columns present in `oco` (data.frame / sf from
# oco2_sif_raster_extent.gpkg; must include VZA) by the band's factor.
correct_vza <- function(oco) {
  if (!apply_vza_correction) return(oco)
  if (is.null(oco$VZA)) stop("correct_vza() needs the VZA column")
  cols <- list("757" = c("Science_SIF_757nm", "Science_SIF_Uncertainty_757nm"),
               "740" = c("SIF_740nm", "SIF_Uncertainty_740nm"))
  for (band in names(cols)) {
    f <- vza_factor(oco$VZA, band)
    for (cl in cols[[band]]) if (!is.null(oco[[cl]])) oco[[cl]] <- oco[[cl]] / f
  }
  oco
}

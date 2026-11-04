# Sensitivity of the grid inversion to the footprint's spatial response.
# downscale_sif_inversion.R treats each sounding as uniformly sensitive across its
# footprint polygon (row of K = fraction of the polygon in each cell). Here the
# response is instead a 2-D Gaussian centred on the polygon centroid, with axes
# along the polygon's principal axes (from its exact second moments), in two widths:
#   gauss_var   sigma_i = L_i / sqrt(12): same second moment as the uniform polygon
#               (isolates the effect of shape: peaked vs flat, same overall size)
#   gauss_fwhm  sigma_i = L_i / 2.355: FWHM = footprint length/width (half
#               sensitivity at the edges; ~40% of the weight outside the polygon)
# where L_i = sqrt(12 * lambda_i) are the equivalent-rectangle side lengths. Rows
# are the Gaussian integrated over grid cells on a 25-m raster, renormalized over
# the grid; soundings with < 50% of their weight in the grid are dropped (as before).
# Everything else (prior, sigma_a, settings) is the inversion's own code.
# Output (footprint_response/): test_days.csv, cells_<date>.csv, printed comparison.

library(sf); library(terra); library(data.table); library(Matrix)
proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
out_dir <- file.path(proj_dir, "footprint_response"); dir.create(out_dir, showWarnings = FALSE)
args <- commandArgs(TRUE)
test_days <- if (length(args)) as.Date(args) else
  as.Date(c("2019-04-19", "2015-01-25", "2020-03-11", "2015-08-21", "2025-11-22", "2017-12-09"))
px <- 25                                                     # integration raster (m)

# ---- the inversion's own setup and functions (everything before "# ---- Run") ----
src <- readLines(file.path(proj_dir, "downscale_sif_inversion.R"))
src <- src[seq_len(grep("^# ---- Run", src) - 1)]
src <- sub("^make_maps .*", "make_maps <- FALSE", src)
src <- sub("^save_A_days .*", 'save_A_days <- as.Date("1900-01-01")', src)
eval(parse(text = src))
build_uniform <- build_jacobian

# ---- grid cell id raster -------------------------------------------------------
gv <- vect(grid)
r0 <- rast(ext(gv) + 4000, resolution = px, crs = crs(gv))
idr <- rasterize(gv, r0, field = "cell_id")
idm <- as.matrix(idr, wide = TRUE)                           # rows from top
e <- ext(idr); x0 <- e[1]; y1 <- e[4]; nr <- nrow(idm); nc <- ncol(idm)

# exact polygon moments -> centroid, principal axes, equivalent side lengths
poly_axes <- function(g) {
  v <- st_coordinates(g)[, 1:2]; v <- v[-nrow(v), , drop = FALSE]
  x <- v[, 1]; y <- v[, 2]; x2 <- c(x[-1], x[1]); y2 <- c(y[-1], y[1])
  cr <- x * y2 - x2 * y; A <- sum(cr) / 2
  cx <- sum((x + x2) * cr) / (6 * A); cy <- sum((y + y2) * cr) / (6 * A)
  Ixx <- sum((y^2 + y * y2 + y2^2) * cr) / 12 / A - cy^2     # per-area second moments about centroid
  Iyy <- sum((x^2 + x * x2 + x2^2) * cr) / 12 / A - cx^2
  Ixy <- sum((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr) / 24 / A - cx * cy
  ev <- eigen(matrix(c(Iyy, Ixy, Ixy, Ixx), 2), symmetric = TRUE)   # cov of (x, y)
  list(c = c(cx, cy), vec = ev$vectors, L = sqrt(12 * pmax(ev$values, 1)))
}

make_builder <- function(width) function(fp) {
  fp$obs_id <- seq_len(nrow(fp))
  rows <- vector("list", nrow(fp)); cov <- numeric(nrow(fp))
  for (i in seq_len(nrow(fp))) {
    a <- poly_axes(st_geometry(fp)[[i]])
    s <- if (width == "var") a$L / sqrt(12) else a$L / 2.3548
    R <- 3.5 * max(s)
    c0 <- max(1, floor((a$c[1] - R - x0) / px)); c1 <- min(nc, ceiling((a$c[1] + R - x0) / px))
    r1 <- min(nr, ceiling((y1 - (a$c[2] - R)) / px)); r0_ <- max(1, floor((y1 - (a$c[2] + R)) / px))
    if (c1 < c0 || r1 < r0_) next
    xs <- x0 + (c0:c1 - 0.5) * px; ys <- y1 - (r0_:r1 - 0.5) * px
    dx <- outer(rep(1, length(ys)), xs - a$c[1]); dy <- outer(ys - a$c[2], rep(1, length(xs)))
    u <- dx * a$vec[1, 1] + dy * a$vec[2, 1]; w <- dx * a$vec[1, 2] + dy * a$vec[2, 2]
    g <- exp(-0.5 * ((u / s[1])^2 + (w / s[2])^2))
    ids <- idm[r0_:r1, c0:c1]
    tot <- (2 * pi * s[1] * s[2]) / px^2                         # full Gaussian integral, in pixels
    inside <- !is.na(ids)
    cov[i] <- sum(g[inside]) / tot
    if (any(inside)) rows[[i]] <- data.table(obs_id = i, cell_id = ids[inside], w = g[inside])[, .(w = sum(w)), by = .(obs_id, cell_id)]
  }
  keep <- which(cov >= min_coverage)
  p <- rbindlist(rows)[obs_id %in% keep]
  K <- sparseMatrix(i = match(p$obs_id, keep), j = p$cell_id, x = p$w, dims = c(length(keep), n_cells))
  K <- Diagonal(x = 1 / rowSums(K)) %*% K
  list(K = K, keep = keep, coverage = cov[keep])
}

run <- function(d, builder, label) {
  build_jacobian <<- builder
  t0 <- Sys.time(); r <- invert_day(d)
  r$cells[, `:=`(method = label, secs = as.numeric(Sys.time() - t0, units = "secs"))]
  r
}
res <- list(); cells <- list()
for (k in seq_along(test_days)) {
  d <- test_days[k]
  for (m in list(list("uniform", build_uniform), list("gauss_var", make_builder("var")), list("gauss_fwhm", make_builder("fwhm")))) {
    r <- run(d, m[[2]], m[[1]])
    c <- r$cells; ok <- c$avg_kernel > 0.5
    interior_area <- sum(cell_area[c$res_area_pr_km2 < 2.9 & c$ak_row_sum > 0.5])
    res[[length(res) + 1]] <- data.table(date = d, method = m[[1]], n_obs = r$summary$n_obs, dofs = r$summary$dofs,
      cells_ak05 = sum(ok), res_med_constrained = median(c$res_area_pr_km2[ok]),
      area_finer_2.9 = r$summary$area_finer_than_footprint, prior_sd = r$summary$prior_sd,
      chi2 = r$summary$chi2_reduced, secs = c$secs[1])
    cells[[length(cells) + 1]] <- c[, .(date, cell_id, method, sif757, avg_kernel, res_area_pr_km2, ak_row_sum)]
    cat(sprintf("%s %-10s DOFS %6.1f  AK>0.5 %4d  res %.2f km2  (%.0f s)\n", d, m[[1]], r$summary$dofs, sum(ok),
                median(c$res_area_pr_km2[ok]), c$secs[1]))
  }
}
res <- rbindlist(res); cells <- rbindlist(cells)
fwrite(res, file.path(out_dir, "test_days.csv")); fwrite(cells, file.path(out_dir, "test_cells.csv"))

# ---- comparison against uniform ------------------------------------------------------
w <- dcast(cells, date + cell_id ~ method, value.var = c("sif757", "avg_kernel"))
cmp <- w[, .(r_sif_var = cor(sif757_uniform, sif757_gauss_var), r_sif_fwhm = cor(sif757_uniform, sif757_gauss_fwhm),
             rms_diff_var = sqrt(mean((sif757_gauss_var - sif757_uniform)^2)),
             rms_diff_fwhm = sqrt(mean((sif757_gauss_fwhm - sif757_uniform)^2)),
             sd_sif_uniform = sd(sif757_uniform)), by = date]
cat("\nPosterior SIF per cell vs uniform:\n"); print(cmp, digits = 2)
cat("\nSummary by method (medians over test days):\n")
print(res[, .(dofs = median(dofs), cells_ak05 = median(cells_ak05), res_constrained = median(res_med_constrained),
              area_finer_2.9 = median(area_finer_2.9), chi2 = median(chi2)), by = method], digits = 3)
rel <- dcast(res, date ~ method, value.var = "dofs")
cat("\nDOFS relative to uniform:\n"); print(rel[, .(date, var = gauss_var / uniform, fwhm = gauss_fwhm / uniform)], digits = 3)

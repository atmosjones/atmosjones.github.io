# Figure 2 for the manuscript: downscaling geometry and effective resolution on one
# target day (panels built and saved by make_fig1.R in fig1_panels_inversion.rds).
#   a  footprints and their SIF          b  one footprint = one row of A (fraction of the footprint in each cell)
#   c  overlapping soundings per cell    d  posterior SIF per cell
#   e  averaging-kernel diagonal         f  effective resolution per cell
# Restyled here: no panel titles, no equation boxes in b, larger text (as in
# Figures 1, 3, 4), and horizontal colour bars below each map so the maps are wider.
# Output: fig2_manuscript.png (600 dpi, at printed size)

library(ggplot2)
library(patchwork)
source(file.path(path.expand("~/Documents/cropland_downscaling_paper"), "fig_fonts.R"))   # fig_font

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
p <- readRDS(file.path(proj_dir, "fig1_panels_inversion.rds"))
# Saved at the printed size (Word text width 6.5 in), so these sizes are the sizes in the
# manuscript: 9.5 pt titles, 8.5 pt tick and in-map labels, in Montserrat.
fig_w <- 6.5; fig_h <- 5.2
fs <- 9.5

restyle <- function(g, legend = TRUE) {
  g <- g + labs(title = NULL)
  # text drawn in the map (scale-bar labels, cell fractions): scale up with the base size
  for (i in seq_along(g$layers)) {
    L <- g$layers[[i]]
    if (inherits(L$geom, c("GeomText", "GeomLabel")) && !is.null(L$aes_params$size)) {
      L$aes_params$size <- (fs - 1) / .pt
      L$aes_params$family <- fig_font
      # scale bar: keep only the "5 km" label (the "0" collides with it at print size)
      if (!is.null(L$aes_params$label) && any(L$aes_params$label == "5 km")) {
        # one label, to the right of the bar and centred on it (bar: y0 .. y0 + 350 m; labels sat at y0 + 950)
        L$aes_params$label <- c("", "5 km")
        L$data$x <- max(L$data$x) + 400
        L$data$y <- L$data$y - 775
        L$aes_params$hjust <- 0; L$aes_params$vjust <- 0.5
      }
      if (inherits(L$data, "sf")) L$aes_params$size <- 7 / .pt                             # cell fractions in b
    }
    g$layers[[i]] <- L
  }
  g + theme(text = element_text(size = fs, family = fig_font), panel.border = element_blank(),
            legend.position = if (legend) "bottom" else "none", legend.direction = "horizontal",
            legend.title = element_text(size = fs), legend.text = element_text(size = fs - 1),
            legend.title.position = "top", legend.key.width = unit(1, "null"), legend.key.height = unit(0.22, "cm"),
            legend.margin = margin(1, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
            plot.margin = margin(2, 7, 2, 7))
}
# b: drop the two text boxes (equation and footprint note); keep the cell fractions
p$d$layers <- Filter(function(l) !inherits(l$geom, "GeomLabel"), p$d$layers)
# a: show its SIF colour bar too (same scale as d)
sif_title <- expression(SIF[757] ~ (mW ~ m^-2 ~ sr^-1 ~ nm^-1))
# rename a panel's colour/fill scale (the scales carry their own names, which override labs())
retitle <- function(g, aes_name, title, ...) {
  sc <- g$scales$get_scales(aes_name)$clone()
  sc$name <- title
  for (k in names(list(...))) sc[[k]] <- list(...)[[k]]
  g$scales$scales <- c(Filter(function(x) !aes_name %in% x$aesthetics, g$scales$scales), list(sc))
  g
}
pa <- restyle(retitle(p$c, "colour", sif_title) + guides(colour = guide_colourbar()))
pb <- restyle(retitle(p$d, "fill", "Fraction of footprint", limits = c(0, 0.25),
                      breaks = c(0, 0.1, 0.2), labels = c("0", "0.1", "0.2")) + guides(fill = guide_colourbar()))
pc <- restyle(retitle(p$g, "fill", "Overlapping soundings", breaks = c(0, 50, 100)))
pd <- restyle(retitle(p$e, "fill", sif_title))
pe <- restyle(retitle(p$f, "fill", "Averaging-kernel diagonal", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")))
pf <- restyle(retitle(p$h, "fill", expression(Effective ~ resolution ~ (km^2)),
                      breaks = c(0.64, 2.9, 10), labels = c("0.64", "2.9", "10")))

fig <- ((pa | pb | pc) + plot_layout(widths = c(1, 1, 1))) / ((pd | pe | pf) + plot_layout(widths = c(1, 1, 1))) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") &
  theme(plot.tag = element_text(size = fs + 0.5, face = "bold", family = fig_font))
ggsave(file.path(proj_dir, "fig2_manuscript.png"), fig, width = fig_w, height = fig_h, dpi = 600, bg = "white",
       device = ragg::agg_png)

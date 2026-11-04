# Colour the facet strips of a facet_wrap ggplot by crop class, so the strips
# match the crop colours used throughout the paper. Returns a gtable; wrap it
# with patchwork::wrap_elements(full = ...) to assemble it with other panels.
#
#   p       ggplot with facet_wrap(~crop_class, ncol = ncol)
#   colors  named vector: crop class -> fill colour
#   labels  facet labels in panel order (levels of the facet variable present)
# Text is black on light fills and white on dark ones.

library(grid)

color_strips <- function(p, colors, labels, ncol) {
  g <- ggplot2::ggplotGrob(p)
  text_col <- function(fill) {
    rgb <- grDevices::col2rgb(fill) / 255
    if (sum(c(0.299, 0.587, 0.114) * rgb) > 0.6) "black" else "white"
  }
  recolor <- function(x, fill, col) {
    set_gp <- function(x, ...) {            # keep gp a gpar object even when it starts NULL
      x$gp <- do.call(gpar, utils::modifyList(unclass(x$gp %||% gpar()), list(...)))
      x
    }
    if (inherits(x, "rect")) x <- set_gp(x, fill = fill)
    if (inherits(x, "text")) x <- set_gp(x, col = col)
    if (!is.null(x$children)) {            # keep names: grid draws children by childrenOrder
      kids <- lapply(x$children, recolor, fill, col)
      names(kids) <- names(x$children)
      class(kids) <- "gList"
      x$children <- kids
    }
    if (!is.null(x$grobs)) x$grobs <- lapply(x$grobs, recolor, fill, col)
    x
  }
  for (k in seq_along(labels)) {
    nm <- sprintf("strip-t-%d-%d", (k - 1) %% ncol + 1, (k - 1) %/% ncol + 1)
    i <- which(g$layout$name == nm)
    if (!length(i)) next
    fill <- colors[[labels[k]]]
    g$grobs[[i]] <- recolor(g$grobs[[i]], fill, text_col(fill))
  }
  g
}
`%||%` <- function(a, b) if (is.null(a)) b else a

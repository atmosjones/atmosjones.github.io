# Register the figure font (Google Fonts, SIL OFL; files in fonts/) for the figure scripts.
# PNG via ragg finds it through systemfonts. (R's cairo_pdf does not: PDFs would need the
# fonts installed in ~/Library/Fonts.) Montserrat is also in fonts/ if wanted again.
library(systemfonts)
fig_font <- "Barlow"
if (!fig_font %in% system_fonts()$family) local({            # use the installed font if there is one
  d <- file.path(path.expand("~/Documents/cropland_downscaling_paper"), "fonts")
  register_font(fig_font, plain = file.path(d, paste0(fig_font, "-400.ttf")), bold = file.path(d, paste0(fig_font, "-700.ttf")),
                italic = file.path(d, paste0(fig_font, "-400-italic.ttf")), bolditalic = file.path(d, paste0(fig_font, "-700.ttf")))
})

# Shared settings for figures saved at their printed width (Word text width 6.5 in)
fig_width <- 6.5
fs_title <- 9.5; fs_tick <- 8.5                                   # pt
ggplot2::update_geom_defaults("text", list(family = fig_font))
ggplot2::update_geom_defaults("label", list(family = fig_font))
font_theme <- ggplot2::theme(text = ggplot2::element_text(family = fig_font))
save_fig <- function(fig, stem, height, width = fig_width) {    # PNG, 600 dpi (PDF devices cannot use fig_font here)
  ggplot2::ggsave(file.path(path.expand("~/Documents/cropland_downscaling_paper"), paste0(stem, ".png")), fig,
                  width = width, height = height, dpi = 600, bg = "white", device = ragg::agg_png)
}

# Short crop-class names for figure text at print size (same in every figure)
crop_short <- c("Winter Wheat" = "Winter wheat", "Dbl Crop WinWht/Soybeans" = "Wheat/soy (double)",
                "Soybeans" = "Soybeans", "Corn" = "Corn", "Sorghum" = "Sorghum", "Cotton" = "Cotton",
                "Grassland/Pasture" = "Grassland/pasture", "Fallow/Idle Cropland" = "Fallow/idle",
                "Trees/Wetlands" = "Trees/wetlands", "Developed" = "Developed", "Other" = "Other")

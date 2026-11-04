# Distribution of the Drought Severity and Coverage Index (DSCI) on the target
# days, against the full weekly record, to judge how to bin drought severity.
#   a  DSCI on the 104 crop-class target days, stacked by season, with the
#      DSCI class cut-offs used in drought_analysis.R (drought_dsci_classes.csv)
#   b  every weekly USDM map, 2014-2025 (drought_index_study_area.csv)
# Output: plot_dsci_histogram.png and a printed summary.

library(data.table)
library(ggplot2)
library(patchwork)

proj_dir <- path.expand("~/Documents/cropland_downscaling_paper")
seasons <- c("Spring", "Summer", "Autumn", "Winter")
season_cols <- c(Spring = "#1b9e77", Summer = "#d95f02", Autumn = "#a6761d", Winter = "#7570b3")
tg <- unique(fread(file.path(proj_dir, "drought_residuals.csv"))[, .(date, dsci, season, year)])
tg[, season := factor(season, levels = seasons)]
wk <- fread(file.path(proj_dir, "drought_index_study_area.csv"))[as.Date(date) >= as.Date("2014-01-01") &
                                                                  as.Date(date) <= as.Date("2025-12-31")]
cls <- fread(file.path(proj_dir, "drought_dsci_classes.csv"))
cuts <- c(cls$low_upper, cls$high_lower)

bw <- 20
pa <- ggplot(tg, aes(dsci, fill = season)) +
  geom_histogram(binwidth = bw, boundary = 0, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = cuts, linetype = 2) +
  annotate("text", x = cuts, y = Inf, label = paste("DSCI", cuts), angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
  scale_fill_manual(values = season_cols) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  labs(x = "DSCI on the target day", y = "Targets", fill = NULL,
       title = sprintf("Target days (n = %d), class cut-offs dashed", nrow(tg))) +
  theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "inside",
                     legend.position.inside = c(0.85, 0.75))
pb <- ggplot(wk, aes(dsci)) +
  geom_histogram(binwidth = bw, boundary = 0, fill = "grey55", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = cuts, linetype = 2) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  labs(x = "DSCI", y = "Weeks", title = sprintf("All weekly Drought Monitor maps, 2014-2025 (n = %d)", nrow(wk))) +
  theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"))
ggsave(file.path(proj_dir, "plot_dsci_histogram.png"), pa / pb, width = 8, height = 7, dpi = 300, bg = "white")

q <- function(x) round(quantile(x, c(0, .1, .25, .33, .5, .67, .75, .9, 1)))
cat("\ntarget DSCI quantiles:\n"); print(q(tg$dsci))
cat("weekly DSCI quantiles:\n"); print(q(wk$dsci))
cat(sprintf("targets with DSCI = 0: %d; < 10: %d; >= 300: %d (of which 2022-2023: %d)\n",
            tg[dsci == 0, .N], tg[dsci < 10, .N], tg[dsci >= 300, .N], tg[dsci >= 300 & year %in% 2022:2023, .N]))
cat(sprintf("weeks with DSCI = 0: %.0f%%\n", 100 * wk[, mean(dsci == 0)]))
print(tg[, .N, keyby = .(band = cut(dsci, c(-1, 0, 50, 100, 200, 300, 500)))])

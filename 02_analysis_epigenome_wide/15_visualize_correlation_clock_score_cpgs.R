### Title: "Evaluating epigenetic clocks and scores: Cross-tissue correlations of clock/score CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2026-08-04"
### Purpose: Bisualize distribution of blood-saliva correlations for clock/score CpGs relative to their weights as raincloud plot

### Setup
## General
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
# libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

## Source functions and utilities
source("./00_functions/01_functions.R")
source("./00_functions/02_utilities.R")

## Load data
# epigenome-wide results
load("03_results/cross_tissue_correlations_clock_score_cpgs.Rdata")

## Prepare data
# extract correlation ranges
range(filter(crosstissue_clockscore_cpgs, adjustment == "unadjusted")$correlation)
range(filter(crosstissue_clockscore_cpgs, adjustment == "cell type adjusted")$correlation)
# add manual jitter for plotting 
crosstissue_clockscore_cpgs <- crosstissue_clockscore_cpgs %>%
  mutate(clock_num = as.numeric(clock),
         y_rain = clock_num - 0.42 + runif(n(), min = -0.06, max = 0.06))
# define levels
clock_levels <- levels(crosstissue_clockscore_cpgs$clock)

### Plot correlation distributions as raincloid plot
raincloudplot <- ggplot(crosstissue_clockscore_cpgs, 
                        aes(x = correlation, y = clock_num, group = clock_num)) +
  # distribution violin plot
  stat_slab(side = "top",
            scale = 0.45,
            fill = GyPr_palette[160],
            color = NA,
            alpha = 0.5,
            show.legend = FALSE) +
  # distribution boxplot
  geom_boxplot(width = 0.24,
               outlier.shape = NA,
               position = position_nudge(y = -0.16),
               color = "grey30",
               fill = "white",
               alpha = 0.85,
               linewidth = 0.4) +
  # dotplot colored by relative CpG weight
  geom_point(mapping = aes(x = correlation, y = y_rain, color = abs_weight_scaled),
             data = arrange(crosstissue_clockscore_cpgs, abs_weight_scaled),
             inherit.aes = FALSE,
             size = 1,
             alpha = 0.75) +
  scale_color_viridis_c(name = "Relative CpG weight", option = "D") +
  scale_y_continuous(breaks = seq_along(clock_levels),
                     labels = clock_levels,
                     expand = expansion(mult = 0, add = c(0.45, 0.55))) +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.25)) +
  facet_wrap(~adjustment, nrow = 1) +
  labs(x = "Spearman correlation (blood vs. saliva)", y = NULL) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.justification = "left")

raincloudplot
# export
ggsave(raincloudplot, filename = "supplementary_figure_7_raincloudplot_clockscore_correlation.png",
       path = "./04_figures", device = 'png', height = 10, width = 7, dpi = 700)

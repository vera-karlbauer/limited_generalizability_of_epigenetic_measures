### Title: "Limited generalizability epigenome-wide: Enrichment visualization for clock CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-11"
### Purpose: Visualize results from enrichment analyses for clock CpGs
### Purpose: Generate panel B of figure 8

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)

## Source functions and utilities
source("./00_functions/01_functions.R")
source("./00_functions/02_utilities.R")

## Load results
load("03_results/enrichment_clock_cpgs.Rdata")

### Prepare data
# # add 'all' category for bio and chrono clocks
# all_chrono <- as.data.frame(t(c("chronological", "all", NA_real_, 1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
# colnames(all_chrono) <- colnames(results_enrichment)
# all_bio <- as.data.frame(t(c("biological", "all", NA_real_, 1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
# colnames(all_bio) <- colnames(results_enrichment)
# results_enrichment <- rbind(results_enrichment, all_chrono, all_bio)

# add significance labels based on p_fdr
results_enrichment <- results_enrichment %>%
  mutate(ratio_overlap = as.numeric(ratio_overlap),
         or = as.numeric(or),
         confint_or_lower = as.numeric(confint_or_lower),
         confint_or_upper = as.numeric(confint_or_upper),
         p = as.numeric(p),
         p_fdr = as.numeric(p_fdr)) %>%
  mutate(significance = case_when(p_fdr >= 0.05 | is.na(p_fdr) ~ " ",
                                  p_fdr < 0.05 & p_fdr >= 0.01 ~ "*",
                                  p_fdr < 0.01 & p_fdr >= 0.001 ~ "**",
                                  p_fdr < 0.001 & p_fdr >= 0.0001 ~ "***",
                                  p_fdr < 0.0001 ~ "****",
                                  TRUE ~ NA))
# convert clock/score & enrichment type to factors for plotting
results_enrichment <- results_enrichment %>%
  mutate(clock_type = replace_values(clock_type, "chronological" ~ "chronological clocks",
                                     "biological" ~ "biological clocks",
                                     "episcores" ~ "epigenetic scores",
                                     "pc" ~ "PC clocks"),
         clock_type = factor(clock_type, levels = c("PC clocks", "epigenetic scores", 
                                                   "biological clocks", "chronological clocks")),
         enrichment_type = replace_values(enrichment_type, "variable_cpg" ~ "variable",
                                          "variable_and_correlated" ~ "varbl. & corr",
                                          "variable_and_correlated_celltype_adjusted" ~ "varbl. & corr celltype"),
         enrichment_type = factor(enrichment_type, levels = c("varbl. & corr celltype", "varbl. & corr", "variable")))

### Create plots
## forest plot with odds ratios
forestplot <- ggplot(data = results_enrichment,
                     aes(y = clock_type,
                         x = or, 
                         fill = enrichment_type, ,
                         color = enrichment_type)
                         ) +
  geom_errorbar(aes(xmin = confint_or_lower, xmax = confint_or_upper), width = 0.6, position = position_dodge(width = 1)) +
  geom_point(shape = 15, size = 3, alpha = 0.6, position=position_dodge(width = 1)) +
  geom_text(aes(label = significance, color = "black"), position = position_dodge(width = 1), vjust = -0.2, hjust = 0.5, size = 3) +
  geom_vline(xintercept = 1, linetype = "longdash") +
  geom_hline(yintercept = c(1.5, 2.5, 3.5), color = "grey", linetype = "longdash") +
  xlab("Enrichment OR (95% CI)") +
  # ylab("Clock/score type") +
  ggtitle("Enrichment of clock/score CpGs") +
  scale_color_discrete(type = rev(c(GyPr_palette[50], GyPr_palette[200], GyPr_palette[160])),
                      name = "Enrichment type", limits = levels(results_enrichment$enrichment_type)) +
  scale_fill_discrete(type = rev(c(GyPr_palette[50], GyPr_palette[200], GyPr_palette[160])),
                       name = "Enrichment type", limits = levels(results_enrichment$enrichment_type)) +
  scale_x_log10(limits=c(0.8,6)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  theme(legend.position = "inside", 
    legend.position.inside = c(0.82, 0.87),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y = element_text(margin = margin(t = 50, b = 50), angle = 40, hjust = 0.95, vjust = -0.5))
# axis.text.y = element_text(angle = 40, hjust=0.95),
forestplot
ggsave(filename = "figure_8b_forestplot_clock_enrichment.png", path = "./04_figures", device = 'png', height = 3.5, width = 5, dpi = 700)

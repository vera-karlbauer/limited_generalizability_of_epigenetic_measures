### Title: "Limited generalizability epigenome-wide: Enrichment visualization for clock CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-11"
### Purpose: Visualize results from enrichment analyses for clock CpGs
### Purpose: Generate panel B of figure 5

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
# add 'all' category for bio and chrono clocks
all_chrono <- as.data.frame(t(c("chronological", "all", NA_real_, 1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
colnames(all_chrono) <- colnames(results_clock_cpgs)
all_bio <- as.data.frame(t(c("biological", "all", NA_real_, 1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
colnames(all_bio) <- colnames(results_clock_cpgs)
results_clock_cpgs <- rbind(results_clock_cpgs, all_chrono, all_bio)

# add significance labels based on p_fdr
results_clock_cpgs <- results_clock_cpgs %>%
  mutate(ratio_overlap = as.numeric(ratio_overlap),
         or = as.numeric(or),
         confint_or_lower = as.numeric(confint_or_lower),
         confint_or_upper = as.numeric(confint_or_upper),
         p = as.numeric(p),
         p_fdr = as.numeric(p_fdr)) %>%
  mutate(significance = case_when(p_fdr > 0.05 | is.na(p_fdr) ~ "",
                                  p_fdr < 0.05 & p_fdr > 0.01 ~ "*",
                                  p_fdr < 0.01 & p_fdr > 0.001 ~ "**",
                                  p_fdr < 0.001 & p_fdr > 0.0001 ~ "***",
                                  p_fdr < 0.0001 ~ "****",
                                  TRUE ~ NA))

### Create plots
## forest plot with odds ratio
results_clock_cpgs <- results_clock_cpgs %>% filter(enrichment_type != "all")

forestplot <- ggplot(data = results_clock_cpgs,
                     aes(y = factor(clock_type, levels = rev(c("chronological", "biological"))),
                         x = or, 
                         fill = factor(enrichment_type, 
                                       levels = rev(c("variable_cpg", "variable_and_correlated", 
                                                      "variable_and_correlated_celltype_adjusted"))),
                         color = factor(enrichment_type, 
                                        levels = rev(c( "variable_cpg", "variable_and_correlated", 
                                                       "variable_and_correlated_celltype_adjusted")))
                         )) +
  geom_errorbar(aes(xmin = confint_or_lower, xmax = confint_or_upper), width = 0.6, position = position_dodge(width = 1)) +
  geom_point(shape = 15, size = 3, alpha = 0.6, position=position_dodge(width = 1)) +
  geom_text(aes(label = significance, color = "black"), position = position_dodge(width = 1), vjust = -0.7, hjust = 0.5) +
  geom_vline(xintercept = 1, linetype = "longdash") +
  geom_hline(yintercept = 1.5, color = "grey", linetype = "longdash") +
  xlab("Enrichment OR (95% CI)") +
  ylab("Clock type") +
  ggtitle("Enrichment of clock CpGs") +
  scale_color_discrete(type = rev(c(GyPr_palette[50], GyPr_palette[200], GyPr_palette[160])),
                      name = "Enrichment type",
                      limits = rev(c("variable_cpg", "variable_and_correlated",
                                     "variable_and_correlated_celltype_adjusted")),
                      labels = rev(c("variable", "varbl. & corr.", 
                                     "varbl. & corr. celltype")), drop = FALSE) +
  scale_fill_discrete(type = rev(c(GyPr_palette[50], GyPr_palette[200], GyPr_palette[160])),
                       name = "Enrichment type",
                       limits = rev(c("variable_cpg", "variable_and_correlated",
                                      "variable_and_correlated_celltype_adjusted")),
                       labels = rev(c("variable", "varbl. & corr.", 
                                      "varbl. & corr. celltype")), drop = FALSE) +
  scale_x_log10(limits=c(0.9,4.7)) +
  theme_bw() +
  theme(legend.position = "inside", 
    legend.position.inside = c(0.82, 0.78),
    axis.text.y = element_text(margin = margin(t = 50, b = 50), angle = 90, hjust = 0.5, vjust = 1.5))
forestplot
ggsave(filename = "figure_5b_forestplot_clock_enrichment.png", path = "./04_figures", device = 'png', height = 3.5, width = 5, dpi = 700)

### Title: "Limited generalizability epigenome-wide: Visualize results"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-02-20"
### Purpose: visualize & summarize results from epigenome-wide analyses
### Purpose: generate Panel A of figure 5

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DescTools)
library(RColorBrewer)
library(psych)
library(tidyr)

## Source functions and utlities
source("./00_functions/01_functions.R")
source("./00_functions/02_utilities.R")

## Import results
load("03_results/results_epigenome_wide.Rdata")

### Analysis - unadjusted
## All CpGs
# Mean within-person Spearman correlation
mean(results_epigenome$cor_spearman, na.rm = TRUE)
max(results_epigenome$cor_spearman, na.rm = TRUE)
min(results_epigenome$cor_spearman, na.rm = TRUE)

# N and % of significantly correlated CpGs (unadjusted)
PercTable(results_epigenome$significant_before_FDR)

# N and % of significantly correlated CpGs (after FDR correction)
PercTable(results_epigenome$significant_after_FDR)

## Variable CpGs only
results_epigenome_variable <- results_epigenome %>%
  filter(variable_cpg == TRUE)

# Mean within-person Spearman correlation
mean(results_epigenome_variable$cor_spearman, na.rm = TRUE)
max(results_epigenome_variable$cor_spearman)
min(results_epigenome_variable$cor_spearman)

# N and % of significantly correlated CpGs (before correction)
PercTable(results_epigenome_variable$significant_before_FDR)

# N and % of significantly correlated CpGs (after FDR correction)
PercTable(results_epigenome_variable$significant_after_FDR)

# N and % of highly correlated (rho>.5) and variable CpGs out of all CpGs
PercTable(results_epigenome$variable_and_correlated)

### Analysis - celltype-adjusted
## All CpGs
# Mean within-person Spearman correlation 
mean(results_epigenome$cor_spearman_celltype_adjusted, na.rm = TRUE)
max(results_epigenome$cor_spearman_celltype_adjusted)
min(results_epigenome$cor_spearman_celltype_adjusted)

# N and % of significantly correlated CpGs (unadjusted)
PercTable(results_epigenome$significant_before_FDR_celltype_adjusted)

# N and % of significantly correlated CpGs (after FDR correction)
PercTable(results_epigenome$significant_after_FDR_celltype_adjusted)

## Variable CpGs only
# Number and percentage of variable CpGs
PercTable(results_epigenome$variable_cpg)

# Mean within-person Spearman correlation 
mean(results_epigenome_variable$cor_spearman_celltype_adjusted, na.rm = TRUE)
max(results_epigenome_variable$cor_spearman_celltype_adjusted)
min(results_epigenome_variable$cor_spearman_celltype_adjusted)

# N and % of significantly correlated CpGs (unadjusted)
PercTable(results_epigenome_variable$significant_before_FDR_celltype_adjusted)

# N and % of significantly correlated CpGs (after FDR correction)
PercTable(results_epigenome_variable$significant_after_FDR_celltype_adjusted)

# N and % of highly correlated (rho>.5) and variable CpGs out of all CpGs
PercTable(results_epigenome$variable_and_correlated_celltype_adjusted)

## Distribution of within-person correlations for variable CpGs only
# convert to long format 
results_epigenome_nocelltype <- results_epigenome_variable %>%
  select(!(ends_with("_celltype_adjusted")))
results_epigenome_nocelltype$adjustment <- rep("unadjusted", nrow(results_epigenome_nocelltype))
results_epigenome_celltype <- results_epigenome_variable %>%
  select(-c("mean_blood", "mean_saliva", "sd_blood", "sd_saliva", "cor_spearman", 
            "p", "p_fdr", "significant_before_FDR", "significant_after_FDR", "icc", "icc_lower", "icc_upper",
            "variable_and_correlated"))
results_epigenome_celltype$adjustment <- rep("celltype_adjusted", nrow(results_epigenome_celltype))
colnames(results_epigenome_celltype) <- colnames(results_epigenome_nocelltype)
results_epigenome_long <- rbind(results_epigenome_nocelltype, results_epigenome_celltype)

# create histogram
hist <- ggplot(results_epigenome_long, aes(x = cor_spearman, fill = adjustment)) +
  geom_histogram(color = GyPr_palette[30],
                 breaks = seq(from = -0.6, to = 1, by = 0.05)) +
  scale_x_continuous(breaks = round(seq(from = -0.6, to = 1, by = 0.1), digits = 1),
                     minor_breaks = seq(from = -0.6, to = 1, by = 0.1)) +
  scale_y_continuous(breaks = round(seq(from = 0, to = 5000000, by = 50000), digits = 1)) +
  scale_fill_manual(values = GyPr_palette[c(120,160)]) +
  # add abline at corr >0.5
  geom_vline(xintercept = 0.5, linetype = 2, color = brewer.pal(3, "Set1")[1]) + 
  # add abline label
  annotate("label", x=0.7, y=40000, size = 4,
           label = "rho > 0.5",
           fill = brewer.pal(3, "Set1")[1]) +
  # add abline at mean unadjusted
  geom_vline(xintercept = mean(filter(results_epigenome_long, adjustment == "unadjusted")$cor_spearman, 
                               na.rm = TRUE),
             linetype = 2, 
             color = GyPr_palette[200]) +
  # add abline at mean adjusted
  geom_vline(xintercept = mean(filter(results_epigenome_long, adjustment == "celltype_adjusted")$cor_spearman, 
                               na.rm = TRUE),
             linetype = 2,
             color = GyPr_palette[160]) +
  # add abline labels
  annotate("label", x=0.3, y=80000, size = 4,
           label=paste0("Mean correlation (unadjusted) = ",
                        round(mean(filter(results_epigenome_long, adjustment == "unadjusted")$cor_spearman, na.rm = TRUE), 
                              digits =2)),
           fill = GyPr_palette[160]) +
  annotate("label", x=0.4, y=110000, size = 4,
           label=paste0("Mean correlation (celltype-adjusted) = ",
                        round(mean(filter(results_epigenome_long, adjustment == "celltype_adjusted")$cor_spearman, na.rm = TRUE), 
                              digits =2)),
           fill = GyPr_palette[120]) +
  xlab("Blood-saliva correlation (Spearman's rho)") +
  ylab("Observed frequency") +
  ggtitle("Cross-tissue correlations for variable CpGs") +
  theme_bw() +
  theme(legend.position = "none", axis.text.y = element_text(angle = 90, hjust = 0.5))
hist
ggsave(hist, filename = "04_figures/figure_5a_histogram_epigenome_wide_within_variable.png", device = "png",
       width = 5, height = 3.5, units = "in", dpi = 700) 

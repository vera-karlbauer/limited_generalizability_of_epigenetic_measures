### Title: "Limited generalizability epigenome-wide: mean cross-tissue correlation of EWAS hits"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-11"
### Purpose: For different categories from EWAS catalog, correlate mean cross-tissue correlation (CpGs weighted by betas)
### Purpose: Generate panel C of figure 5

### Setup
## general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(data.table)
library(writexl)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(Hmisc)
library(circhelp)

## Source functions and utilities
source("./00_functions/01_functions.R")
source("./00_functions/02_utilities.R")

## Load data
# epigenome-wide results
load("03_results/results_epigenome_wide.Rdata")
# EWAS catalog (https://www.ewascatalog.org/download/)
ewascatalog <- read.table("../../../../tools/ewascatalog/ewascatalog-results.txt", header = TRUE, sep = "\t")
# EWAS catalog trait annotation
ewas_studies <- fread("../../../../tools/ewascatalog/ewascatalog_studies_annotated.txt", sep = "\t")

## Prepare data
# subset epigenome-wide blood-saliva correlation results to variable cpgs
results_epigenome_variable <- filter(results_epigenome, variable_cpg == TRUE)
# subset EWAS catalog for CpGs that were tested and are variable in epigenome-wide blood-saliva correlation analyses
ewascatalog <- ewascatalog %>%
  filter(CpG %in% results_epigenome_variable$CpG_name)
# link CpG associations to trait supercategories
ewas_studies <- ewas_studies %>%
  select(StudyID, trait_category)
ewascatalog <- left_join(ewascatalog, ewas_studies, by = "StudyID")
# show number of unique trait categories
length(unique(ewascatalog$trait_category, na.rm = TRUE))

# remove all CpGs without reported betas
ewascatalog <- ewascatalog %>%
  filter(!is.na(Beta))
# only keep trait categories with >100 traits
ewascatalog <- ewascatalog %>%
  group_by(trait_category) %>%
  filter(n() > 100) %>%
  ungroup()
# show number of unique trait categories after filtering
length(unique(ewascatalog$trait_category,  na.rm = TRUE))
# list unique traits
unique(ewascatalog$trait_category, na.rm = TRUE)
# show number of unique studies after filtering
length(unique(ewascatalog$StudyID, na.rm = TRUE))

### Calculate weighted mean cross-tissue correlation of each EWAS trait with and without cell type adjustment
## prepare loop
# distinct traits
traits <- na.exclude(unique(ewascatalog$trait_category))
results_ewas_correlation <- as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(results_ewas_correlation) <- c("trait", "number_cpgs", "mean_cor_spearman", "se_cor_spearman",
                                   "mean_cor_spearman_celltype_adjusted", "se_cor_spearman_celltype_adjusted")
# run loop over all traits, calculate both with and without cell type adjustment
for(i in 1:length(traits)){
  trait = traits[i]
  print(paste0("Trait: ", trait))
  # set up empty results object
  current_results <- as.data.frame(matrix(nrow = 1, ncol = 0))
  current_results$trait = trait 
  # get weighted betas for each CpG corresponding to trait 
  weights <- ewascatalog %>%
    filter(trait_category == trait) %>%
    select(CpG, Beta) %>%
    group_by(CpG) %>%
    summarise(mean_abs_beta = mean(abs(Beta), na.rm = TRUE), .groups = "drop")
  current_results$number_cpgs = nrow(weights)
  ## without cell typ adjustment
  # get correlations
  cor_nocelltype <- results_epigenome %>%
    filter(CpG_name %in% weights$CpG) %>%
    select(CpG_name, cor_spearman) %>%
    arrange(match(CpG_name, weights$CpG))
  # calculate weighted mean
  current_results$mean_cor_spearman = weighted.mean(cor_nocelltype$cor_spearman, weights$mean_abs_beta)
  # calculated weighted standard error
  current_results$se_cor_spearman = weighted_sem(cor_nocelltype$cor_spearman, weights$mean_abs_beta)
  ## with cell type adjustment
  # get correlations
  cor_celltype <- results_epigenome %>%
    filter(CpG_name %in% weights$CpG) %>%
    select(CpG_name, cor_spearman_celltype_adjusted) %>%
    arrange(match(CpG_name, weights$CpG))
  # calculate weighted mean
  current_results$mean_cor_spearman_celltype_adjusted = weighted.mean(cor_celltype$cor_spearman_celltype_adjusted, weights$mean_abs_beta)
  # calculated weighted standard error 
  current_results$se_cor_spearman_celltype_adjusted = weighted_sem(cor_celltype$cor_spearman, weights$mean_abs_beta)
  ## append results
  results_ewas_correlation <- rbind(results_ewas_correlation, current_results)
}


test <- results_ewas_correlation

### Visualize results
# convert data to long format
results_ewas_correlation_long <- results_ewas_correlation %>%
  pivot_longer(cols = c(mean_cor_spearman, mean_cor_spearman_celltype_adjusted, se_cor_spearman, se_cor_spearman_celltype_adjusted),
             names_to = c(".value", "adjustment"),
             names_pattern = "(.*?)(_celltype_adjusted)?$",
             values_drop_na = FALSE) %>%
  mutate(adjustment = case_when(adjustment == "_celltype_adjusted" ~ "celltype-adjusted",
                                TRUE ~ "unadjusted"))  %>%
  group_by(trait) %>%
  arrange(desc(mean_cor_spearman)) %>%
  ungroup() %>%
  mutate(se_lower = mean_cor_spearman - se_cor_spearman,
         se_upper = mean_cor_spearman + se_cor_spearman)
# create barplot
ewasplot <- ggplot(data = results_ewas_correlation_long,
                     aes(x = mean_cor_spearman,
                         y = reorder(trait, mean_cor_spearman),
                         fill = adjustment,
                         color = adjustment)) +
  geom_errorbar(aes(xmin = se_lower, xmax = se_upper), width = 0.6, position = position_dodge(width = 0.9)) +
  # geom_col(color = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.6) +
  # geom_point(shape = 15, size = 2, alpha = 0.6, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin = 0, xmax = mean_cor_spearman), size = 3, position = position_dodge(width = 0.9), alpha = 0.6) +
  ylab("Traits from EWAS catalog") +
  xlab("Weighted mean cross-tissue correlation of trait CpGs") +
  xlim(0,0.4) +
  ggtitle("Cross-tissue correlations of EWAS traits") +
  scale_fill_discrete(type = c(GyPr_palette[160], GyPr_palette[200]),
                       name = "Adjustment",
                       breaks = c("celltype-adjusted", "unadjusted")) +
  scale_color_discrete(type = c(GyPr_palette[160], GyPr_palette[200]),
                      name = "Adjustment",
                      breaks = c("celltype-adjusted", "unadjusted")) +
  # add hlines to separate trait categories
  geom_hline(yintercept = seq(1.5, 10, 1), color = "grey", linewidth = 0.25) +
  theme_bw() + 
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.15),
        axis.text.y = element_text(angle = 40, hjust=0.95),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill="white", color = "white"))
ewasplot
ggsave(filename = "figure_5c_forestplot_ewas_correlation.png", path = "./04_figures", device = 'png', height = 7, width = 5, dpi = 700)

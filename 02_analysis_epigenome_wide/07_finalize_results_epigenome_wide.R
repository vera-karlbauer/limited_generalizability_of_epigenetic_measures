### Title: "Limited generalizability epigenome-wide: Combine uncorrected and cell type corrected results
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-02-20"
### Purpose: combine epigenome-wide correlation and ICC results and export (eTable 11)

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(writexl)

### Load data
# no cell type correction
load("03_results/results_nocelltype_combined_annotated.Rdata")

# cell type correction
load("03_results/results_celltype_combined_annotated.Rdata")

# variable CpGs across tissues
load("../../../01_data/02_proc/09_cross_tissue/variable_cpgs_crosstissue.Rdata")

### Join results
# remove annotation and residualized means/sds from cell type corrected df
results_celltype_combined <- results_celltype_combined %>%
  select(probe_ID, cor_spearman, p, p_fdr, icc, icc_lower, icc_upper)

# rename columns with cell type correction
colnames(results_celltype_combined) <- paste0(colnames(results_celltype_combined), "_celltype_adjusted")
results_celltype_combined <- results_celltype_combined %>%
  rename(probe_ID = probe_ID_celltype_adjusted)

# join by probe_ID
results_epigenome <- full_join(results_nocelltype_combined, results_celltype_combined,
                               join_by("probe_ID"))
results_epigenome <- results_epigenome %>%
  relocate(ends_with("_celltype_adjusted"), .before = genesUniq)

# check whether all CpGs have complete data
results_epigenome %>%
  filter(is.na(cor_spearman) | is.na(cor_spearman_celltype_adjusted)) %>%
  nrow()

# remove geneNames column (unnecessary)
results_epigenome <- results_epigenome %>%
  select(-geneNames)

# add significance yes/no columns
# before FDR, no celltype
results_epigenome <- results_epigenome %>%
  mutate(significant_before_FDR = case_when(p < 0.05 ~ TRUE, 
                                            p >= 0.05 ~ FALSE,
                                            TRUE ~ NA),
         # after FDR, no celltype
         significant_after_FDR = case_when(p_fdr < 0.05 ~ TRUE, 
                                           p_fdr >= 0.05 ~ FALSE,
                                           TRUE ~ NA),
         # before FDR, celltype
        significant_before_FDR_celltype_adjusted = case_when(p_celltype_adjusted < 0.05 ~ TRUE, 
                                                                     p_celltype_adjusted >= 0.05 ~ FALSE,
                                                                     TRUE ~ NA),
         # after FDR, celltype
        significant_after_FDR_celltype_adjusted = case_when(p_fdr_celltype_adjusted < 0.05 ~ TRUE, 
                                                                   p_fdr_celltype_adjusted >= 0.05 ~ FALSE,
                                                                   TRUE ~ NA)
         )

# add variable cpg yes/no column
results_epigenome <- results_epigenome %>%
  mutate(variable_cpg = case_when(
    probe_ID %in% variable_cpgs_crosstissue ~ TRUE, TRUE ~ FALSE))

## indicate relevant correlation with threshold of >0.5 amongst variable CpGs
results_epigenome <- results_epigenome %>%
  mutate(variable_and_correlated = case_when(variable_cpg == TRUE & cor_spearman > 0.5 ~ TRUE,
                                           TRUE ~ FALSE),
       variable_and_correlated_celltype_adjusted = case_when(variable_cpg == TRUE & cor_spearman_celltype_adjusted > 0.5 ~ TRUE,
                                                             TRUE ~ FALSE))

## round correlations and ICC estimates to 4 digits (leave p-values untouched)
results_epigenome <- results_epigenome %>%
  mutate(across(c(mean_blood, sd_blood, mean_saliva, sd_saliva, cor_spearman, icc, icc_lower, icc_upper,
                  cor_spearman_celltype_adjusted, icc_celltype_adjusted, 
                  icc_lower_celltype_adjusted, icc_upper_celltype_adjusted), 
                ~ round(.x, digits = 4)))

### Export
save(results_epigenome, file = "03_results/results_epigenome_wide.Rdata")
write.csv(results_epigenome, file = "05_tables_for_publication/etable_11_blood_saliva_epigenome_wide.csv", row.names = FALSE)

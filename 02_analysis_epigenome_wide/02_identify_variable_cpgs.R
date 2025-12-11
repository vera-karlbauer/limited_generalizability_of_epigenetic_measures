### Title: "Limited generalizability epigenome-wide: Select variable CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-04-01"
### Purpose: Select variable CpGs for corrected and uncorrected data
### Note: to be run on HPC cluster

# Variable CpGs identified according to Hannon et al. (2015): calculate
# DNA methylation difference between 10th and 90th percentile across 
# all samples, then select sites where this is > 5%

## Step 0: load packages
library(dplyr)
library(stringr)

## Step 1: load data
# EPICv2 filtered data: blood
load("../../../01_data/02_proc/01_epic/kids2health_sp1_2_4_blood_new/processed_data/final_data/Betas_clean_filtered_quantile_bmiq_combated_final.Rdata")
betas_filtered_blood <- Betas_clean_filtered_quantile_bmiq_combated
remove(Betas_clean_filtered_quantile_bmiq_combated)

# EPICv2 filtered data: saliva
load("../../../01_data/02_proc/01_epic/kids2health_sp1_2_4_saliva/processed_data/final_data/Betas_clean_filtered_quantile_bmiq_combated_final.Rdata")
betas_filtered_saliva <- Betas_clean_filtered_quantile_bmiq_combated
remove(Betas_clean_filtered_quantile_bmiq_combated)

# relevant phenotypes: 
load("../../../01_data/02_proc/09_cross_tissue/all_cross_tissue_inner.Rdata")

# common CpGs in both tissues
load("../../../01_data/02_proc/09_cross_tissue/common_cpgs.Rdata")

## Step 2: subset to common CpGs and individuals
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  filter(timepoint == "T0")
# blood
# convert all dashes to colons in cpg names (affects 1 cpg)
rownames(betas_filtered_blood) <- str_replace_all(rownames(betas_filtered_blood) , "\\-", "\\.")
betas_filtered_blood <- betas_filtered_blood[,all_cross_tissue_inner$arrayid_blood]
betas_filtered_blood <- betas_filtered_blood[common_cpgs,]

# saliva
# convert all dashes to colons in cpg names (affects 1 cpg)
rownames(betas_filtered_saliva) <- str_replace_all(rownames(betas_filtered_saliva) , "\\-", "\\.")
betas_filtered_saliva <- betas_filtered_saliva[,all_cross_tissue_inner$arrayid_saliva]
betas_filtered_saliva <- betas_filtered_saliva[common_cpgs,]

## Step 3: define function for calculating sample differences
calculate_percentile_difference <- function(x) {
  quantiles <- quantile(x, probs = c(0.1, 0.9), na.rm = TRUE)
  return(quantiles[2] - quantiles[1])
}

## Step 4: select variable CpGs in blood 
mdiff_blood <- apply(betas_filtered_blood, 1, calculate_percentile_difference)
variable_cpgs_blood <- betas_filtered_blood[mdiff_blood > 0.05, ]
print(paste0(nrow(variable_cpgs_blood), " variable CpGs identified in blood."))

## Step 5: select variable CpGs in saliva
mdiff_saliva <- apply(betas_filtered_saliva, 1, calculate_percentile_difference)
variable_cpgs_saliva <- betas_filtered_saliva[mdiff_saliva > 0.05, ]
print(paste0(nrow(variable_cpgs_saliva), " variable CpGs identified in saliva"))

## Step 6: subset to variable CpGs in both tissues
variable_cpgs_crosstissue <- base::intersect(rownames(variable_cpgs_blood), rownames(variable_cpgs_saliva))
print(paste0(length(variable_cpgs_crosstissue), " variable CpGs identified in both tissues."))

## Step 7: export list of variable CpGs
save(variable_cpgs_crosstissue, file = "../../../01_data/02_proc/09_cross_tissue/variable_cpgs_crosstissue.Rdata")

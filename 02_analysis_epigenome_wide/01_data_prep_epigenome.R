### Title: "Limited generalizability epigenome-wide: data preparation"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-01-28"
### Purpose: create datasets used for epigenome-wide cross-tissue correlation analyses
### Note: to be run on HPC cluster

## Step 0: Load packages

# load packages
library(rmarkdown)
library(dplyr)
library(tibble)
library(janitor)
library(stringr)

## Step 1: Load data
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

## Step 2: Extract common individuals in both datasets (T0 only)
# filter common individuals at T0
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  filter(timepoint == "T0")
print(paste0("Number of individuals at baseline with both tissues: ", nrow(all_cross_tissue_inner)))

# subset blood
betas_filtered_blood <- betas_filtered_blood[,all_cross_tissue_inner$arrayid_blood]
rownames(betas_filtered_blood) <- str_replace_all(rownames(betas_filtered_blood), "\\-", "\\.")

# subset saliva
betas_filtered_saliva <- betas_filtered_saliva[,all_cross_tissue_inner$arrayid_saliva]
rownames(betas_filtered_saliva) <- str_replace_all(rownames(betas_filtered_saliva), "\\-", "\\.")

## Step 3: Extract common CpGs in both datasets
# generate list of common CpGs
common_cpgs <- intersect(rownames(betas_filtered_blood), rownames(betas_filtered_saliva))
print(paste0("Number of common CpGs post-QC in both tissues: ", length(common_cpgs)))
# convert all dashes to colons in cpg_names (affects 1 cpg)
common_cpgs <- str_replace_all(common_cpgs, "\\-", "\\.")
# export common cpgs
save(common_cpgs, file = "../../../01_data/02_proc/09_cross_tissue/common_cpgs.Rdata")

# subset blood
betas_filtered_blood <- betas_filtered_blood[common_cpgs,]

# subset saliva
betas_filtered_saliva <- betas_filtered_saliva[common_cpgs,]

## Step 4: Reshape data
# rows = participants, 2 columns per CpG, format = "CpGname_blood", "CpGname_saliva"
# transpose & rename blood & append child_id (= cross-tissue subject identifier)
betas_final_blood <- data.frame(t(betas_filtered_blood))
remove(betas_filtered_blood)
colnames(betas_final_blood) <- paste0(colnames(betas_final_blood), "_blood")
betas_final_blood <- betas_final_blood %>%
  rownames_to_column(var = "arrayid_blood")
betas_final_blood <- right_join(select(all_cross_tissue_inner, arrayid_blood, child_id), betas_final_blood, by = "arrayid_blood")
# transpose & rename saliva & append child_id (= cross-tissue subject identifier)
betas_final_saliva <- data.frame(t(betas_filtered_saliva))
remove(betas_filtered_saliva)
colnames(betas_final_saliva) <- paste0(colnames(betas_final_saliva), "_saliva")
betas_final_saliva <- betas_final_saliva %>%
  rownames_to_column(var = "arrayid_saliva")
betas_final_saliva <- right_join(select(all_cross_tissue_inner, arrayid_saliva, child_id), betas_final_saliva, by = "arrayid_saliva")
# join datasets
betas_filtered_crosstissue <- inner_join(betas_final_blood, betas_final_saliva, by = "child_id")
# remove datasets from environment when no longer needed
remove(betas_final_blood)
remove(betas_final_saliva)

## Step 5: Extract cell types (needed for residualization of each CpG)
# Blood: residualize for everything except neutrophils.
# Saliva: residualize for epithelial cells
celltypes_crosstissue <- all_cross_tissue_inner %>%
  filter(arrayid_blood %in% betas_filtered_crosstissue$arrayid_blood) %>%
  select(child_id, arrayid_blood, arrayid_saliva,
         CD4Tnv_blood, Baso_blood, CD4Tmem_blood, Bmem_blood, Bnv_blood, Treg_blood, 
         CD8Tmem_blood, CD8Tnv_blood, Eos_blood, NK_blood, Mono_blood,
         Epithelial_saliva) 

## Step 6: Separate methylation data into 70 chunks (CpG-wise) for parallel analysis & export
# define chunks
chunks <- split(common_cpgs, rep_len(1:70, length(common_cpgs)))

# partition and export datasets
for(i in 1:length(chunks)) {
  colnames = c(paste0(unlist(chunks[i]), "_blood"), paste0(unlist(chunks[i]), "_saliva"))
  betas = select(betas_filtered_crosstissue, child_id, arrayid_blood,
                arrayid_saliva, any_of(colnames)
  )
  # save chunk
  filename = paste0("../../../01_data/02_proc/09_cross_tissue/70chunks/betas_filtered_crosstissue_", i, ".Rdata")
  save(betas, file = filename)
  # remove object from environment
  remove(betas)
  # save CpG names
  cpg_names = unlist(chunks[i])
  filename = paste0("../../../01_data/02_proc/09_cross_tissue/70chunks/common_cpgs_", i, ".Rdata")
  save(cpg_names, file = filename)
}

## Step 7: Export additional files
# # cell types
save(celltypes_crosstissue, file = "../../../01_data/02_proc/09_cross_tissue/celltypes_crosstissue.Rdata")

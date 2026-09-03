### Title: "Evaluating epigenetic clocks and scores: Score refinement"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2026-08-13"
### Purpose: Explore characteristics of epigenetic scores when filtered for tissue-stable CpGs

### Setup
## General
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
# libraries
library(dplyr)
library(writexl)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(meffonym, verbose = FALSE)
library(readxl)
library(janitor)
library(tibble)

## Source functions and utilities
source("./00_functions/01_functions.R")
source("./00_functions/02_utilities.R")

## Load data
# clock/score CpGs with weights and cross-tissue correlations
load("03_results/cross_tissue_correlations_clock_score_cpgs.Rdata")
# original score/clock data
load("../../../01_data/02_proc/09_cross_tissue/all_cross_tissue_inner.Rdata")
# EPICv2 data: blood
load("../../../01_data/02_proc/01_epic/kids2health_sp1_2_4_blood_new/processed_data/final_data/Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1_final.Rdata")
betas_blood <- Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1
remove(Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1)
# EPICv2 data: saliva
load("../../../01_data/02_proc/01_epic/kids2health_sp1_2_4_saliva/processed_data/final_data/Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1_final.Rdata")
betas_saliva <- Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1
remove(Betas_clean_unfiltered_quantile_bmiq_combated_pseudo_v1)

### Process data
episcore_cpgs <- crosstissue_clockscore_cpgs %>%
  filter(clock %in% c("CRP (Ligthart, 2016)", "Epigenetic-g") & adjustment == "unadjusted") 

### Refine scores: Epigenetic CRP (Lighthart et al.)
## get number of available CpGs in score when filtering for ICC > 0.10
n_refined <- episcore_cpgs %>%
  filter(clock == "CRP (Ligthart, 2016)") %>%
  filter(icc > 0.10) %>%
  nrow()
n_original <- episcore_cpgs %>%
  filter(clock == "CRP (Ligthart, 2016)") %>%
  nrow()
print(paste0(n_refined, " out of ", n_original, " CpGs remaining (", round(n_refined/n_original*100, digits = 2), "%)"))

## refine score - blood
weights <- episcore_cpgs %>%
  filter(clock == "CRP (Ligthart, 2016)") %>%
  filter(icc > 0.10) %>%
  filter(cpg %in% rownames(betas_blood))
weighted_betas_blood <- as.data.frame(betas_blood) %>%
  rownames_to_column(var = "CpG_name") %>%
  filter(CpG_name %in% weights$cpg) %>%
  arrange(match(CpG_name, weights$cpg)) %>%
  mutate(across(where(is.numeric), ~ .x * weights$weight)) %>%
  remove_rownames() %>%
  column_to_rownames(var = "CpG_name")
# sanity check: cpg order identical?
weights$cpg == rownames(weighted_betas_blood)
sumscore_blood <- weighted_betas_blood %>%
  dplyr::summarise_all(.funs = sum) %>%
  t() %>%
  as.data.frame() %>%
  mutate(V1 = as.numeric(V1)) %>%
  rename("refined_epigenetic_crp_score_ligthart_2016_blood" = "V1") %>%
  rownames_to_column(var = "arrayid_blood")
## append and create celltype-residualized version
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  left_join(sumscore_blood, by = "arrayid_blood") %>%
  mutate(refined_epigenetic_crp_score_ligthart_2016_resid_blood = 
           lm(refined_epigenetic_crp_score_ligthart_2016_blood 
              ~ CD4Tnv_blood + Baso_blood + CD4Tmem_blood + Bmem_blood + Bnv_blood + Treg_blood 
              + CD8Tmem_blood + CD8Tnv_blood + Eos_blood + NK_blood + Mono_blood)$residuals)

## refine score - saliva
weights <- episcore_cpgs %>%
  filter(clock == "CRP (Ligthart, 2016)") %>%
  filter(icc > 0.10) %>%
  filter(cpg %in% rownames(betas_saliva))
weighted_betas_saliva <- as.data.frame(betas_saliva) %>%
  rownames_to_column(var = "CpG_name") %>%
  filter(CpG_name %in% weights$cpg) %>%
  arrange(match(CpG_name, weights$cpg)) %>%
  mutate(across(where(is.numeric), ~ .x * weights$weight)) %>%
  remove_rownames() %>%
  column_to_rownames(var = "CpG_name")
# sanity check: cpg order identical?
weights$cpg == rownames(weighted_betas_saliva)
sumscore_saliva <- weighted_betas_saliva %>%
  dplyr::summarise_all(.funs = sum) %>%
  t() %>%
  as.data.frame() %>%
  mutate(V1 = as.numeric(V1)) %>%
  rename("refined_epigenetic_crp_score_ligthart_2016_saliva" = "V1") %>%
  rownames_to_column(var = "arrayid_saliva")
## append and create celltype-residualized version
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  left_join(sumscore_saliva, by = "arrayid_saliva") %>%
  mutate(refined_epigenetic_crp_score_ligthart_2016_resid_saliva = 
           lm(refined_epigenetic_crp_score_ligthart_2016_saliva 
              ~ Epithelial_saliva)$residuals)

### Refine scores: epigenetic-g
## get number of available CpGs in score when filtering for ICC > 0.10
n_refined <- episcore_cpgs %>%
  filter(clock == "Epigenetic-g") %>%
  filter(icc > 0.10) %>%
  nrow()
n_original <- episcore_cpgs %>%
  filter(clock == "Epigenetic-g") %>%
  nrow()
print(paste0(n_refined, " out of ", n_original, " CpGs remaining (", round(n_refined/n_original*100, digits = 2), "%)"))
## refine score - blood
weights <- episcore_cpgs %>%
  filter(clock == "Epigenetic-g") %>%
  filter(icc > 0.10) %>%
  filter(cpg %in% rownames(betas_blood))
weighted_betas_blood <- as.data.frame(betas_blood) %>%
  rownames_to_column(var = "CpG_name") %>%
  filter(CpG_name %in% weights$cpg) %>%
  arrange(match(CpG_name, weights$cpg)) %>%
  mutate(across(where(is.numeric), ~ .x * weights$weight)) %>%
  remove_rownames() %>%
  column_to_rownames(var = "CpG_name")
# sanity check: cpg order identical?
weights$cpg == rownames(weighted_betas_blood)
sumscore_blood <- weighted_betas_blood %>%
  dplyr::summarise_all(.funs = sum) %>%
  t() %>%
  as.data.frame() %>%
  mutate(V1 = as.numeric(V1)) %>%
  rename("refined_epigenetic_g_blood" = "V1") %>%
  rownames_to_column(var = "arrayid_blood")
## append and create celltype-residualized version
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  left_join(sumscore_blood, by = "arrayid_blood") %>%
  mutate(refined_epigenetic_g_resid_blood = 
           lm(refined_epigenetic_g_blood
              ~ CD4Tnv_blood + Baso_blood + CD4Tmem_blood + Bmem_blood + Bnv_blood + Treg_blood 
              + CD8Tmem_blood + CD8Tnv_blood + Eos_blood + NK_blood + Mono_blood)$residuals)

## refine score - saliva
weights <- episcore_cpgs %>%
  filter(clock == "Epigenetic-g") %>%
  filter(icc > 0.10) %>%
  filter(cpg %in% rownames(betas_saliva))
weighted_betas_saliva <- as.data.frame(betas_saliva) %>%
  rownames_to_column(var = "CpG_name") %>%
  filter(CpG_name %in% weights$cpg) %>%
  arrange(match(CpG_name, weights$cpg)) %>%
  mutate(across(where(is.numeric), ~ .x * weights$weight)) %>%
  remove_rownames() %>%
  column_to_rownames(var = "CpG_name")
# sanity check: cpg order identical?
weights$cpg == rownames(weighted_betas_saliva)
sumscore_saliva <- weighted_betas_saliva %>%
  dplyr::summarise_all(.funs = sum) %>%
  t() %>%
  as.data.frame() %>%
  mutate(V1 = as.numeric(V1)) %>%
  rename("refined_epigenetic_g_saliva" = "V1") %>%
  rownames_to_column(var = "arrayid_saliva")
## append and create celltype-residualized version
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  left_join(sumscore_saliva, by = "arrayid_saliva") %>%
  mutate(refined_epigenetic_g_resid_saliva = 
           lm(refined_epigenetic_g_saliva 
              ~ Epithelial_saliva)$residuals)

### Examine refined score characteristics: Epigenetic CRP Lighthart
## Correlation with original score
# blood
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_blood, 
         all_cross_tissue_inner$epigenetic_crp_score_ligthart_2016_blood)
# saliva
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_saliva, 
         all_cross_tissue_inner$epigenetic_crp_score_ligthart_2016_saliva)
## Cross-tissue correlation
all_cross_tissue_inner <- all_cross_tissue_inner %>%
  filter(timepoint == "T0")
# without cell type adjustment
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_blood, 
         all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_saliva,
         method = "pearson")
# cell type adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_resid_blood, 
         all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_resid_saliva,
         method = "pearson")

## Performance
# blood performance
# unadjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_blood, 
         all_cross_tissue_inner$crp_pgML_ln)
# celltype adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_resid_blood, 
         all_cross_tissue_inner$crp_pgML_ln)
# saliva performance
# unadjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_saliva, 
         all_cross_tissue_inner$crp_pgML_ln)
# celltype adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_crp_score_ligthart_2016_resid_saliva, 
         all_cross_tissue_inner$crp_pgML_ln)

### Examine refined score characteristics: Epigenetic-g
## Correlation with original score
# blood
cor.test(all_cross_tissue_inner$refined_epigenetic_g_blood, 
         all_cross_tissue_inner$epigenetic_g_blood)
# saliva
cor.test(all_cross_tissue_inner$refined_epigenetic_g_saliva, 
         all_cross_tissue_inner$epigenetic_g_saliva)
## Cross-tissue correlation
# without cell type adjustment
cor.test(all_cross_tissue_inner$refined_epigenetic_g_blood, 
         all_cross_tissue_inner$refined_epigenetic_g_saliva,
         method = "pearson")
# cell type adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_g_resid_blood, 
         all_cross_tissue_inner$refined_epigenetic_g_resid_saliva,
         method = "pearson")

## Performance
# blood performance
# unadjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_g_blood, 
         all_cross_tissue_inner$iq_sonr_nonverbal)
# celltype adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_g_resid_blood, 
         all_cross_tissue_inner$iq_sonr_nonverbal)
# saliva performance
# unadjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_g_saliva, 
         all_cross_tissue_inner$iq_sonr_nonverbal)
# celltype adjusted
cor.test(all_cross_tissue_inner$refined_epigenetic_g_resid_saliva, 
         all_cross_tissue_inner$iq_sonr_nonverbal)

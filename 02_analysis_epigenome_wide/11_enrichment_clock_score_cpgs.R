### Title: "Limited generalizability epigenome-wide: Enrichment of clock CpGs for cross-tissue CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-12"
### Purpose: Test CpGs on chronological clocks, biological clocks, epigenetic scores, and PC clocks/scores for enrichment for subset of variable and tissue-conserved CpGs

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(methylclockData)
library(DunedinPACE)

## Source functions
source("./00_functions/01_functions.R")

## Load data
# clock/score CpGs with weights and cross-tissue correlations
load("03_results/cross_tissue_correlations_clock_score_cpgs.Rdata")
load("03_results/results_epigenome_wide.Rdata")

## Prepare data
cpgs <- list()
# extract chronological cpgs
cpgs[["chronological"]] <- unique(filter(crosstissue_clockscore_cpgs, clock %in% c("Horvath", "SkinBlood", "Hannum", "Wu", "PedBE"))$cpg)
# extract biological cpgs
cpgs[["biological"]] <- unique(filter(crosstissue_clockscore_cpgs, clock %in% c("PhenoAge", "GrimAge", "GrimAge2", "DunedinPACE"))$cpg)
# extract score cpgs
cpgs[["episcores"]] <- unique(filter(crosstissue_clockscore_cpgs, clock %in% c("BMI (Do, 2023)", "BMI (Wahl, 2017)", "CRP (Wielscher, 2022)", "CRP (Ligthart, 2016)", "Epigenetic-g", "DNAmTL"))$cpg)
# extract PC clock cpgs
cpgs[["pc"]] <- unique(filter(crosstissue_clockscore_cpgs, clock %in% c("PCHorvath", "PCSkinBlood", "PCHannum", "PCPhenoAge", "PCGrimAge", "PCDNAmTL"))$cpg)

### Test for enrichment
# test enrichment of a) chronological clock cpgs, b) biological clock cpgs, c) epigenetic score cpgs
# 1) variable cpgs, 2) variable and highly tissue-correlated cpgs, 3) variable and highly tissue-correlated cpgs with celltype adjustment
# set up empty results object
results_enrichment <- as.data.frame(matrix(nrow = 0, ncol = 9))
colnames(results_enrichment) <- c("clock_type", "enrichment_type", "size_overlap", "ratio_overlap", "or", "or_confint_lower", "or_confint_upper", "p")
types <- c("variable_cpg", "variable_and_correlated", "variable_and_correlated_celltype_adjusted")
cgset <- c("chronological", "biological", "episcores", "pc")

# loop over chrono, bio, score, pc CpGs
for(i in 1: length(cgset)) {
  current_cgset = cpgs[[cgset[i]]]
  print(paste0("Testing: ", cgset[i], " clock/score CpGs"))
  # loop over variable, variable & correlated, variable and correlated /w celltype adjustement
  for(k in 1: length(types)){
    type = types[k]
    print(paste0("Testing enrichment for: ", type))
    # select set of CpGs
    # set up empty results object
    current_results <- as.data.frame(matrix(nrow = 1, ncol = 0))
    current_results$clock_type = cgset[i]
    current_results$enrichment_type = type
    current_data = select(results_epigenome, CpG_name, ends_with(type))
    colnames(current_data) = c("CpG_name", "type")
    # set up 2x2 table for Fisher's exact test
    # clock cpg == TRUE, enrichment category == TRUE
    a = nrow(filter(current_data, CpG_name %in% current_cgset & type == TRUE))
    # clock cpg == TRUE, enrichment category == FALSE
    b = nrow(filter(current_data, CpG_name %in% current_cgset & type == FALSE))
    # clock cpg == FALSE, enrichment category == TRUE
    c = nrow(filter(current_data, !(CpG_name %in% current_cgset) & type == TRUE))
    # clock cpg == FALSE, enrichment category == FALSE
    d = nrow(filter(current_data, !(CpG_name %in% current_cgset) & type == FALSE))
    # Fisher's exact test
    test = fisher.test(matrix(c(a,b,c,d), ncol=2, byrow=T))
    # extract coefficients 
    current_results$size_overlap = a
    current_results$ratio_overlap = a/(a+b)
    current_results$or = test$estimate
    current_results$confint_or_lower = test$conf.int[1]
    current_results$confint_or_upper = test$conf.int[2]
    current_results$p = test$p.value
    ## append to results
    results_enrichment <- rbind(results_enrichment, current_results)
  }
}

## FDR-correct & round values
results_enrichment$p_fdr <- p.adjust(results_enrichment$p, method = "fdr")
results_enrichment <- round_values(results_enrichment)

## Export results
save(results_enrichment, file = "03_results/enrichment_clock_cpgs.Rdata")

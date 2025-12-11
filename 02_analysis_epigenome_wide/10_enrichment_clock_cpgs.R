### Title: "Limited generalizability epigenome-wide: Enrichment of clock CpGs for cross-tissue CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-12"
### Purpose: Test CpGs on chronological & biological clocks for enrichment for subset of variable and tissue-conserved CpGs

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

## Load data
# epigenome-wide results
load("03_results/results_epigenome_wide.Rdata")
# grimage cpgs 
# retrieve from: https://github.com/bio-learn/biolearn/raw/refs/heads/master/biolearn/data/GrimAgeV1.csv
grimage_link <- "https://github.com/bio-learn/biolearn/raw/refs/heads/master/biolearn/data/GrimAgeV1.csv"
grimage_data <- read.csv(grimage_link)

### Prepare data
## get chronological clock CpGs (Horvath, SkinBlood, Hannum, Wu, PedBE)
horvath <- filter(methylclockData::get_coefHorvath(), CpGmarker != "(Intercept)")$CpGmarker
skinblood <- filter(methylclockData::get_coefSkin(), CpGmarker != "(Intercept)")$CpGmarker
hannum <- filter(methylclockData::get_coefHannum(), CpGmarker != "(Intercept)")$CpGmarker
wu <- filter(methylclockData::get_coefWu(), CpGmarker != "(Intercept)")$CpGmarker
pedbe <- filter(methylclockData::get_coefPedBE(), CpGmarker != "(Intercept)")$CpGmarker

chrono_cpgs <- intersect(unique(c(horvath, skinblood, hannum, wu, pedbe)), results_epigenome$CpG_name)

## get biological clock CpGs (PhenoAge, GrimAge, DunedinPACE)
phenoage <- filter(methylclockData::get_coefLevine(), CpGmarker != "(Intercept)")$CpGmarker
grimage <- filter(grimage_data, grepl("cg", var) == TRUE)$var
dunedin <- DunedinPACE::getRequiredProbes()$DunedinPACE

bio_cpgs <- intersect(unique(c(phenoage, grimage, dunedin)), results_epigenome$CpG_name)

### Test for enrichment
# test a) chronological clock cpg and b) biological clock cpg enrichment for:
# 1) variable cpgs, 2) variable and highly tissue-correlated cpgs, 3) variable and highly tissue-correlated cpgs with celltype adjustment

## Chronological clock CpGs
# set up empty results object
results_chrono <- as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(results_chrono) <- c("clock_type", "enrichment_type", "size_overlap", "ratio_overlap", "or", "or_confint_lower", "or_confint_upper", "p")
types <- c("variable_cpg", "variable_and_correlated", "variable_and_correlated_celltype_adjusted")

for(i in 1: length(types)){
  type = types[i]
  print(paste0("Testing for: ", type))
  # set up empty results object
  current_results <- as.data.frame(matrix(nrow = 1, ncol = 0))
  current_results$clock_type = "chronological"
  current_results$enrichment_type = type
  current_data = select(results_epigenome, CpG_name, ends_with(type))
  colnames(current_data) = c("CpG_name", "type")
  # set up 2x2 table for Fisher's exact test
  # clock cpg == TRUE, enrichment category == TRUE
  a = nrow(filter(current_data, CpG_name %in% chrono_cpgs & type == TRUE))
  # clock cpg == TRUE, enrichment category == FALSE
  b = nrow(filter(current_data, CpG_name %in% chrono_cpgs & type == FALSE))
  # clock cpg == FALSE, enrichment category == TRUE
  c = nrow(filter(current_data, !(CpG_name %in% chrono_cpgs) & type == TRUE))
  # clock cpg == FALSE, enrichment category == FALSE
  d = nrow(filter(current_data, !(CpG_name %in% chrono_cpgs) & type == FALSE))
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
  results_chrono <- rbind(results_chrono, current_results)
}

## Biological clock CpGs
results_bio <- as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(results_bio) <- c("clock_type", "enrichment_type", "size_overlap", "ratio_overlap", "or", "or_confint_lower", "or_confint_upper", "p")
types <- c("variable_cpg", "variable_and_correlated", "variable_and_correlated_celltype_adjusted")

for(i in 1: length(types)){
  type = types[i]
  print(paste0("Testing for: ", type))
  # set up empty results object
  current_results <- as.data.frame(matrix(nrow = 1, ncol = 0))
  current_results$clock_type = "biological"
  current_results$enrichment_type = type
  current_data = select(results_epigenome, CpG_name, ends_with(type))
  colnames(current_data) = c("CpG_name", "type")
  # set up 2x2 table for Fisher's exact test
  # clock cpg == TRUE, enrichment category == TRUE
  a = nrow(filter(current_data, CpG_name %in% bio_cpgs & type == TRUE))
  # clock cpg == TRUE, enrichment category == FALSE
  b = nrow(filter(current_data, CpG_name %in% bio_cpgs & type == FALSE))
  # clock cpg == FALSE, enrichment category == TRUE
  c = nrow(filter(current_data, !(CpG_name %in% bio_cpgs) & type == TRUE))
  # clock cpg == FALSE, enrichment category == FALSE
  d = nrow(filter(current_data, !(CpG_name %in% bio_cpgs) & type == FALSE))
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
  results_bio <- rbind(results_bio, current_results)
}

## Combine results and FDR-correct
results_clock_cpgs <- rbind(results_chrono, results_bio)
results_clock_cpgs$p_fdr <- p.adjust(results_clock_cpgs$p, method = "fdr")

## Export results
save(results_clock_cpgs, file = "03_results/enrichment_clock_cpgs.Rdata")

### Title: "Limited generalizability epigenome-wide: Blood-saliva correlations without cell type correction"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-01-30"
### Purpose: calculate Spearman correlations and ICCs between blood and saliva for each CpG, 
### test significance, export results
### ICC measure: ICC2 (random effects, absolute agreement)
### Note: to be run on HPC cluster

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(psych)
# arguments & name assignment
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()
chunk <- args[[1]]
dataset <- paste0("../../../01_data/02_proc/09_cross_tissue/70chunks/betas_filtered_crosstissue_", chunk, ".Rdata")
cpg_names <- paste0("../../../01_data/02_proc/09_cross_tissue/70chunks/common_cpgs_", chunk, ".Rdata")
print(paste0("Cross-tissue correlation analysis without cell type correction for chunk ", chunk, " , ", Sys.Date()))

### Load data
load(dataset)
load(cpg_names)
# load("../../../01_data/02_proc/09_cross_tissue/betas_filtered_crosstissue_test.Rdata")
# assign("betas", betas_filtered_crosstissue_test)
# load("../../../01_data/02_proc/09_cross_tissue/common_cpgs_test.Rdata")
# assign("cpg_names", common_cpgs_test)


print("Setup & data import complete---------------------------------------------")

### Correlation analysis
print(paste0("Correlation analysis for " , length(cpg_names), " CpGs"))
# create empty results object to append to
resultsobj <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(resultsobj) <- c("cpg_name", "mean_blood", "sd_blood", "mean_saliva", "sd_saliva", 
                          "cor_spearman", "p")

# loop cor.test over all cpgs and append results
for(i in 1:length(cpg_names)){
  # subset
  current_cpg = cpg_names[i]
  current_dat = select(betas, starts_with(current_cpg))
  # correlate
  cor = cor.test(current_dat[,1], current_dat[,2],
                 alternative = "two.sided", method = "spearman",
                 use = "pairwise", conf.level = 0.95)
  # extract results
  current_results = data.frame(matrix(nrow = 1, ncol = 0))
  current_results$cpg_name = current_cpg
  current_results$mean_blood = mean(select(current_dat, ends_with("_blood"))[,1], na.rm = TRUE)
  current_results$sd_blood = sd(select(current_dat, ends_with("_blood"))[,1], na.rm = TRUE)
  current_results$mean_saliva = mean(select(current_dat, ends_with("_saliva"))[,1], na.rm = TRUE)
  current_results$sd_saliva = sd(select(current_dat, ends_with("_saliva"))[,1], na.rm = TRUE)
  current_results$cor_spearman = cor$estimate
  current_results$p = cor$p.value
  # append results to main results object
  resultsobj <- rbind(resultsobj, current_results)
}

### ICC analysis
print(paste0("ICC analysis for " , length(cpg_names), " CpGs"))
# create empty results object to append to
icc_results <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(icc_results) <- c("cpg_name", "icc", "icc_lower", "icc_upper")
# loop ICC over all cpgs and append results
for(i in 1:length(cpg_names)){
  # subset
  current_cpg = cpg_names[i]
  current_dat = select(betas, starts_with(current_cpg))
  icc = ICC(current_dat)
  current_results = data.frame(matrix(nrow = 1, ncol = 0))
  current_results$cpg_name = current_cpg
  current_results$icc = filter(icc$results, type == "ICC2")$ICC
  current_results$icc_lower = filter(icc$results, type == "ICC2")$"lower bound"
  current_results$icc_upper = filter(icc$results, type == "ICC2")$"upper bound"
  icc_results <- rbind(icc_results, current_results)
  }

### Combine results
resultsobj <- full_join(resultsobj, icc_results, by = "cpg_name")

print("Analysis complete---------------------------------------------")

### Export results
filename =  paste0("03_results/epigenome_wide_crosstissue_nocelltype_", chunk, ".Rdata")
save(resultsobj, file = filename)
print(paste0("Results have been saved under ", filename))

### Total runtime
end_time <- Sys.time()
total_sec <- end_time - start_time
print(paste0("Runtime: ", total_sec))

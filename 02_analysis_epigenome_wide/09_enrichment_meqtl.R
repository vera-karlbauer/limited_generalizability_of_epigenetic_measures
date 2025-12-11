### Title: "Limited generalizability epigenome-wides: Enrichment (meQTLs)"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-04-10"
### Purpose: Test highly correlated & variable CpGs for meQTL enrichment 
### Note: to be run on HPC cluster

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(ggplot2)

### Load data
# epigenome-wide results
load("03_results/results_epigenome_wide.Rdata")
# meQTL database (http://mqtldb.godmc.org.uk)
meqtl <- read.csv("../../../01_data/02_proc/09_cross_tissue/assoc_meta_all.csv")
# hg38 Illumina manifest: 
manifest_hg38 <- read.csv("../../../01_data/02_proc/01_epic/illumina_epic_v2_info/EPIC-8v2-0_A1.csv", header = TRUE, skip = 7)

### Enrichment analysis

## prep meQTL data
# filter meqtls for loci in results
print(paste0(nrow(meqtl), " meQTLs pre filtering"))
meqtl_filtered <- meqtl %>%
  filter(cpg %in% results_epigenome$CpG_name)
print(paste0(nrow(meqtl_filtered), " meQTLs post filtering"))
# assess meqtl significance according to pval threshold of 1e-8 for cis and 1e-14 for trans
meqtl_filtered <- meqtl_filtered %>%
  mutate(signif_meqtl = case_when(cistrans == TRUE & pval < 1e-8 ~ TRUE, 
                                  cistrans == FALSE & pval < 1e-14 ~ TRUE,
                                  TRUE ~ FALSE))

## prep results data
# filter results for loci tested in meQTL analysis (excludes all EPICv2 only loci)
results_meqtl_filtered <- results_epigenome %>%
  filter(CpG_name %in% meqtl_filtered$cpg)

## add meqtl info
results_meqtl_filtered <- results_meqtl_filtered %>%
  mutate(is_signif_meqtl = case_when(CpG_name %in% filter(meqtl_filtered, signif_meqtl == TRUE)$cpg ~ TRUE,
         TRUE ~ FALSE))
print(paste0(nrow(results_meqtl_filtered), " CpGs tested for enrichment"))

### Test for enrichment (Fisher's exact test)
print("----------- Enrichment for meQTLs -----------")
print("Enrichment against variable CpGs only (no cell type adjustment)")
results_meqtl_filtered_variable <- results_meqtl_filtered %>%
  filter(variable_cpg == TRUE)
print(paste0(nrow(results_meqtl_filtered_variable), " CpGs tested for enrichment"))
# without cell type correction
print("Enrichment cross-table (no adjustment):")
table(results_meqtl_filtered_variable$variable_and_correlated, results_meqtl_filtered_variable$is_signif_meqtl)
print("Fisher's exact test for meQTL enrichment without cell type correction")
fisher.test(results_meqtl_filtered_variable$variable_and_correlated, results_meqtl_filtered_variable$is_signif_meqtl,
            alternative = "two.sided")
# with cell type correction
print("Enrichment cross-table (cell type adjustment):")
table(results_meqtl_filtered_variable$variable_and_correlated_celltype_adjusted, results_meqtl_filtered_variable$is_signif_meqtl)
print("Fisher's exact test for meQTL enrichment with cell type correction")
fisher.test(results_meqtl_filtered_variable$variable_and_correlated_celltype_adjusted, results_meqtl_filtered_variable$is_signif_meqtl,
            alternative = "two.sided")

### Export results
save(results_meqtl_filtered, file = "03_results/results_epigenome_meqtl.Rdata")

### Title: "Limited generalizability epigenome-wide: Combine and annotate results without cell type correction"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-02-20"
### Purpose: combine results for each chunk of data, add annotation, perform FDR p-value correction, export
### Note: to be run on HPC cluster

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
# set chunk size
chunks <- 70

### Import annotations 
# hg38 Illumina manifest: 
manifest_hg38 <- read.csv("../../../01_data/02_proc/01_epic/illumina_epic_v2_info/EPIC-8v2-0_A1.csv", header = TRUE, skip = 7)
# hg19 Gencodev26 annotation: https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg19.manifest.gencode.v26lift37.tsv.gz
annotation_hg19 <- read.table("../../../01_data/02_proc/01_epic/illumina_epic_v2_info/Wanding_Zhou_EPICv2_annotation/EPICv2.hg19.manifest.gencode.v26lift37.tsv", header = TRUE)

### Filter annotations
# select probe names on EPICv1 & 450k
manifest_hg38 <- manifest_hg38 %>%
  select(IlmnID, Name, Methyl450_Loci, EPICv1_Loci)
# select relevant hg19 annotation
annotation_hg19 <- annotation_hg19 %>%
  select(probeID, CpG_chrm, CpG_beg, CpG_end, genesUniq, geneNames)
  
### Import & append results
# import first chunk
filename <- "03_results/epigenome_wide_crosstissue_nocelltype_1.Rdata"
load(filename)
results_nocelltype_combined <- resultsobj
# import & append rest of chunks
for(chunk in 2:chunks){
  filename = paste0("03_results/epigenome_wide_crosstissue_nocelltype_", chunk, ".Rdata")
  load(filename)
  results_nocelltype_combined <- rbind(results_nocelltype_combined, resultsobj)
}

### FDR-correct p-values
results_nocelltype_combined$p_fdr <- p.adjust(results_nocelltype_combined$p, method = "fdr")

### Append annotation (hg19)
results_nocelltype_combined <- left_join(results_nocelltype_combined, annotation_hg19, join_by("cpg_name" == "probeID"))
results_nocelltype_combined <- left_join(results_nocelltype_combined, manifest_hg38, join_by("cpg_name" == "IlmnID"))
results_nocelltype_combined <- results_nocelltype_combined %>%
  dplyr::rename(probe_ID = cpg_name,
                CpG_name = Name) %>%
  relocate(p_fdr, .after = p) %>%
  relocate(CpG_name, .after = probe_ID) %>%
  relocate(CpG_chrm, .after = CpG_name) %>%
  relocate(CpG_beg, .after = CpG_chrm) %>%
  relocate(CpG_end, .after = CpG_beg)

### Export
save(results_nocelltype_combined, file = "03_results/results_nocelltype_combined_annotated.Rdata")

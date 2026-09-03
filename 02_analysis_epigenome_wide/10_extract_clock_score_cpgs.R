### Title: "Evaluating epigenetic clocks and scores: Cross-tissue correlations of clock/score CpGs"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2026-08-28"
### Purpose: Extract CpGs and weights for all tested epigenetic clocks/scores, and combine with cross-tissue correlation results for each CpG

### Setup
## General
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
# libraries
library(dplyr)
library(data.table)
library(writexl)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(methylclockData)
library(DunedinPACE)
library(methylCIPHER)
library(meffonym, verbose = FALSE)
library(readxl)
library(janitor)
library(ggdist)

## Load data
# epigenome-wide results
load("03_results/results_epigenome_wide.Rdata")

### Extract clock/score CpGs and weights
## Setup
# define clock order
clockorder <- c("Horvath", "PCHorvath", "SkinBlood", "PCSkinBlood", "Hannum", "PCHannum", 
                "Wu", "PedBE", "PhenoAge", "PCPhenoAge", "GrimAge", "PCGrimAge", "GrimAge2",
                "DunedinPACE", "BMI (Do, 2023)", "BMI (Wahl, 2017)",
                "CRP (Wielscher, 2022)", "CRP (Ligthart, 2016)", "Epigenetic-g",
                "DNAmTL", "PCDNAmTL", "Smoking exposure")
# empty clock list
clock_weights <- list()

## extract from methylclock (Horvath, SkinBlood, Hannum, Wu, PedBE, PhenoAge, DNAmTL)
clock_weights[["Horvath"]] <- get_coefHorvath() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["SkinBlood"]] <- get_coefSkin() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["Hannum"]] <- get_coefHannum() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["Wu"]] <- get_coefWu() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["PedBE"]] <- get_coefPedBE() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["PhenoAge"]] <- get_coefLevine() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)
clock_weights[["DNAmTL"]] <- get_coefTL() %>% rename(cpg = CpGmarker, weight = CoefficientTraining)

## DunedinPACE
clock_weights[["DunedinPACE"]] <- as.data.frame(cbind(DunedinPACE::mPACE_Models$model_probes$DunedinPACE, 
                                                      DunedinPACE::mPACE_Models$model_weights$DunedinPACE)) %>%
  rename(cpg = V1, weight = V2) %>%
  mutate(weight = as.numeric(weight))

## PC clocks
# load PC clock data (from author source code: https://github.com/MorganLevineLab/PC-Clocks)
load("../../../../tools/pc_clocks_no_imputation/CalcAllPCClocks.RData")
# list object names
pc_object_candidates <- c(PCHorvath = "CalcPCHorvath1",
                          PCSkinBlood = "CalcPCHorvath2",
                          PCHannum = "CalcPCHannum",
                          PCPhenoAge = "CalcPCPhenoAge",
                          PCGrimAge = "CalcPCGrimAge",
                          PCDNAmTL = "CalcPCDNAmTL")
# build extractor function for all PC clocks except PCGrimAge
extract_pc_clock_weights <- function(pc_obj, clock_name = "") {
  stopifnot(!is.null(pc_obj$model), !is.null(pc_obj$rotation))
  rot_raw <- pc_obj$rotation
  if (is.list(rot_raw) && !is.matrix(rot_raw) && !is.data.frame(rot_raw) && length(rot_raw) == 2) {
    cpg_names <- rot_raw[[1]]
    loadings  <- as.matrix(rot_raw[[2]])
    stopifnot(length(cpg_names) == nrow(loadings))
    rownames(loadings) <- cpg_names
    rot <- loadings
  } else if (is.data.frame(rot_raw)) {
    name_col <- intersect(c("cpg", "CpG", "CpGmarker", "cg", "probeID"), colnames(rot_raw))
    stopifnot(length(name_col) >= 1)
    cpg_names <- rot_raw[[name_col[1]]]
    rot <- as.matrix(rot_raw[, setdiff(colnames(rot_raw), name_col[1])])
    rownames(rot) <- cpg_names
  } else if (is.matrix(rot_raw) && !is.null(rownames(rot_raw))) {
    rot <- rot_raw
  } else {
    stop(sprintf("%s: unrecognized rotation structure — inspect str(pc_obj$rotation) manually", clock_name))
  }
  mod <- pc_obj$model
  if (is.matrix(mod) || inherits(mod, "dgCMatrix")) {
    coef_vec <- as.numeric(mod)
    names(coef_vec) <- rownames(mod)
  } else {
    coef_vec <- mod
  }
  coef_vec <- coef_vec[!grepl("intercept", names(coef_vec), ignore.case = TRUE)]
  common_pcs <- intersect(colnames(rot), names(coef_vec))
  if (length(common_pcs) == 0) {
    stop(sprintf("%s: no overlapping PC names between rotation and model", clock_name))
  }
  if (length(common_pcs) < length(coef_vec)) {
    message(sprintf("%s: using %d/%d PCs present in both rotation and model",
                    clock_name, length(common_pcs), length(coef_vec)))
  }
  rot      <- rot[, common_pcs, drop = FALSE]
  coef_vec <- coef_vec[common_pcs]
  eff_weight <- as.numeric(rot %*% coef_vec)
  names(eff_weight) <- rownames(rot)
  tibble(cpg = names(eff_weight), weight = eff_weight)
}
# build extractor for PCGrimAge (considering biomarker composites)
extract_pc_grimage_weights <- function(pc_obj) {
  rot <- pc_obj$rotation
  if (is.null(rownames(rot))) {
    dn <- dimnames(rot)
    rownames(rot) <- dn[[1]]
    colnames(rot) <- dn[[2]]
  }
  final_coefs <- pc_obj[["PCGrimAge.model"]]
  components  <- pc_obj[["components"]]
  submodel_components <- components[paste0(components, ".model") %in% names(pc_obj)]
  eff_weight_total <- setNames(numeric(nrow(rot)), rownames(rot))
  for (comp in submodel_components) {
    submodel <- pc_obj[[paste0(comp, ".model")]]
    all_terms <- names(submodel)
    pcs         <- intersect(all_terms, colnames(rot))
    non_pc_terms <- setdiff(all_terms, colnames(rot))
    if (length(non_pc_terms) > 0) {
      message(sprintf("Component '%s': dropping %d non-CpG covariate term(s): %s",
                      comp, length(non_pc_terms), paste(non_pc_terms, collapse = ", ")))
    }
    if (length(pcs) == 0) {
      warning(sprintf("Component '%s': no PC terms matched rotation matrix — skipping.", comp))
      next
    }
    eff_weight_comp <- as.numeric(rot[, pcs, drop = FALSE] %*% submodel[pcs])
    final_name <- paste0("DNAm", sub("^PC", "", comp))
    
    if (!final_name %in% names(final_coefs)) {
      warning(sprintf("Component '%s': no matching final coefficient '%s' found.", comp, final_name))
      next
    }
    final_coef <- final_coefs[[final_name]]
    eff_weight_total <- eff_weight_total + final_coef * eff_weight_comp
    
    message(sprintf("PCGrimAge: added component '%s' (%d/%d PC terms used, final coef = %.4g)",
                    comp, length(pcs), length(all_terms), final_coef))
  }
  tibble(cpg = names(eff_weight_total), weight = eff_weight_total)
}
# extract PC clock weights
pc_clock_weights <- purrr::imap(pc_object_candidates, function(obj_name, display_name) {
  if (!exists(obj_name)) {
    warning(sprintf("Object '%s' not found for %s — skipping.", obj_name, display_name))
    return(NULL)
  }
  obj <- get(obj_name)
  
  if (display_name == "PCGrimAge") {
    extract_pc_grimage_weights(obj)
  } else {
    extract_pc_clock_weights(obj, clock_name = display_name)
  }
})
names(pc_clock_weights) <- names(pc_object_candidates)
# append PC clock weights to weight list
for (nm in names(pc_clock_weights)) {
  clock_weights[[nm]] <- pc_clock_weights[[nm]]
  }

## epigenetic clocks & scores included in meffonym
meffonym_scores <- c("grimage", "grimagev2", "wahl-bmi", "do.bmi", "gs-cognitive")
for(i in 1:length(meffonym_scores)) {
  clock = meffonym_scores[i]
  cpg = meffonym.get.model(clock)$vars[!(meffonym.get.model(clock)$vars %in% c("intercept"))]
  weight = meffonym.get.model(clock)$coefs
  clock_weights[[clock]] <- as.data.frame(cbind(cpg, weight)) %>%
    dplyr::filter(!(cpg %in% c("female", "age"))) %>%
    mutate(weight = as.numeric(weight))
}
names(clock_weights)[names(clock_weights) == "grimage"] <- "GrimAge"
names(clock_weights)[names(clock_weights) == "grimagev2"] <- "GrimAge2"
names(clock_weights)[names(clock_weights) == "wahl-bmi"] <- "BMI (Wahl, 2017)"
names(clock_weights)[names(clock_weights) == "do.bmi"] <- "BMI (Do, 2023)"
names(clock_weights)[names(clock_weights) == "gs-cognitive"] <- "Epigenetic-g"

## epigenetic CRP
# Wielscher
crp_wielscher <- read_excel("../../01_preproc/epi_scores/crp_wielscher_2022/Wielscher_et_al_2022_CRP_EWAS_supplementary.xlsx",
                            sheet = 3) %>%
  filter(indepenent_loci == "1")
clock_weights[["CRP (Wielscher, 2022)"]] <- dplyr::select(crp_wielscher, ID, Effect_metaAnalysis) %>%
  rename("cpg" = ID, "weight" = Effect_metaAnalysis) %>%
  filter(grepl("intercept", cpg) == FALSE) %>%
  mutate(weight = as.numeric(weight))
# Ligthart
crp_ligthart <- read_excel("../../01_preproc/epi_scores/crp_ligthart_2016/Ligthart_et_al_2016_supplement_CRP.xlsx") %>%
  row_to_names(row_number = 1, remove_row = TRUE, remove_rows_above = TRUE) %>%
  clean_names()
clock_weights[["CRP (Ligthart, 2016)"]] <- dplyr::select(crp_ligthart, marker_name, effect) %>%
  rename("cpg" = marker_name, "weight" = effect) %>%
  filter(grepl("intercept", cpg) == FALSE) %>%
  mutate(weight = as.numeric(weight))

## epigenetic prenatal smoking score
smoking <- read.csv2("../../01_preproc/epi_scores/prenatal_smoking/weights_prenatal_smoking_exposure_richmond.csv")
clock_weights[["Smoking exposure"]] <- dplyr::select(smoking, CG, Weight) %>%
  rename("cpg" = CG, "weight" = Weight) %>%
  filter(grepl("intercept", cpg) == FALSE) %>%
  mutate(weight = as.numeric(weight))

## export list of clock/score CpGs and corresponding weights
save(clock_weights, file = "03_results/clock_score_cpg_weights.Rdata")

### Combine extracted weights with cross-tissue correlations
## build extractor/combination function
extract_clock_cpgs <- function(clock_weights, results_epigenome, normalize_weight = TRUE) {
  purrr::imap_dfr(clock_weights, function(tbl, clock_name) {
    if (is.null(tbl)) return(NULL)
    tbl <- tbl %>% distinct(cpg, .keep_all = TRUE)
    n_cpgs_val <- nrow(tbl)
    joined <- tbl %>%
      inner_join(results_epigenome, by = c("cpg" = "CpG_name"))
    n_missing <- n_cpgs_val - nrow(joined)
    if (n_missing > 0) {
      message(sprintf("%s: %d/%d CpGs not found in results_epigenome",
                      clock_name, n_missing, n_cpgs_val))
    }
    joined %>%
      mutate(
        abs_weight = abs(weight),
        n_cpgs = n_cpgs_val,  
        weight_share = abs_weight / sum(abs_weight),
        relative_weight = weight_share * n_cpgs
      ) %>%
      { if (normalize_weight)
        mutate(., abs_weight_scaled = ((relative_weight - min(relative_weight))
                                       /(max(relative_weight) - min(relative_weight))))
        else mutate(., abs_weight_scaled = relative_weight) } %>%
      select(cpg, weight, abs_weight, n_cpgs, relative_weight, abs_weight_scaled,
             cor_spearman, cor_spearman_celltype_adjusted, icc) %>%
      pivot_longer(cols = c(cor_spearman, cor_spearman_celltype_adjusted),
                   names_to = "adjustment", values_to = "correlation") %>%
      mutate(adjustment = dplyr::recode(adjustment,
                                        cor_spearman = "unadjusted",
                                        cor_spearman_celltype_adjusted = "cell type adjusted"),
             clock = clock_name)
  })
}

## for each epigenetic clock/score CpGs, extract blood-saliva correlations
crosstissue_clockscore_cpgs <- extract_clock_cpgs(clock_weights, results_epigenome)
crosstissue_clockscore_cpgs <- crosstissue_clockscore_cpgs %>%
  mutate(clock = factor(clock, levels = rev(clockorder)),
         adjustment = factor(adjustment, levels = c("unadjusted", "cell type adjusted")))

## export combined results
save(crosstissue_clockscore_cpgs, file = "03_results/cross_tissue_correlations_clock_score_cpgs.Rdata")

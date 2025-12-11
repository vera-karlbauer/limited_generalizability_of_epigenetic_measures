### Title: "Cross-tissue correlations: Define functions"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-25"
### Purpose: Define functions for cross-tissue project

### Function 1: Rename epigenetic clock from variable name format to plotting format (e.g. "epiage_horvath_blood" to "Horvath")

#' rename_clocks
#'
#' @param names a vector of clock names in variable name format
#'
#' @return vector of clock names in plotting/publication format
#' @export

rename_clocks <- function(names){
  new_names = str_replace_all(names, 
                              c("_blood" = "",
                                "_saliva" = "",
                                "_resid" = "",
                                "_T0" = "",
                                "_T2" = "",
                                "epiage_" = "", "_blood" = "", "_saliva" = "",
                                "pc_horvath" = "PCHorvath",
                                "horvath" = "Horvath",
                                "pc_skinblood" = "PCSkinBlood",
                                "skinblood" = "SkinBlood",
                                "pc_hannum" = "PCHannum",
                                "hannum" = "Hannum",
                                "pedbe" = "PedBE",
                                "pc_phenoage" = "PCPhenoAge",
                                "phenoage" = "PhenoAge",
                                "wu" = "Wu",
                                "grimage2" = "GrimAge2",
                                "pc_grimage" = "PCGrimAge",
                                "grimage" = "GrimAge",
                                "dunedin_pace" = "DunedinPACE"))
  return(new_names)
}

### Function 2: Rename age acceleration marker from variable name format to plotting format (e.g. "epiage_accel_horvath_blood" to "HorvathAccel")

#' rename_accel
#'
#' @param names a vector of clock names in variable name format
#'
#' @return vector of clock names in plotting/publication format
#' @export
rename_accel <- function(names){
  new_names = str_replace_all(names, 
                              c("_blood" = "",
                                "_saliva" = "",
                                "_resid" = "",
                                "_T0" = "",
                                "_T2" = "",
                                "epiage_accel_" = "", "_blood" = "", "_saliva" = "",
                                "pc_horvath" = "PCHorvathAccel",
                                "horvath" = "HorvathAccel",
                                "pc_skinblood" = "PCSkinBloodAccel",
                                "skinblood" = "SkinBloodAccel",
                                "pc_hannum" = "PCHannumAccel",
                                "hannum" = "HannumAccel",
                                "pedbe" = "PedBEAccel",
                                "pc_phenoage" = "PCPhenoAgeAccel",
                                "phenoage" = "PhenoAgeAccel",
                                "wu" = "WuAccel",
                                "grimage2" = "GrimAge2Accel",
                                "pc_grimage" = "PCGrimAgeAccel",
                                "grimage" = "GrimAgeAccel"))
  return(new_names)
}


### Function 3: Rename age delta marker from variable name format to plotting format (e.g. "epiage_delta_horvath_blood" to "HorvathDelta")

#' rename_delta
#'
#' @param names a vector of clock names in variable name format
#'
#' @return vector of clock names in plotting/publication format
#' @export
rename_delta <- function(names){
  new_names = str_replace_all(names, 
                              c("_blood" = "",
                                "_saliva" = "",
                                "_resid" = "",
                                "_T0" = "",
                                "_T2" = "",
                                "epiage_delta_" = "", "_blood" = "", "_saliva" = "",
                                "pc_horvath" = "PCHorvathAccel",
                                "horvath" = "HorvathAccel",
                                "pc_skinblood" = "PCSkinBloodAccel",
                                "skinblood" = "SkinBloodAccel",
                                "pc_hannum" = "PCHannumAccel",
                                "hannum" = "HannumAccel",
                                "pedbe" = "PedBEAccel",
                                "pc_phenoage" = "PCPhenoAgeAccel",
                                "phenoage" = "PhenoAgeAccel",
                                "wu" = "WuAccel",
                                "grimage2" = "GrimAge2Accel",
                                "pc_grimage" = "PCGrimAgeAccel",
                                "grimage" = "GrimAgeAccel"))
  return(new_names)
}


### Function 4: Rename epigenetic scores from variable name format to plotting format (e.g. "epigenetic_bmi_score_do_2023" to "BMI (Do, 2023)")

#' rename_scores
#'
#' @param names a vector of score names in variable name format
#'
#' @return vector of score names in plotting/publication format
#' @export
rename_scores <- function(names){
  new_names = str_replace_all(names, 
                              c("_blood" = "",
                                "_saliva" = "",
                                "_resid" = "",
                                "_T0" = "",
                                "_T2" = "",
                                "epigenetic_bmi_score_do_2023" = "BMI (Do, 2023)",
                                "epigenetic_bmi_score_wahl_2017" = "BMI (Wahl, 2017)",
                                "epigenetic_crp_score_ligthart_2016" = "CRP (Ligthart, 2016)",
                                "epigenetic_crp_score_wielscher_2022" = "CRP (Wielscher, 2022)",
                                "epigenetic_g" = "Epigenetic-g",
                                "epiage_pc_telo" = "PC DNAmTL",
                                "epiage_telo" = "DNAmTL",
                                "smoking_exposure_score" = "Smoking exposure"))
  return(new_names)
}

### Function 5: Rename phenotypes from variable name format to plotting format (e.g. "epiage_delta_horvath_blood" to "HorvathDelta")

#' rename_phenos
#'
#' @param names a vector of clock names in variable name format
#'
#' @return vector of clock names in plotting/publication format
#' @export
rename_phenos <- function(names){
  new_names = str_replace_all(names, 
                              c("ku_bmi" = "BMI",
                                "crp_pgML_ln" = "log-transformed salivary CRP",
                                "iq_sonr_nonverbal" = "nonverbal IQ (SON-R)",
                                "telomere_length_blood_ts_ratio" = "blood T/S ratio"))
  return(new_names)
}

### Function 6: Round all numeric values to 2 digits except for p-values. 
# P-values are rounded according to JAMA standards: p>0.99 ="p>.99), 2 digits up to p = 0.01; 3 digits up to 0.001, p<0.001 = "p<.001"
#' round_values
#'
#' @param df a data frame of numeric and/or non-mumeric values
#'
#' @return a data frame with all numeric values rounded according to JAMA standards
#' @export
round_values <- function(df) {
  format_p <- function(p) {
    out <- rep(NA_character_, length(p))
    na <- is.na(p)
    
    gt_099  <- !na & p > 0.99
    gt_001  <- !na & !gt_099 & p > 0.01
    gt_0001 <- !na & !gt_099 & !gt_001 & p > 0.001
    rest    <- !na & !gt_099 & !gt_001 & !gt_0001
    
    out[gt_099]  <- ">.99"
    out[gt_001]  <- round(p[gt_001], digits = 2)
    out[gt_0001] <- round(p[gt_0001], digits = 3)
    out[rest]    <- "<.001"
    out
    }
  res <- df
  for (nm in names(res)) {
    x <- res[[nm]]
    if (is.numeric(x)) {
      if (nm %in% c("p", "p_fdr", "p_celltype_adjusted", "p_fdr_celltype_adjusted",
                    "p_chrono_age", "p_chrono_age_celltype_adjusted", "p_pheno", "p_pheno_celltype_adjusted", 
                    "p_tissue", "p_tissue_celltype_adjusted", "p_interaction", "p_interaction_celltype_adjusted", 
                    "p_interaction_fdr", "p_interaction_fdr_celltype_adjusted")) {
        res[[nm]] <- format_p(x)          # character with p-formatting
      } else {
        res[[nm]] <- round(x, 2)          # numeric rounded to 2 decimals
      }
    }
  }
  res
}


### Function 7: extract model coefficients from coxme(lmekin) output 
# Based on: https://stackoverflow.com/questions/43720260/how-to-extract-p-values-from-lmekin-objects-in-coxme-package
#' extract_coxme_table
#'
#' @param mod a model derived from a coxme function
#'
#' @return a table of summary statistics derived from coxme output
#' @export
extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2)
  p <- signif(1 - pchisq((beta/se)^2, 1), 2)
  table = data.frame(cbind(beta, se, z, p))
  return(table)
}

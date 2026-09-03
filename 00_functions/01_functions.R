### Title: "Evaluating epigenetic clocks and scores: Define functions"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-25"
### Purpose: Define functions for analysis and data handling

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
#' @param names a vector of names in variable name format
#'
#' @return vector of clock in plotting/publication format
#' @export
rename_phenos <- function(names){
  new_names = str_replace_all(names, 
                              c("ku_bmi" = "BMI",
                                "crp_pgML_ln" = "log-transformed salivary CRP",
                                "iq_sonr_nonverbal" = "nonverbal IQ (SON-R)",
                                "telomere_length_blood_ts_ratio" = "blood T/S ratio"))
  return(new_names)
}

### Function 6: Round all numeric values to 2 digits unless zero, then round to first non-zero digit 
#' round_values
#'
#' @param df a data frame of numeric and/or non-numeric values
#'
#' @return a data frame with rounded numeric values
#' @export
round_values <- function(df) {
  round_one <- function(x) {
    ifelse(
      is.na(x) | x == 0 | is.infinite(x),
      x,
      ifelse(
        round(x, 2) != 0,
        round(x, 2),
        round(x, ceiling(-log10(abs(x))))
      )
    )
  }
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) round_one(col) else col
  })
  
  df
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

### Function 8: Rename cell types from variable name format to plotting/table format (e.g. "NK_blood" to "Natural killer cells (blood)")

#' rename_celltypes
#'
#' @param names a vector of names in variable name format
#'
#' @return vector of names in plotting/publication format
#' @export
rename_celltypes <- function(names){
  new_names = str_replace_all(names, 
                              c("Baso_blood" = "Basophils (blood)",
                                "Bmem_blood" = "Memory B-cells (blood)",
                                "Bnv_blood" = "Naïve B-cells (blood)",
                                "CD4Tmem_blood" = "Memory CD4+ T-cells (blood)",
                                "CD4Tnv_blood" = "Naïve CD4+ T-cells (blood)",
                                "CD8Tmem_blood" = "Memory CD8+ T-cells (blood)",
                                "CD8Tnv_blood" = "Naïve CD8+ T-cells (blood)",
                                "Eos_blood" = "Eosinophils (blood)",
                                "Mono_blood" = "Monocytes (blood)",
                                "Neu_blood" = "Neutrophils (blood)",
                                "NK_blood" = "Natural killer cells (blood)",
                                "Treg_blood" = "T-regulatory cells (blood)",
                                "Epithelial_saliva" = "Buccal epithelial cells (saliva)"))
  return(new_names)
}

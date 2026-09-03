### Title: "Limited generalizability epigenome-wide: Annotate EWAS catalog traits"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-05"
### Purpose: Manually EWAS catalog traits to trait categories used for enrichment analysis

### Setup
# general
rm(list=ls())
# output in non-scientific notation, 4 digits
options(scipen = 6, digits = 4)
set.seed(123)
# libraries
library(dplyr)
library(data.table)
library(stringr)

### Load data
# epigenome-wide results
load("03_results/results_epigenome_wide.Rdata")
# EWAS catalog studies (https://www.ewascatalog.org/download/)
ewascatalog_studies <- fread("../../../../tools/ewascatalog/ewascatalog-studies.txt")

### Annotate EWAS catalog
## manual annotation of trait categories
# show number of unique traits 
length(unique(ewascatalog_studies$Trait))
# annotate traits
ewascatalog_studies <- ewascatalog_studies %>%
  mutate(trait_category = case_when(
           str_detect(Trait, regex("\\b(smoking|cigarette|tobacco|cotinine)s?\\b", ignore_case = TRUE)) 
           ~ "Smoking",
           str_detect(Trait, regex("\\b(alcohol)\\b", ignore_case = TRUE)) 
           ~ "Alcohol",
           str_detect(Trait, regex("\\b(education|socioeconomic|economic|marital|insurance|income|migration|residence|educational attainment|single parent|residential|disadvantage)\\b", ignore_case = TRUE)) 
           ~ "Socioeconomic",
           str_detect(Trait, regex("\\b(pregnancy|birth|prenatal|pregnant|perinatal|infant|gestation|gestational|birthweight|fetal|in vitro|reproductive|breastfeeding|neonatal|feeding|subscapular|placenta|folate|maternal|trimester|skinfold|parity|nicu|preeclampsia|intrauterine)\\b", ignore_case = TRUE)) 
           ~ "Gestation/Birth",
           str_detect(Trait, regex("\\b(trauma|abuse|traumatic|stress|victimization|maltreatment|neglect|violence|deprivation|bullying|adversity|discrimination|instability)\\b", ignore_case = TRUE)) 
           ~ "Trauma/Stress",
           str_detect(Trait, regex("\\b(schizophrenia|psychosis|psychotic|bipolar|mdd|depression|depressive|anxiety|phobia|anorexia|ocd|adhd|autism|autistic|ptsd|post-traumatic|tic|aggressive|conduct|deficit|Attention deficit hyperactivity disorder)\\b", ignore_case = TRUE)) 
           ~ "Mental disorders",
           str_detect(Trait, regex("\\b(substance use|cannabis|cocaine|opioid|drug use)\\b", ignore_case = TRUE)) 
           ~ "Substance use",
           str_detect(Trait, regex("\\b(wellbeing|mental health)\\b", ignore_case = TRUE)) 
           ~ "Mental health/wellbeing",
           str_detect(Trait, regex("\\b(education|social|socioeconomic|economic|marital|insurance|income|migration|residence|educational attainment|single parent|residential|disadvantage)\\b", ignore_case = TRUE)) 
           ~ "Socioeconomic",
           str_detect(Trait, regex("\\b(dementia|alzheimer|parkinson|amyloid|palsy|huntington|neurodegenerative|down syndrome|downs syndrome|multiple sclerosis|braak|seizures|apoe)s?\\b", ignore_case = TRUE)) 
           ~ "Neurological",
           str_detect(Trait, regex("\\b(bmi|body mass index|weight|glucose|insulin|lipid|cholesterol|triglyceride|waist|diabetes|obese|obesity|fat|lipoprotein|postprandial lipemia|metabolic|ldl|plasma resistin|homa-ir|adiponectin|adipositas|apolipoprotein|proinsulin|vldl|hdl|idl|lp)s?\\b", ignore_case = TRUE)) 
           ~ "Metabolic",
           str_detect(Trait, regex("\\b(cognitive|cognition|intelligence|iq|memory|executive)\\b", ignore_case = TRUE)) 
           ~ "Cognitive",
           str_detect(Trait, regex("\\b(cardiovascular|heart|blood pressure|systolic|diastolic|cardio|hypertension|hypotension|infarction|stroke|cardiac|atrial|vascular|carotid|lvmi|rwt|pulse rate|cvd)\\b", ignore_case = TRUE)) 
           ~ "Cardiovascular",
           str_detect(Trait, regex("\\b(crp|c-reactive protein|inflammation|inflammatory|interleukin|tumour necrosis factor alpha|tumor necrosis factor alpha|immunoglobulin|ige|IgG glycosylation)\\b", ignore_case = TRUE)) 
           ~ "Immune markers",
           str_detect(Trait, regex("\\b(immune|autoimmune|rheuma|rheumatoid|arthritis|inflammatory bowel|chrohn|asthma|allergic|allergy|lupus|ulcerative cholitis|psoriasis|sjogrens syndrome|ulcerative colitis|atopy|eosinophilia)\\b", ignore_case = TRUE)) 
           ~ "Autoimmune disease",
           str_detect(Trait, regex("\\b(hiv|human immunodeficiency virus|covid)\\b", ignore_case = TRUE)) 
           ~ "Infectious disease",
           str_detect(Trait, regex("\\b(lung|copd|respiratory|respiration|fev|fvc|expiratory|forced vital capacity|bronchodilator|exhaled|fev1)\\b", ignore_case = TRUE)) 
           ~ "Pulmonary",
           str_detect(Trait, regex("\\b(diet|food|nutrient|micronutrient|consumption|intake|eating|linoleic acid|fatty acid|choline|betaine)s?\\b", ignore_case = TRUE)) 
           ~ "Diet/Nutrition",
           str_detect(Trait, regex("\\b(cancer|carcinoma|adenoma|tumour|tumor|adenocarcinoma|melanoma)s?\\b", ignore_case = TRUE)) 
           ~ "Cancer",
           str_detect(Trait, regex("\\b(expression)\\b", ignore_case = TRUE)) 
           ~ "Gene expression",
           str_detect(Trait, regex("\\b(protein levels)\\b", ignore_case = TRUE)) 
           ~ "Protein",
           str_detect(Trait, regex("\\b(mir)\\b", ignore_case = TRUE)) 
           ~ "microRNA expression",
           str_detect(Trait, regex("\\b(pain)\\b", ignore_case = TRUE)) 
           ~ "Pain",
           str_detect(Trait, regex("\\b(mortality)\\b", ignore_case = TRUE)) 
           ~ "Mortality",
           str_detect(Trait, regex("\\b(ancestry|ethnicity)\\b", ignore_case = TRUE)) 
           ~ "Ancestry/Ethnicity",
           str_detect(Trait, regex("\\b(telomere)\\b", ignore_case = TRUE)) 
           ~ "Telomere length",
           str_detect(Trait, regex("\\b(exposure|pollution|arsenic|sunlight|altitude|atmospheric|particulate|immunoprophylaxis|cadmium|pm2.5|nox|dicamba|heptachlor|mesotrione|picloram|acetochlor|diphenyl|exposure|dichlorodiphenyldichloroethylene|copper|dichlorophenoxyacetic|aldrin|atrazine|chlordane|ddt|dieldrin|glyphosate|lindane|malathion|metolachlor|toxaphen|polution|toxaphene|MEP|EtP|PrP|n-BuP|BPA|BPF|MiBP|MnBP|MEHP|MEHHP|MEOHP|MECPP|MBzP)\\b", ignore_case = TRUE)) 
           ~ "Environmental exposure",
           str_detect(Trait, regex("\\b(vegf|hba1c|growth differentiation factor 15|conjugated dienes|glutathione|homocysteine|trimethylaminenoxide|transferase|citrate|aminotransferase|diacylglycerol|ctx|hormone|vitamin|selenium|coagulation|triiodothyronine|thyrotropin)\\b", ignore_case = TRUE)) 
           ~ "Blood biomarker",
           str_detect(Trait, regex("\\b(kidney|filtration|creatinine|IgA nephropathy)\\b", ignore_case = TRUE)) 
           ~ "Kidney",
           str_detect(Trait, regex("\\b(puberty|pubertal|menarche)\\b", ignore_case = TRUE)) 
           ~ "Puberty",
           str_detect(Trait, regex("\\b(cortical|thalamus|hippocampus)\\b", ignore_case = TRUE)) 
           ~ "Brain imaging",
           str_detect(Trait, regex("\\b(mass|density|circumference|height|length|bone|area|body)\\b", ignore_case = TRUE)) 
           ~ "Body composition",
           str_detect(Trait, regex("\\b(physical activity|sitting)\\b", ignore_case = TRUE)) 
           ~ "Physical activity",
           PMID == "24014485" ~ "Blood biomarker",
           PMID == "33413638" ~ "Blood biomarker",
           Date == "2020-08-20" ~ "Blood biomarker",
           str_detect(Trait, regex("\\b(disease|disorder|cirrhosis|osteoarthritis|chronic|cleft)\\b", ignore_case = TRUE)) 
           ~ "Other disease/disorder",
           str_detect(Trait, regex("\\b(age|aging|ageing)\\b", ignore_case = TRUE)) 
           ~ "Age",
           str_detect(Trait, regex("\\b(sex)\\b", ignore_case = TRUE)) 
           ~ "Sex",
           TRUE ~ NA))

# show number of studies without trait annotation
ewascatalog_studies %>% filter(is.na(trait_category)) %>% nrow()

# show number of unique trait categories
length(unique(ewascatalog_studies$trait_category))

## Only keep blood-based and salivary EWAS
# show all unique tissues
unique(ewascatalog_studies$Tissue)
# convert blood cell types or BECs to saliva
ewascatalog_studies <- ewascatalog_studies %>%
  mutate(Tissue = case_when(Tissue %in% c("Lymphocytes", "CD4+ T-cells", "Monocytes", "T cells", 
                                          "Whole blood, CD4+ T cells", "CD4+ T-cells, whole blood", "CD4+ T-cells, Whole blood",
                                          "CD4+ T cells, monocytes", "CD14+ monocytes", "Whole blood, cord blood", "Maternal whole blood",
                                          "B cells", "Peripheral blood mononuclear cells", "Whole blood, breast tissue", "Cord blood, whole blood",
                                          "Leukocytes", "Whole blood, CD4+ T cells, CD14+ monocytes", "whole blood", "Whole blood, CD4+ T-cells, CD14+ monocytes",
                                          "Whole blood, heel prick blood spots", "CD19+ B cells", "Peripheral blood", "blood") ~ "Blood",
                            Tissue == "Whole blood" ~ "Blood", 
                            Tissue %in% c("Buccal") ~ "Saliva",
                            TRUE ~ Tissue))
# filter and only keep blood&saliva
ewascatalog_studies <- ewascatalog_studies %>%
  filter(Tissue == "Blood" | Tissue == "Saliva")

# restrict to studies in children
ewascatalog_studies <- ewascatalog_studies %>%
  filter(Age == "Children")

# export catalogue with annotation
write.table(ewascatalog_studies, file = "../../../../tools/ewascatalog/ewascatalog_studies_annotated.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

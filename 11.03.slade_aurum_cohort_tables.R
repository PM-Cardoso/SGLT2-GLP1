####################
## Description:
##  - In this file we:
##    - Produce tables to describe cohorts
####################


library(tidyverse)
library(tableone)


source("11.02.slade_aurum_set_data.R")

####

# Inclusion criteria

set_up_data_sglt2_glp1(dataset.type = "diagnostics")

###

vars <- c("agetx", "t2dmduration", "sex", "ethnicity", "drug_canagliflozin", "drug_dapagliflozin", "drug_empagliflozin", "drug_ertugliflozin",
          "drug_dulaglutide", "drug_exenatide_short", "drug_exenatide_long", "drug_liraglutide", "drug_lixisenatide", "deprivation", "smoke", "drugline", 
          "ncurrtx", "MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1", "prehba1c", "prehba1c_na", "prebmi", "prebmi_na", "preegfr", "preegfr_na", 
          "prehdl", "prehdl_na", "prealt", "prealt_na", "prealbuminblood", "prealbuminblood_na", "prebilirubin", "prebilirubin_na", 
          "pretotalcholesterol", "pretotalcholesterol_na", "premap", "premap_na", "prediabeticnephropathy", "preneuropathy", "preretinopathy", "preangina",
          "ASCVD", "preaf", "prerevasc", "preheartfailure", "prehypertension", "preihd", "premyocardialinfarction", "prepad", "prestroke",
          "pretia", "preckd", "precld", "posthba1cfinal", "posthba1cfinal_na", "hba1cmonth", "hba1cmonth_na")

####

# Characteristics of individuals predicted benefit >5
hba1c.train.benefit <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  # left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  #             select(patid, pated, ethnicity), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_)),
         CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No")) %>%
  mutate(benefit = factor(benefit),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))



tab.benefits <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "benefit", data = hba1c.train.benefit, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


#####################
# Full cohort
full.cohort <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))
  

## Construct a table
tab.full.cohort <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = full.cohort, test = FALSE)
## Show table with SMD
print(tab.full.cohort, smd = TRUE)

# summary(tab.full.cohort)


#####################
#### PS model train
ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train") %>%
  select(patid, pated) %>%
  left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort"), by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))


## Construct a table
tab.ps.model.train <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = ps.model.train, test = FALSE)
## Show table with SMD
print(tab.ps.model.train, smd = TRUE)

# summary(tab.ps.model.train)


#####################
#### PS model test
ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test") %>%
  select(patid, pated) %>%
  left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort"), by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))


## Construct a table
tab.ps.model.test <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = ps.model.test, test = FALSE)
## Show table with SMD
print(tab.ps.model.test, smd = TRUE)

# summary(tab.ps.model.test)


#####################
#### HbA1c model train
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  select(patid, pated) %>%
  left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort"), by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))


## Construct a table
tab.hba1c.train <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = hba1c.train, test = FALSE)
## Show table with SMD
print(tab.hba1c.train, smd = TRUE)

# summary(tab.hba1c.train)



#####################
#### HbA1c model test
hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  select(patid, pated) %>%
  left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort"), by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
         preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
         drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
         drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
         drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
         drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
         drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
         drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
         drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
         drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
         drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
         prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
         prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
         preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
         prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
         prealt_na = ifelse(is.na(prealt), "Yes", "No"),
         prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
         prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
         pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
         premap_na = ifelse(is.na(premap), "Yes", "No"),
         posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
         hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))


## Construct a table
tab.hba1c.test <- CreateTableOne(vars = vars, factorVars = c("MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1"), includeNA = TRUE, strata = "drugclass", data = hba1c.test, test = FALSE)
## Show table with SMD
print(tab.hba1c.test, smd = TRUE)

# summary(tab.hba1c.test)

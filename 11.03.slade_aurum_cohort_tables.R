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


####

# Characteristics of individuals predicted benefit >5
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  # left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  #             select(patid, pated, ethnicity), by = c("patid", "pated")) %>%
  mutate(benefit = ifelse(effects < -5, "SGLT2i", ifelse(effects > 5, "GLP1-RA", NA_real_)),
         CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No")) %>%
  mutate(benefit = factor(benefit))


vars <- c("drugsubstances", "agetx", "sex", "t2dmduration", "ethnicity", "deprivation",
          "smoke", "prehospitalisation", "drugline", "ncurrtx", "hba1cmonth", "preacr", "prealbuminblood",
          "prealt", "preast", "prebilirubin", "prebmi", "prehaematocrit", "prehaemoglobin",
          "prehba1c", "prehdl", "premap", "preegfr", "pretotalcholesterol",
          "pretriglyceride", "preangina", "precld", "prediabeticnephropathy", 
          "preheartfailure", "prehypertension", "preihd", "premyocardialinfarction",
          "preneuropathy", "prepad", "preretinopathy", "prerevasc", "prestroke",
          "pretia", "preaf")
factorvars <- c("drugsubstances", "sex", "ethnicity", "deprivation", "smoke", "drugline", "prehospitalisation", "ncurrtx", "preangina", 
                "precld", "prediabeticnephropathy", "preheartfailure", "prehypertension", "preihd",
                "premyocardialinfarction", "preneuropathy", "prepad", "preretinopathy", "prerevasc",
                "prestroke", "pretia", "preaf")


tab.benefits <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "benefit", data = hba1c.train, test = FALSE)
## Show table with SMD
print(tab.benefits, smd = TRUE)


############################
vars <- c("drugsubstances", "agetx", "sex", "t2dmduration", "ethnicity", "deprivation",
          "smoke", "prehospitalisation", "drugline", "ncurrtx", "hba1cmonth", "preacr", "prealbuminblood",
          "prealt", "preast", "prebilirubin", "prebmi", "prehaematocrit", "prehaemoglobin",
          "prehba1c", "prehdl", "premap", "preegfr", "pretotalcholesterol",
          "pretriglyceride", "preangina", "precld", "prediabeticnephropathy", 
          "preheartfailure", "prehypertension", "preihd", "premyocardialinfarction",
          "preneuropathy", "prepad", "preretinopathy", "prerevasc", "prestroke",
          "pretia", "preaf")
factorvars <- c("drugsubstances", "sex", "ethnicity", "deprivation", "smoke", "drugline", "prehospitalisation", "ncurrtx", "preangina", 
                "precld", "prediabeticnephropathy", "preheartfailure", "prehypertension", "preihd",
                "premyocardialinfarction", "preneuropathy", "prepad", "preretinopathy", "prerevasc",
                "prestroke", "pretia", "preaf")


#####################
# Full cohort
full.cohort <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")))
  

## Construct a table
tab.full.cohort <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), factorVars = c(factorvars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "drugclass", data = full.cohort, test = FALSE)
## Show table with SMD
print(tab.full.cohort, smd = TRUE)

# summary(tab.full.cohort)


#####################
#### PS model train
ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train") %>%
  select(patid, pated) %>%
  left_join(full.cohort, by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")))


## Construct a table
tab.ps.model.train <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), factorVars = c(factorvars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "drugclass", data = ps.model.train, test = FALSE)
## Show table with SMD
print(tab.ps.model.train, smd = TRUE)

# summary(tab.ps.model.train)


#####################
#### PS model test
ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test") %>%
  select(patid, pated) %>%
  left_join(full.cohort, by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")))


## Construct a table
tab.ps.model.test <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), factorVars = c(factorvars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "drugclass", data = ps.model.test, test = FALSE)
## Show table with SMD
print(tab.ps.model.test, smd = TRUE)

# summary(tab.ps.model.test)


#####################
#### HbA1c model train
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  select(patid, pated) %>%
  left_join(full.cohort, by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")))


## Construct a table
tab.hba1c.train <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), factorVars = c(factorvars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "drugclass", data = hba1c.train, test = FALSE)
## Show table with SMD
print(tab.hba1c.train, smd = TRUE)

# summary(tab.hba1c.train)



#####################
#### HbA1c model test
hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  select(patid, pated) %>%
  left_join(full.cohort, by = c("patid", "pated")) %>%
  mutate(CV_problems = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
         microvascular_complications = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
         ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
         deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")))


## Construct a table
tab.hba1c.test <- CreateTableOne(vars = c(vars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "posthba1cfinal", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), factorVars = c(factorvars, "preckd", "CV_problems", "microvascular_complications", "ASCVD", "MFN",  "DPP4", "GLP1", "SGLT2", "SU", "TZD"), includeNA = TRUE, strata = "drugclass", data = hba1c.test, test = FALSE)
## Show table with SMD
print(tab.hba1c.test, smd = TRUE)

# summary(tab.hba1c.test)

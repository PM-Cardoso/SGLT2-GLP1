####################
## Description:
##  - This file includes the framework for patient selection.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])

# load libraries
require(tidyverse)


###############################################################################
###############################################################################
############################### Functions #####################################
###############################################################################
###############################################################################


set_up_data <- function() {
  
  # load original dataset # name - t2d_dataset_converted
  load("/slade/CPRD_data/mastermind_2022/20221019_t2d_dataset.Rda")
  
  cprd <- t2d_dataset_converted
  
  
  ################################################
  ##### Select only SGLT-2 and GLP-1
  ################################################
  
  #SGLT2 vs GLP1
  cprd <- cprd %>% 
    filter(drugclass == "GLP1" | drugclass == "SGLT2") # nrow = 160744
  
  # table(cprd$drugclass)
  # # GLP1  SGLT2
  # # 60627 100117
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  cprd <- cprd %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # table(cprd$timeprevcombo_less61) 
  # # removing 26132
  
  cprd <- cprd %>% 
    filter(is.na(timeprevcombo_less61)) %>%
    select(-timeprevcombo_less61)

  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # table(cprd$INS) 
  # # removing 23406
  
  cprd <- cprd %>% 
    filter(INS == 0)      
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # table(cprd$drugline)
  
  # Not needed since there are no patients with drugline = 1
  
  
  ################################################
  ##### Drop if within 1 year of diagnosis
  ################################################
  
  cprd <- cprd %>%
    mutate(dstart_1year = difftime(dstartdate, dm_diag_date, units = "days")) %>%
    mutate(dstart_1year = ifelse(dstart_1year < 365, 1, NA_real_))
  
  # table(cprd$dstart_1year)
  # # removing 2822
  
  cprd <- cprd %>%
    filter(is.na(dstart_1year)) %>%
    select(-dstart_1year)
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # cprd <- cprd %>%
  #   mutate(pated = paste(as.character(patid), as.character(dstartdate), sep = "."))
  # 
  # duplicated.pated <- unique(cprd$pated[duplicated(cprd$pated)])
  # 
  # cprd <- cprd %>%
  #   mutate(duplicated.pated = ifelse(pated %in% duplicated.pated, 1, NA_real_))
  # 
  # # table(cprd$duplicated.pated)
  # # # removing 220
  # 
  # cprd <- cprd %>%
  #   filter(is.na(duplicated.pated)) %>%
  #   select(-duplicated.pated)
  
  
  ################################################
  ##### Drop patients initiating before 1/1/2015
  ################################################
  
  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < "2015-01-01", 1 , NA_real_))
  
  # table(cprd$dstartdate_cutoff)
  # # removing 21666
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff)) %>%
    select(-dstartdate_cutoff)
  
  
  ################################################
  ##### Drop patients initiating within 91 days of registration
  ################################################
  
  # table(cprd$dm_diag_flag)
  # # removing 4306
  
  cprd <- cprd %>%
    filter(dm_diag_flag == 0)
  
  
  ################################################
  ##### Drop patients initiating drug whilst on the opposite
  ################################################
  
  # table(cprd$drugclass, cprd$SGLT2) # removing 3927
  # #            0     1
  # # GLP1  14829  3927
  # # SGLT2     0 63656
  

  # table(cprd$drugclass, cprd$GLP1) # removing 3445
  # #           0     1
  # # GLP1      0 18756
  # # SGLT2 60211  3445
  
  
  cprd <- cprd %>%
    mutate(opposite_drug = ifelse( (drugclass == "GLP1" & SGLT2 == 1) | (drugclass == "SGLT2" & GLP1 == 1) , 1, NA_real_)) %>%
    filter(is.na(opposite_drug)) %>%
    select(-opposite_drug)
  
  
  ################################################
  ##### Drop if HbA1c not in 53 - 120 at baseline + missing initial hba1c
  ################################################
  
  cprd <- cprd %>% 
    mutate(hb_extreme= ifelse(prehba1c < 53 | prehba1c > 120 | is.na(prehba1c), 1, 0))
  
  # table(cprd$hb_extreme)
  # # removing 7220
  
  cprd <- cprd %>% 
    filter(hb_extreme==0) %>%
    select(-hb_extreme)
  
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  ################################################
  ##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
  #####   - posthba1c12m but if missing
  #####     - posthba1c6m
  ################################################
  
  cprd <- cprd %>%
    mutate(posthba1c_final = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%
    mutate(posthba1c_final = as.numeric(posthba1c_final))
    
  
  ################################################
  ##### Sociodemographic variables
  #####   - Age: agetx (new var)
  ################################################
  
  cprd <- cprd %>%
    mutate(agetx = difftime(dstartdate, dob, units = "days") / 365.25) %>%
    mutate(agetx = as.numeric(agetx))
  
  ################################################
  ##### Sociodemographic variables
  #####   - Sex: malesex
  ################################################
  
  cprd <- cprd %>%
    mutate(malesex = ifelse(gender == 1, 1, 0)) %>%
    mutate(malesex = as.factor(malesex))
  
  ################################################
  ##### Sociodemographic variables
  #####   - Duration of diabetes: t2dmduration
  ################################################
  
  cprd <- cprd %>%
    mutate(t2dmduration = difftime(dstartdate, dm_diag_date, units = "days") / 365.25) %>%
    mutate(t2dmduration = as.numeric(t2dmduration))
  
  ################################################
  ##### Sociodemographic variables
  #####   - Ethnicity: Ethnicity
  ################################################
  
  cprd <- cprd %>%
    mutate(Ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed")))
  
  ################################################
  ##### Sociodemographic variables
  #####   - Deprivation: Deprivation
  ################################################
  
  cprd <- cprd %>%
    mutate(Deprivation = factor(imd2015_10))
  
  ################################################
  ##### Sociodemographic variables
  #####   - Smoking Status: Smoke
  ################################################
  
  ## Variable not existant
  
  ################################################
  ##### Diabetes treatment
  #####   - Line Therapy: drugline
  ################################################
  
  ## turn all > 4 to 5+
  cprd <- cprd %>%
    mutate(drugline = ifelse(drugline > 4, 5, drugline)) %>%
    mutate(drugline = as.factor(drugline))
  
  ################################################
  ##### Diabetes treatment
  #####   - Drugs taken alongside treatment
  #####     - SU
  #####     - MFN
  #####     - DPP4
  #####     - TZD
  ################################################
  
  cprd <- cprd %>%
    mutate(SU = as.factor(SU),
           MFN = as.factor(MFN),
           DPP4 = as.factor(DPP4),
           TZD = as.factor(TZD))
  
  ################################################
  ##### Diabetes treatment
  #####   - Outcome month: hba1cmonth
  ################################################
  
  cprd <- cprd %>%
    mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
    mutate(hba1cmonth = as.numeric(hba1cmonth)) %>%
    select(-hba1cmonth_12, - hba1cmonth_6)
  
  
  ################################################
  ##### Biomarkers
  #####   - hba1c: prehba1c
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - BMI: prebmi
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - eGFR: preegfr
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Albumin:Creatine ratio: preacr
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Serum albumin: prealbumin_blood
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Alanine aminotransferase: prealt
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Aspartate aminotransferase: preast
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Bilirubin: prebilirubin
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Fasting glucose: prefastingglucose
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Fasting haematocrit: prehaematocrit
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Fasting haemoglobin: prehaemoglobin
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - High-density lipoprotein (HDL): prehdl
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Mean arterial BP: premap
  ################################################
  
  cprd <- cprd %>%
    mutate(premap = predbp + ((presbp - predbp) / 3))
  
  ################################################
  ##### Biomarkers
  #####   - Total cholesterol: pretotalcholesterol
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Biomarkers
  #####   - Triglycerides: pretriglyceride
  ################################################
  
  # Nothing to do 
  
  ################################################
  ##### Comorbidities
  #####   - Angina: predrug_earliest_angina
  ################################################
  
  # Nothing to do 
  
}


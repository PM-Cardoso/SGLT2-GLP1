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
  
  # This is wrong, needs to be changed. Katie's email on 21/10/22 at 11:15
  
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
  ##### Drug of interest
  #####   - GLP1 and SGLT2: drugclass
  ################################################
  
  cprd <- cprd %>%
    mutate(drugclass = factor(drugclass))
  
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
    mutate(SU = factor(SU, levels = c(0, 1), labels = c("No", "Yes")),
           MFN = factor(MFN, levels = c(0, 1), labels = c("No", "Yes")),
           DPP4 = factor(DPP4, levels = c(0, 1), labels = c("No", "Yes")),
           TZD = factor(TZD, levels = c(0, 1), labels = c("No", "Yes")))
  
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
  #####   - Angina: predrug_latest_angina
  ################################################
  
  # any angina medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_angina = ifelse(is.na(predrug_latest_angina), "No", "Yes")) %>%
    mutate(pre_angina = factor(pre_angina))
  
  # any angina medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_angina_5yrecent = difftime(dstartdate, predrug_latest_angina, units = "days") / 365.25) %>%
    mutate(pre_angina_5yrecent = ifelse(pre_angina_5yrecent > 5 | is.na(pre_angina_5yrecent), "No", "Yes")) %>%
    mutate(pre_angina_5yrecent = factor(pre_angina_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Chronic Liver Disease: predrug_latest_cld
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_cld = ifelse(is.na(predrug_latest_cld), "No", "Yes")) %>%
    mutate(pre_cld = factor(pre_cld))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_cld_5yrecent = difftime(dstartdate, predrug_latest_cld, units = "days") / 365.25) %>%
    mutate(pre_cld_5yrecent = ifelse(pre_cld_5yrecent > 5 | is.na(pre_cld_5yrecent), "No", "Yes")) %>%
    mutate(pre_cld_5yrecent = factor(pre_cld_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Diabetic Nephropathy: predrug_latest_diabeticnephropathy
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_diabeticnephropathy = ifelse(is.na(predrug_latest_diabeticnephropathy), "No", "Yes")) %>%
    mutate(pre_diabeticnephropathy = factor(pre_diabeticnephropathy))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_diabeticnephropathy_5yrecent = difftime(dstartdate, predrug_latest_diabeticnephropathy, units = "days") / 365.25) %>%
    mutate(pre_diabeticnephropathy_5yrecent = ifelse(pre_diabeticnephropathy_5yrecent > 5 | is.na(pre_diabeticnephropathy_5yrecent), "No", "Yes")) %>%
    mutate(pre_diabeticnephropathy_5yrecent = factor(pre_diabeticnephropathy_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Heart failure: predrug_latest_heartfailure
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_heartfailure = ifelse(is.na(predrug_latest_heartfailure), "No", "Yes")) %>%
    mutate(pre_heartfailure = factor(pre_heartfailure))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_heartfailure_5yrecent = difftime(dstartdate, predrug_latest_heartfailure, units = "days") / 365.25) %>%
    mutate(pre_heartfailure_5yrecent = ifelse(pre_heartfailure_5yrecent > 5 | is.na(pre_heartfailure_5yrecent), "No", "Yes")) %>%
    mutate(pre_heartfailure_5yrecent = factor(pre_heartfailure_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Hypertension: predrug_latest_hypertension
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_hypertension = ifelse(is.na(predrug_latest_heartfailure), "No", "Yes")) %>%
    mutate(pre_hypertension = factor(pre_hypertension))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_hypertension_5yrecent = difftime(dstartdate, predrug_latest_heartfailure, units = "days") / 365.25) %>%
    mutate(pre_hypertension_5yrecent = ifelse(pre_hypertension_5yrecent > 5 | is.na(pre_hypertension_5yrecent), "No", "Yes")) %>%
    mutate(pre_hypertension_5yrecent = factor(pre_hypertension_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Ischaemic Heart Disease: predrug_latest_ihd
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_ihd = ifelse(is.na(predrug_latest_ihd), "No", "Yes")) %>%
    mutate(pre_ihd = factor(pre_ihd))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_ihd_5yrecent = difftime(dstartdate, predrug_latest_ihd, units = "days") / 365.25) %>%
    mutate(pre_ihd_5yrecent = ifelse(pre_ihd_5yrecent > 5 | is.na(pre_ihd_5yrecent), "No", "Yes")) %>%
    mutate(pre_ihd_5yrecent = factor(pre_ihd_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Myocardial Infarction: predrug_latest_myocardialinfarction
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_myocardialinfarction = ifelse(is.na(predrug_latest_myocardialinfarction), "No", "Yes")) %>%
    mutate(pre_myocardialinfarction = factor(pre_myocardialinfarction))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_myocardialinfarction_5yrecent = difftime(dstartdate, predrug_latest_myocardialinfarction, units = "days") / 365.25) %>%
    mutate(pre_myocardialinfarction_5yrecent = ifelse(pre_myocardialinfarction_5yrecent > 5 | is.na(pre_myocardialinfarction_5yrecent), "No", "Yes")) %>%
    mutate(pre_myocardialinfarction_5yrecent = factor(pre_myocardialinfarction_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Neuropathy: predrug_latest_neuropathy
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_neuropathy = ifelse(is.na(predrug_latest_neuropathy), "No", "Yes")) %>%
    mutate(pre_neuropathy = factor(pre_neuropathy))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_neuropathy_5yrecent = difftime(dstartdate, predrug_latest_neuropathy, units = "days") / 365.25) %>%
    mutate(pre_neuropathy_5yrecent = ifelse(pre_neuropathy_5yrecent > 5 | is.na(pre_neuropathy_5yrecent), "No", "Yes")) %>%
    mutate(pre_neuropathy_5yrecent = factor(pre_neuropathy_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Peripheral Arterial Disease: predrug_latest_pad
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_pad = ifelse(is.na(predrug_latest_pad), "No", "Yes")) %>%
    mutate(pre_pad = factor(pre_pad))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_pad_5yrecent = difftime(dstartdate, predrug_latest_pad, units = "days") / 365.25) %>%
    mutate(pre_pad_5yrecent = ifelse(pre_pad_5yrecent > 5 | is.na(pre_pad_5yrecent), "No", "Yes")) %>%
    mutate(pre_pad_5yrecent = factor(pre_pad_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Retinopathy: predrug_latest_retinopathy
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_retinopathy = ifelse(is.na(predrug_latest_retinopathy), "No", "Yes")) %>%
    mutate(pre_retinopathy = factor(pre_retinopathy))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_retinopathy_5yrecent = difftime(dstartdate, predrug_latest_retinopathy, units = "days") / 365.25) %>%
    mutate(pre_retinopathy_5yrecent = ifelse(pre_retinopathy_5yrecent > 5 | is.na(pre_retinopathy_5yrecent), "No", "Yes")) %>%
    mutate(pre_retinopathy_5yrecent = factor(pre_retinopathy_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Cardiac Revascularisation: predrug_latest_revasc
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_revasc = ifelse(is.na(predrug_latest_revasc), "No", "Yes")) %>%
    mutate(pre_revasc = factor(pre_revasc))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_revasc_5yrecent = difftime(dstartdate, predrug_latest_revasc, units = "days") / 365.25) %>%
    mutate(pre_revasc_5yrecent = ifelse(pre_revasc_5yrecent > 5 | is.na(pre_revasc_5yrecent), "No", "Yes")) %>%
    mutate(pre_revasc_5yrecent = factor(pre_revasc_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Stroke: predrug_latest_stroke
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_stroke = ifelse(is.na(predrug_latest_stroke), "No", "Yes")) %>%
    mutate(pre_stroke = factor(pre_stroke))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_stroke_5yrecent = difftime(dstartdate, predrug_latest_stroke, units = "days") / 365.25) %>%
    mutate(pre_stroke_5yrecent = ifelse(pre_stroke_5yrecent > 5 | is.na(pre_stroke_5yrecent), "No", "Yes")) %>%
    mutate(pre_stroke_5yrecent = factor(pre_stroke_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Transient Ischaemic Attack: predrug_latest_tia
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_tia = ifelse(is.na(predrug_latest_tia), "No", "Yes")) %>%
    mutate(pre_tia = factor(pre_tia))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_tia_5yrecent = difftime(dstartdate, predrug_latest_tia, units = "days") / 365.25) %>%
    mutate(pre_tia_5yrecent = ifelse(pre_tia_5yrecent > 5 | is.na(pre_tia_5yrecent), "No", "Yes")) %>%
    mutate(pre_tia_5yrecent = factor(pre_tia_5yrecent))
  
  ################################################
  ##### Comorbidities
  #####   - Atrial fibrillation: predrug_latest_af
  ################################################
  
  # any Chronic Liver Disease medcode in history 
  
  cprd <- cprd %>%
    mutate(pre_af = ifelse(is.na(predrug_latest_af), "No", "Yes")) %>%
    mutate(pre_af = factor(pre_af))
  
  # any Chronic Liver Disease medcode in the last 5 years
  
  cprd <- cprd %>%
    mutate(pre_af_5yrecent = difftime(dstartdate, predrug_latest_af, units = "days") / 365.25) %>%
    mutate(pre_af_5yrecent = ifelse(pre_af_5yrecent > 5 | is.na(pre_af_5yrecent), "No", "Yes")) %>%
    mutate(pre_af_5yrecent = factor(pre_af_5yrecent))
  
  
  
  ###############################################################################
  ###############################################################################
  #################### Final dataset - all patients #############################
  ###############################################################################
  ###############################################################################
  
  final.dataset <- cprd %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1c_final,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, Ethnicity, Deprivation, Smoke, 
      # Diabetes treatment 
      drugline, SU, MFN, DPP4, TZD, hba1cmonth,
      # Biomarkers
      
    )
  
  
}


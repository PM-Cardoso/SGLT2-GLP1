####################
## Description:
##  - This file includes the framework for patient selection.
####################

# load libraries
require(tidyverse)


###############################################################################
###############################################################################
############################### Functions #####################################
###############################################################################
###############################################################################


set_up_data_sglt2_glp1 <- function(dataset.type) {
  ##### Input variables
  # dataset.type: a character string mentioning the type of dataset required
  
  # initial checks
  if (missing(dataset.type)) {stop("'dataset.type' needs to be supplied")}
  if (!is.character(dataset.type)) {stop("'dataset.type' must be a character string")}
  if (!(dataset.type %in% c("full_cohort"))) {stop("'dataset.type' must one of: full_cohort /")}
  
  
  # load original dataset # name - t2d_1stinstance_2ndline
  load("/slade/CPRD_data/mastermind_2022/20221019_t2d_1stinstance_2ndline.Rda")
  
  cprd <- t2d_1stinstance_2ndline
  
  
  ################################################
  ##### Select only SGLT-2 and GLP-1
  ################################################
  
  #SGLT2 vs GLP1
  cprd <- cprd %>% 
    filter(drugclass == "GLP1" | drugclass == "SGLT2") # nrow = 160744
  
  # table(cprd$drugclass)
  
  
  ################################################
  ##### Drop if at least 6 months of follow-up available (treatment start < 01/04/2020) 
  ################################################
  
  cprd <- cprd %>%
    mutate(follow_up = ifelse(dstartdate > "2020-04-01", 1, NA_real_))
  
  # table(cprd$follow_up)
  
  
  cprd <- cprd %>%
    filter(is.na(follow_up))
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  cprd <- cprd %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # table(cprd$timeprevcombo_less61)
  
  
  cprd <- cprd %>% 
    filter(is.na(timeprevcombo_less61))

  
  ################################################
  ##### Drop patients initiating before 1/1/2015
  ################################################
  
  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < "2015-01-01", 1 , NA_real_))
  
  # table(cprd$dstartdate_cutoff)
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))

  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # table(cprd$INS) 
  
  cprd <- cprd %>% 
    filter(INS == 0)      
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # table(cprd$drugline)
  
  cprd <- cprd %>%
    filter(drugline > 1)
  
  ################################################
  ##### Drop if within 1 year of diagnosis
  ################################################
  
  cprd <- cprd %>%
    mutate(dstart_1year = difftime(dstartdate, dm_diag_date, units = "days")) %>%
    mutate(dstart_1year = ifelse(dstart_1year < 365.25, 1, NA_real_))
  
  # table(cprd$dstart_1year)
  
  cprd <- cprd %>%
    filter(is.na(dstart_1year))
  
  
  ################################################
  ##### Drop patients initiating drug whilst on the opposite
  ################################################
  
  # table(cprd$drugclass, cprd$SGLT2)
  
  # table(cprd$drugclass, cprd$GLP1)
  
  
  cprd <- cprd %>%
    mutate(opposite_drug = ifelse( (drugclass == "GLP1" & SGLT2 == 1) | (drugclass == "SGLT2" & GLP1 == 1) , 1, NA_real_)) %>%
    filter(is.na(opposite_drug))
  
  
# ------------------------------------------------------------------------------------------------
  
  # This needs more consideration
  
# ------------------------------------------------------------------------------------------------
  
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # table(cprd$preckdstage)
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if HbA1c not in 53 - 120 at baseline + missing initial hba1c
  ################################################
  
  cprd <- cprd %>% 
    mutate(hb_extreme= ifelse(prehba1c < 53 | prehba1c > 120 | is.na(prehba1c), 1, 0))
  
  # table(cprd$hb_extreme)
  
  cprd <- cprd %>% 
    filter(hb_extreme==0)
  
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  
  ################################################
  ##### Drug of interest
  ################################################
  
  cprd <- cprd %>%
    #####   - GLP1 and SGLT2: drugclass
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
  ################################################
  
  cprd <- cprd %>%
    #####   - Age: agetx (new var)
    mutate(agetx = as.numeric(dstartdate_age)) %>%
    #####   - Sex: malesex
    mutate(malesex = ifelse(gender == 1, 1, 0)) %>%
    mutate(malesex = as.factor(malesex)) %>%
    #####   - Duration of diabetes: t2dmduration
    mutate(t2dmduration = as.numeric(dstartdate_dm_dur)) %>%
    #####   - Ethnicity: ethnicity
    mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed"))) %>%
    #####   - Deprivation: deprivation
    mutate(deprivation = factor(imd2015_10)) %>%
    #####   - Smoking Status: smoke
    mutate(smoke = factor(smoking_cat)) %>%
    #####   - Line Therapy: drugline: turn all > 4 to 5+
    mutate(drugline = ifelse(drugline > 4, 5, drugline)) %>%
    mutate(drugline = as.factor(drugline))
  
  ################################################
  ##### Diabetes treatment
  ################################################
  
  cprd <- cprd %>%
    #####   - Drugs taken alongside treatment
    #####     - SU
    #####     - MFN
    #####     - DPP4
    #####     - TZD
    mutate(SU = factor(SU, levels = c(0, 1), labels = c("No", "Yes")),
           MFN = factor(MFN, levels = c(0, 1), labels = c("No", "Yes")),
           DPP4 = factor(DPP4, levels = c(0, 1), labels = c("No", "Yes")),
           TZD = factor(TZD, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  #####   - Outcome month: hba1cmonth
    mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
    mutate(hba1cmonth = as.numeric(hba1cmonth))
  
  
  ################################################
  ##### Biomarkers
  ################################################
  
    #####   - hba1c: prehba1c (Nothing to do)
    #####   - BMI: prebmi (Nothing to do)
    #####   - eGFR: preegfr (Nothing to do)
    #####   - Albumin:Creatine ratio: preacr (Nothing to do)
    #####   - Serum albumin: prealbumin_blood (Nothing to do)
    #####   - Alanine aminotransferase: prealt (Nothing to do)
    #####   - Aspartate aminotransferase: preast (Nothing to do)
    #####   - Bilirubin: prebilirubin (Nothing to do)
    #####   - Fasting glucose: prefastingglucose (Nothing to do)
    #####   - Fasting haematocrit: prehaematocrit (Nothing to do)
    #####   - Fasting haemoglobin: prehaemoglobin (Nothing to do)
    #####   - High-density lipoprotein (HDL): prehdl (Nothing to do)
  
    #####   - Mean arterial BP: premap
  
  cprd <- cprd %>%
    mutate(premap = predbp + ((presbp - predbp) / 3))
  
    #####   - Total cholesterol: pretotalcholesterol (Nothing to do)
    #####   - Triglycerides: pretriglyceride (Nothing to do)
  
  
  
  ################################################
  ##### Comorbidities
  ################################################
  
  cprd <- cprd %>%
    #####   - Angina: predrug_earliest_angina
    mutate(pre_angina = ifelse(is.na(predrug_earliest_angina) | predrug_earliest_angina > dstartdate, "No", "Yes")) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(pre_cld = ifelse(is.na(predrug_earliest_cld) | predrug_earliest_cld > dstartdate, "No", "Yes")) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(pre_diabeticnephropathy = ifelse(is.na(predrug_earliest_diabeticnephropathy) | predrug_earliest_diabeticnephropathy > dstartdate, "No", "Yes")) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(pre_heartfailure = ifelse(is.na(predrug_earliest_heartfailure) | predrug_earliest_heartfailure > dstartdate, "No", "Yes")) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(pre_hypertension = ifelse(is.na(predrug_earliest_hypertension) | predrug_earliest_hypertension > dstartdate, "No", "Yes")) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(pre_ihd = ifelse(is.na(predrug_earliest_ihd) | predrug_earliest_ihd > dstartdate, "No", "Yes")) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(pre_myocardialinfarction = ifelse(is.na(predrug_earliest_myocardialinfarction) | predrug_earliest_myocardialinfarction > dstartdate, "No", "Yes")) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(pre_neuropathy = ifelse(is.na(predrug_earliest_neuropathy) | predrug_earliest_neuropathy > dstartdate, "No", "Yes")) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(pre_pad = ifelse(is.na(predrug_earliest_pad) | predrug_earliest_pad > dstartdate, "No", "Yes")) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(pre_retinopathy = ifelse(is.na(predrug_earliest_retinopathy) | predrug_earliest_retinopathy > dstartdate, "No", "Yes")) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(pre_revasc = ifelse(is.na(predrug_earliest_revasc) | predrug_earliest_revasc > dstartdate, "No", "Yes")) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(pre_stroke = ifelse(is.na(predrug_earliest_stroke) | predrug_earliest_stroke > dstartdate, "No", "Yes")) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pre_tia = ifelse(is.na(predrug_earliest_tia) | predrug_earliest_tia > dstartdate, "No", "Yes")) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(pre_af = ifelse(is.na(predrug_earliest_af) | predrug_earliest_af > dstartdate, "No", "Yes"))
  
  
  
  ###############################################################################
  ###############################################################################
  #################### Final dataset - all patients #############################
  ###############################################################################
  ###############################################################################
  #
  # Add all variables necessary for ALL analysis in the paper.
  #
  
  final.dataset <- cprd %>%
    select(
      # information regarding patient
      patid, multi_drug_start,
      # response hba1c
      posthba1c_final,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, ethnicity, deprivation, smoke, 
      # Diabetes treatment 
      drugline, SU, MFN, DPP4, TZD, hba1cmonth, dstartdate, dstopdate,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbumin_blood, prealt, preast, prebilirubin, prefastingglucose,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      pre_angina, pre_cld, pre_diabeticnephropathy, pre_heartfailure, pre_hypertension, pre_ihd, pre_myocardialinfarction, 
      pre_neuropathy, pre_pad, pre_retinopathy, pre_revasc, pre_stroke, pre_tia, pre_af
    )
  
  
  if (dataset.type == "full_cohort") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  ###############################################################################
  
  hba1c.model.dataset <- final.dataset

  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(multi_drug_start == 0)
    
  
  ################################################
  ##### Drop patients with only 1 script (eg: dstartdate = dstopdate)
  ################################################
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    mutate(one_script = ifelse(dstartdate == dstopdate, 1, NA_real_))
  
  
  print("you went too long")
  
}


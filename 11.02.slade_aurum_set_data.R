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
  if (!(dataset.type %in% c("full.cohort", "hba1c.train", "hba1c.test", "weight.dataset"))) {
    stop("'dataset.type' must one of: full.cohort / hba1c.train / hba1c.test / weight.dataset")
  }
  
  
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
  ##### Drop patients diagnosed within 3 months of registration
  ################################################
  
  # regstartdate # date of registration
  # 
  # dm_diag_date # date of diagnosis
  
  # cprd <- cprd %>%
  #   mutate(dm_reg_3_months = ifelse(difftime(regstartdate, dm_diag_date, units = "days") > 0 & difftime(regstartdate, dm_diag_date, units = "days") < 90, 1, NA_real_))
  
  
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
  ##### Drop patients initiating drug whilst on the opposite
  ################################################
  
  # table(cprd$drugclass, cprd$SGLT2)
  
  # table(cprd$drugclass, cprd$GLP1)
  
  
  cprd <- cprd %>%
    mutate(opposite_drug = ifelse( (drugclass == "GLP1" & SGLT2 == 1) | (drugclass == "SGLT2" & GLP1 == 1) , 1, NA_real_)) %>%
    filter(is.na(opposite_drug))
  
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # table(cprd$preckdstage)
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  ### Add variable that identifies an individual entry in the data
  
  cprd <- cprd %>%
    mutate(pated = paste(patid, drugclass, dstartdate, sep = ".")) %>%
  
  ################################################
  ##### Drug of interest
  ################################################
  
    #####   - GLP1 and SGLT2: drugclass
    mutate(drugclass = factor(drugclass))
  
  ################################################
  ##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
  #####   - posthba1c12m but if missing
  #####     - posthba1c6m
  ################################################
  
  cprd <- cprd %>%
    mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%
    mutate(posthba1cfinal = as.numeric(posthba1cfinal))
    
  
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
    #####   - Serum albumin: prealbumin_blood (remove the _ from the name)
  
  cprd <- cprd %>%
    rename("prealbuminblood" = "prealbumin_blood",
           "prealbuminblooddate" = "prealbumin_blooddate",
           "prealbuminblooddrugdiff" = "prealbumin_blooddrugdiff")
  
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
    mutate(preangina = ifelse(is.na(predrug_earliest_angina) | predrug_earliest_angina > dstartdate, "No", "Yes")) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = ifelse(is.na(predrug_earliest_cld) | predrug_earliest_cld > dstartdate, "No", "Yes")) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = ifelse(is.na(predrug_earliest_diabeticnephropathy) | predrug_earliest_diabeticnephropathy > dstartdate, "No", "Yes")) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = ifelse(is.na(predrug_earliest_heartfailure) | predrug_earliest_heartfailure > dstartdate, "No", "Yes")) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = ifelse(is.na(predrug_earliest_hypertension) | predrug_earliest_hypertension > dstartdate, "No", "Yes")) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = ifelse(is.na(predrug_earliest_ihd) | predrug_earliest_ihd > dstartdate, "No", "Yes")) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = ifelse(is.na(predrug_earliest_myocardialinfarction) | predrug_earliest_myocardialinfarction > dstartdate, "No", "Yes")) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = ifelse(is.na(predrug_earliest_neuropathy) | predrug_earliest_neuropathy > dstartdate, "No", "Yes")) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = ifelse(is.na(predrug_earliest_pad) | predrug_earliest_pad > dstartdate, "No", "Yes")) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = ifelse(is.na(predrug_earliest_retinopathy) | predrug_earliest_retinopathy > dstartdate, "No", "Yes")) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = ifelse(is.na(predrug_earliest_revasc) | predrug_earliest_revasc > dstartdate, "No", "Yes")) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = ifelse(is.na(predrug_earliest_stroke) | predrug_earliest_stroke > dstartdate, "No", "Yes")) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = ifelse(is.na(predrug_earliest_tia) | predrug_earliest_tia > dstartdate, "No", "Yes")) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = ifelse(is.na(predrug_earliest_af) | predrug_earliest_af > dstartdate, "No", "Yes"))
  
  
  
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
      patid, pated, multi_drug_start, timeprevcombo,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, ethnicity, deprivation, smoke, 
      # Diabetes treatment 
      drugline, SU, MFN, DPP4, TZD, hba1cmonth, dstartdate, dstopdate,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate
    )
  
  
  # if full cohort was requested
  if (dataset.type == "full.cohort") {
    return(final.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  hba1c.model.dataset <- final.dataset
  
  ################################################
  ##### Drop if first-line treatment
  ################################################

  # table(hba1c.model.dataset$drugline)

  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(drugline != 1)
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(hba1c.model.dataset$multi_drug_start)
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################

  hba1c.model.dataset <- hba1c.model.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))

  # table(hba1c.model.dataset$timeprevcombo_less61)

  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################

  hba1c.model.dataset <- hba1c.model.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # table(hba1c.model.dataset$hb_extreme_53)
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(is.na(hb_extreme_53))
  
  ################################################
  ##### Drop if HbA1c >120
  ################################################
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    mutate(hb_extreme_120 = ifelse(prehba1c > 120, 1, NA_real_))
  
  # table(hba1c.model.dataset$hb_extreme_120)
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(is.na(hb_extreme_120))
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # summary(hba1c.model.dataset$prehba1c)
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # summary(hba1c.model.dataset$posthba1cfinal)
  
  hba1c.model.dataset <- hba1c.model.dataset %>%
    filter(!is.na(posthba1cfinal))
    
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.hba1c.model.dataset <- hba1c.model.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, ethnicity, deprivation, smoke, 
      # Diabetes treatment 
      drugline, SU, MFN, DPP4, TZD, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  
  ################################################
  ##### Define training and testing cohorts
  ################################################
  
  # nrow(final.hba1c.model.dataset); table(final.hba1c.model.dataset$drugclass)
  
  # Training dataset
  set.seed(123)
  hba1c.model.dataset.train <- final.hba1c.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6)
  
  if (dataset.type == "hba1c.train") {
    return(hba1c.model.dataset.train)
  }
  
  
  # Testing dataset
  hba1c.model.dataset.test <- subset(final.hba1c.model.dataset, !(pated %in% hba1c.model.dataset.train$pated))
  
  if (dataset.type == "hba1c.test") {
    return(hba1c.model.dataset.test)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Weight population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  weight.dataset <- final.dataset
  
  ################################################
  ##### Drop if only 1 script (eg. dstartdate = dstopdate)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(one_script = ifelse(dstartdate == dstopdate, 1, NA_real_))
  
  # table(weight.dataset$one_script)
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(one_script))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(weight.dataset$multi_drug_start)
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # summary(weight.dataset$preweight)
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # summary(weight.dataset$postweight)
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(postweight))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.weight.dataset <- weight.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, ethnicity, deprivation, smoke, 
      # Diabetes treatment 
      drugline, SU, MFN, DPP4, TZD, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight
    )
  
  
  if (dataset.type == "weight.dataset") {
    return(final.weight.dataset)
  }
  
  
}


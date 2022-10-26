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
  if (!(dataset.type %in% c("full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset"))) {
    stop("'dataset.type' must one of: full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset")
  }
  
  
  # # load original dataset # name - t2d_1stinstance_2ndline
  # load("/slade/CPRD_data/mastermind_2022/20221019_t2d_1stinstance_2ndline.Rda")
  # 
  # cprd <- t2d_1stinstance_2ndline
  
  
  # load original dataset # name - t2d_1stinstance
  load("/slade/CPRD_data/mastermind_2022/20221024_t2d_1stinstance.Rda")
  
  cprd <- t2d_1stinstance
  
  ################################################
  ##### Select only SGLT-2 and GLP-1
  ################################################
  
  #SGLT2 vs GLP1
  cprd <- cprd %>% 
    filter(drugclass == "GLP1" | drugclass == "SGLT2")
  
  # table(cprd$drugclass)
  
  
  
  ################################################
  ##### Drop patients initiating before 1/1/2013
  ################################################
  
  #######################
  # Explore adjusted HbA1c repsonse by calendar year
  
  cprd  <- cprd %>%
    mutate(yrdrugstart = format(dstartdate, format = "%Y")) %>%
    mutate(yrdrugstart = as.numeric(yrdrugstart))
  
  # 
  # # table(cprd$drugclass, cprd$yrdrugstart)
  # #        2005  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017  2018  2019  2020
  # # GLP1      1   231  1711  3355  4810  4674  6074  4773  4359  4851  4674  5221  5845  7634  4822
  # # SGLT2     0     0     0     0     0     0     1  1624  7083 13057 14286 15707 17356 19328 11408
  # 
  # 
  # cprd <- cprd %>%
  #   mutate(adj_posthba1c = ifelse(is.na(posthba1c12m), posthba1c6m - prehba1c, posthba1c12m - prehba1c))
  # 
  # cprd <- cprd %>%
  #   mutate(adj_posthba1c = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m),
  #          ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN) %>%
  #   mutate(ncurrtx = as.factor(ncurrtx))
  # 
  # cprd <- cprd %>%
  #   filter(!str_detect(drugsubstances, "&"))
  # 
  # # change drugclass to drugsubstances
  # 
  # m1 <- ols(adj_posthba1c~prehba1c+rcs(yrdrugstart,3)*drugclass + drugline + ncurrtx,data=cprd,x=TRUE,y=TRUE)
  # 
  # m1.p <- Predict(m1,prehba1c=median(cprd$prehba1c, na.rm = TRUE),drugclass = c("SGLT2", "GLP1"), yrdrugstart=c(seq(2007,2020)), drugline=c(2), ncurrtx = c(1)) %>%
  #   as.data.frame() %>%
  #   mutate(Response = yhat-prehba1c,
  #          ci.l = lower - prehba1c,
  #          ci.u = upper - prehba1c) %>%
  #   mutate(Response = ifelse(drugclass == "SGLT2" & yrdrugstart < 2013, NA, Response),
  #          ci.l = ifelse(drugclass == "SGLT2" & yrdrugstart < 2013, NA, Response),
  #          ci.u = ifelse(drugclass == "SGLT2" & yrdrugstart < 2013, NA, Response))
  # 
  # m1.p %>%
  #   ggplot() + theme_bw() +
  #   geom_line(aes(x = yrdrugstart, y = Response, colour = drugclass)) +
  #   xlab("Year of drug start") + ylab("HbA1c change")
    
  #######################
  
  
  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < "2013-01-01", 1 , NA_real_))
  
  # table(cprd$dstartdate_cutoff); table(cprd$dstartdate_cutoff, cprd$drugclass)
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))

  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # table(cprd$INS); table(cprd$INS, cprd$drugclass)
  
  cprd <- cprd %>% 
    filter(INS == 0)      
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # table(cprd$preckdstage); table(cprd$preckdstage, cprd$drugclass)
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # table(cprd$drugline); table(cprd$drugline, cprd$drugclass)
  
  cprd <- cprd %>%
    filter(drugline != 1)
  
  
  
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
    mutate(drugline = ifelse(drugline > 4, "5+", drugline)) %>%
    mutate(drugline = as.factor(drugline))
  
  ################################################
  ##### Diabetes treatment
  ################################################
  
  cprd <- cprd %>%
    #####   - Drugs taken alongside treatment
    # #####     - SU
    # #####     - MFN
    # #####     - DPP4
    # #####     - TZD
    # mutate(SU = factor(SU, levels = c(0, 1), labels = c("No", "Yes")),
    #        MFN = factor(MFN, levels = c(0, 1), labels = c("No", "Yes")),
    #        DPP4 = factor(DPP4, levels = c(0, 1), labels = c("No", "Yes")),
    #        TZD = factor(TZD, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    mutate(ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN) %>%
    mutate(ncurrtx = ifelse(ncurrtx > 4, "5+", ncurrtx)) %>%
    mutate(ncurrtx = as.factor(ncurrtx)) %>%
    
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
    mutate(preangina = factor(preangina)) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = ifelse(is.na(predrug_earliest_cld) | predrug_earliest_cld > dstartdate, "No", "Yes")) %>%
    mutate(precld = factor(precld)) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = ifelse(is.na(predrug_earliest_diabeticnephropathy) | predrug_earliest_diabeticnephropathy > dstartdate, "No", "Yes")) %>%
    mutate(prediabeticnephropathy = factor(prediabeticnephropathy)) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = ifelse(is.na(predrug_earliest_heartfailure) | predrug_earliest_heartfailure > dstartdate, "No", "Yes")) %>%
    mutate(preheartfailure = factor(preheartfailure)) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = ifelse(is.na(predrug_earliest_hypertension) | predrug_earliest_hypertension > dstartdate, "No", "Yes")) %>%
    mutate(prehypertension = factor(prehypertension)) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = ifelse(is.na(predrug_earliest_ihd) | predrug_earliest_ihd > dstartdate, "No", "Yes")) %>%
    mutate(preihd = factor(preihd)) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = ifelse(is.na(predrug_earliest_myocardialinfarction) | predrug_earliest_myocardialinfarction > dstartdate, "No", "Yes")) %>%
    mutate(premyocardialinfarction = factor(premyocardialinfarction)) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = ifelse(is.na(predrug_earliest_neuropathy) | predrug_earliest_neuropathy > dstartdate, "No", "Yes")) %>%
    mutate(preneuropathy = factor(preneuropathy)) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = ifelse(is.na(predrug_earliest_pad) | predrug_earliest_pad > dstartdate, "No", "Yes")) %>%
    mutate(prepad = factor(prepad)) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = ifelse(is.na(predrug_earliest_retinopathy) | predrug_earliest_retinopathy > dstartdate, "No", "Yes")) %>%
    mutate(preretinopathy = factor(preretinopathy)) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = ifelse(is.na(predrug_earliest_revasc) | predrug_earliest_revasc > dstartdate, "No", "Yes")) %>%
    mutate(prerevasc = factor(prerevasc)) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = ifelse(is.na(predrug_earliest_stroke) | predrug_earliest_stroke > dstartdate, "No", "Yes")) %>%
    mutate(prestroke = factor(prestroke)) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = ifelse(is.na(predrug_earliest_tia) | predrug_earliest_tia > dstartdate, "No", "Yes")) %>%
    mutate(pretia = factor(pretia)) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = ifelse(is.na(predrug_earliest_af) | predrug_earliest_af > dstartdate, "No", "Yes")) %>%
    mutate(preaf = factor(preaf))
  
  
  
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
      drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
      # eGFR analysis
      postegfr12m, postegfr6m
    )
  
  # nrow(final.dataset); table(final.dataset$drugclass)
  
  # if full cohort was requested
  if (dataset.type == "full.cohort") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Propensity score model ###############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  ps.model.dataset <- final.dataset
  
  #:----------------------------------------------------
  # Select variables needed
  
  ps.model.dataset <- ps.model.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, ethnicity, deprivation, smoke, 
      # Diabetes treatment 
      drugline, ncurrtx, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, 
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # nrow(ps.model.dataset); table(ps.model.dataset$drugclass)
  
  # Training dataset
  set.seed(123)
  ps.model.dataset.train <- ps.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6) %>%
    ungroup()
  
  if (dataset.type == "ps.model.train") {
    return(ps.model.dataset.train)
  }
  
  
  # Testing dataset
  ps.model.dataset.test <- subset(ps.model.dataset, !(pated %in% ps.model.dataset.train$pated))
  
  if (dataset.type == "ps.model.test") {
    return(ps.model.dataset.test)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  hba1c.model.dataset.train <- ps.model.dataset.train %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  hba1c.model.dataset.test <- ps.model.dataset.test %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(hba1c.model.dataset.train$multi_drug_start); table(hba1c.model.dataset.train$multi_drug_start, hba1c.model.dataset.train$drugclass) 
  # 
  # table(hba1c.model.dataset.test$multi_drug_start); table(hba1c.model.dataset.test$multi_drug_start, hba1c.model.dataset.test$drugclass)
  
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(multi_drug_start == 0)
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  
  # table(hba1c.model.dataset.train$timeprevcombo_less61); table(hba1c.model.dataset.train$timeprevcombo_less61, hba1c.model.dataset.train$drugclass) 
  # 
  # table(hba1c.model.dataset.test$timeprevcombo_less61); table(hba1c.model.dataset.test$timeprevcombo_less61, hba1c.model.dataset.test$drugclass) 
  
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(timeprevcombo_less61))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(timeprevcombo_less61))
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # table(hba1c.model.dataset.train$hb_extreme_53); table(hba1c.model.dataset.train$hb_extreme_53, hba1c.model.dataset.train$drugclass)
  # 
  # table(hba1c.model.dataset.test$hb_extreme_53); table(hba1c.model.dataset.test$hb_extreme_53, hba1c.model.dataset.test$drugclass)
  
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(hb_extreme_53))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # table(is.na(hba1c.model.dataset.train$prehba1c)); table(is.na(hba1c.model.dataset.train$prehba1c), hba1c.model.dataset.train$drugclass)
  #
  # table(is.na(hba1c.model.dataset.test$prehba1c)); table(is.na(hba1c.model.dataset.test$prehba1c), hba1c.model.dataset.test$drugclass)

  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(prehba1c))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # table(is.na(hba1c.model.dataset.train$posthba1cfinal)); table(is.na(hba1c.model.dataset.train$posthba1cfinal), hba1c.model.dataset.train$drugclass)
  # 
  # table(is.na(hba1c.model.dataset.test$posthba1cfinal)); table(is.na(hba1c.model.dataset.test$posthba1cfinal), hba1c.model.dataset.test$drugclass)
  
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(posthba1cfinal))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(posthba1cfinal))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  # Training dataset
  final.hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # nrow(final.hba1c.model.dataset.train); table(final.hba1c.model.dataset.train$drugclass)
  
  
  if (dataset.type == "hba1c.train") {
    return(final.hba1c.model.dataset.train)
  }
  
  # Testing dataset
  final.hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # nrow(final.hba1c.model.dataset.test); table(final.hba1c.model.dataset.test$drugclass)
  
  if (dataset.type == "hba1c.test") {
    return(final.hba1c.model.dataset.test)
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
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(weight.dataset$multi_drug_start); table(weight.dataset$multi_drug_start, weight.dataset$drugclass)
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # table(weight.dataset$timeprevcombo_less61); table(weight.dataset$timeprevcombo_less61, weight.dataset$drugclass)
  
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # table(weight.dataset$hb_extreme_53); table(weight.dataset$hb_extreme_53, weight.dataset$drugclass)
  
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # table(is.na(weight.dataset$prehba1c)); table(is.na(weight.dataset$prehba1c), weight.dataset$drugclass)
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # table(is.na(weight.dataset$preweight)); table(is.na(weight.dataset$preweight), weight.dataset$drugclass)
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # table(is.na(weight.dataset$postweight)); table(is.na(weight.dataset$postweight), weight.dataset$drugclass)
  
  
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
      agetx, malesex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight
    )

   # nrow(final.weight.dataset); table(final.weight.dataset$drugclass)
  
  if (dataset.type == "weight.dataset") {
    return(final.weight.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Discontinuation population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  discontinuation.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(discontinuation.dataset$multi_drug_start); table(discontinuation.dataset$multi_drug_start, discontinuation.dataset$drugclass)
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # table(discontinuation.dataset$timeprevcombo_less61); table(discontinuation.dataset$timeprevcombo_less61, discontinuation.dataset$drugclass)
  
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # table(discontinuation.dataset$hb_extreme_53); table(discontinuation.dataset$hb_extreme_53, discontinuation.dataset$drugclass)
  
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # table(is.na(discontinuation.dataset$prehba1c)); table(is.na(discontinuation.dataset$prehba1c), discontinuation.dataset$drugclass)
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.discontinuation.dataset <- discontinuation.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # nrow(final.discontinuation.dataset); table(final.discontinuation.dataset$drugclass)
  
  if (dataset.type == "discontinuation.dataset") {
    return(final.discontinuation.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################## eGFR population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  egfr.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # table(egfr.dataset$multi_drug_start); table(egfr.dataset$multi_drug_start, egfr.dataset$drugclass)
  
  egfr.dataset <- egfr.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # table(egfr.dataset$timeprevcombo_less61); table(egfr.dataset$timeprevcombo_less61, egfr.dataset$drugclass)
  
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # table(egfr.dataset$hb_extreme_53); table(egfr.dataset$hb_extreme_53, egfr.dataset$drugclass)
  
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # table(is.na(egfr.dataset$prehba1c)); table(is.na(egfr.dataset$prehba1c), egfr.dataset$drugclass)
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if baseline eGRF is missing
  ################################################
  
  # table(is.na(egfr.dataset$preegfr)); table(is.na(egfr.dataset$preegfr), egfr.dataset$drugclass)
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(preegfr))
  
  
  ################################################
  ##### Drop if post eGRF is missing
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(postegfr = ifelse(is.na(postegfr12m), postegfr6m, postegfr12m))
  
  # table(is.na(egfr.dataset$postegfr)); table(is.na(egfr.dataset$postegfr), egfr.dataset$drugclass)
  
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(postegfr))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.egfr.dataset <- egfr.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, malesex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # eGFR analysis
      postegfr
    )
  
  # nrow(final.egfr.dataset); table(final.egfr.dataset$drugclass)
  
  if (dataset.type == "egfr.dataset") {
    return(final.egfr.dataset)
  }
  
  
}


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
  if (!(dataset.type %in% c("diagnostics", "synthetic", "full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset"))) {
    stop("'dataset.type' must one of: diagnostics / synthetic / full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset")
  }
  
  
  # load original dataset # name - t2d_1stinstance
  load("/slade/CPRD_data/mastermind_2022/20221110_t2d_1stinstance.Rda")
  
  cprd <- t2d_1stinstance
  
  
  ################################################
  ##### Select only SGLT-2 and GLP-1
  ################################################
  
  #SGLT2 vs GLP1
  cprd <- cprd %>% 
    filter(drugclass == "GLP1" | drugclass == "SGLT2")
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Select only SGLT-2 and GLP-1")
    print("################################################")
    print(nrow(cprd))
    print(table(cprd$drugclass))
    
  }
  
  
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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients initiating before 1/1/2013")
    print("################################################")
    print(table(cprd$dstartdate_cutoff))
    print(table(cprd$dstartdate_cutoff, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))

  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if treated with insulin when starting new drug")
    print("################################################")
    print(table(cprd$INS))
    print(table(cprd$INS, cprd$drugclass))
    
  }
  
  cprd <- cprd %>% 
    filter(INS == 0)      
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients with ESRD")
    print("################################################")
    print(table(cprd$preckdstage))
    print(table(cprd$preckdstage, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if first-line treatment")
    print("################################################")
    print(table(cprd$drugline_all))
    print(table(cprd$drugline_all, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(drugline_all != 1)
  
  
  ################################################
  ##### Drop if semaglutide
  ################################################
  
  cprd <- cprd %>%
    mutate(semaglutide_drug = ifelse(str_detect(drugsubstances, "Semaglutide"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if semaglutide")
    print("################################################")
    print(table(cprd$semaglutide_drug))
    
  }
  
  cprd <- cprd %>%
    filter(is.na(semaglutide_drug))
  
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
  # 1 - GLP1; 2 - SGLT2
  
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
    #####   - Sex: sex
    mutate(sex = factor(ifelse(gender == 1, "Male", "Female"))) %>%
    
    #####   - Duration of diabetes: t2dmduration
    mutate(t2dmduration = as.numeric(dstartdate_dm_dur_all)) %>%
    #####   - Ethnicity: ethnicity
    mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed"))) %>%
    
    #####   - Deprivation: deprivation
    mutate(deprivation = factor(imd2015_10)) %>%
    
    #####   - Smoking Status: smoke
    mutate(smoke = factor(smoking_cat)) %>%
    
    #####   - Line Therapy: drugline: turn all > 4 to 5+
    mutate(drugline = ifelse(drugline_all > 4, 5, drugline_all)) %>%
    mutate(drugline = factor(drugline, levels = c(2, 3, 4, 5), labels = c("2", "3", "4", "5+"))) %>%
    
    #####   - Hospitalisations in previous year
    mutate(prehospitalisation = factor(hosp_admission_prev_year, levels = c(0, 1), labels = c("No", "Yes")))
  
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
    mutate(ncurrtx = ifelse(ncurrtx > 4, 5, ncurrtx)) %>%
    mutate(ncurrtx = factor(ncurrtx, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+"))) %>%
    
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
    mutate(preangina = factor(predrug_angina, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = factor(predrug_cld, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = factor(predrug_diabeticnephropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = factor(predrug_hypertension, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = factor(predrug_myocardialinfarction, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = factor(predrug_revasc, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = factor(predrug_stroke, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = factor(predrug_tia, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = factor(predrug_af, levels = c(0, 1), labels = c("No", "Yes")))
  
  
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
      patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
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
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Final dataset - all patients")
    print("################################################")
    print(nrow(final.dataset))
    print(table(final.dataset$drugclass))
    
  }
  
  # if full cohort was requested
  if (dataset.type == "full.cohort") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Synthetic dataset ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # Create synthetic dataset
  if (dataset.type == "synthetic") {
    # load package
    require(synthpop)
    
    set.seed(123)
    syn.dataset <- synthpop::syn(final.dataset %>%
                                   select(
                                     # response hba1c
                                     posthba1cfinal,
                                     # therapies of interest
                                     drugclass,
                                     # Sociodemographic features
                                     agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
                                     # Diabetes treatment 
                                     drugline, ncurrtx, hba1cmonth, yrdrugstart,
                                     # Biomarkers
                                     prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
                                     prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
                                     # Comorbidities
                                     preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
                                     preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
                                   ) %>%
                                   sample_n(520),
                                 print.flag = FALSE
    )
    
    return(syn.dataset$syn)
    
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
      agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
      # Diabetes treatment 
      drugline, ncurrtx, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, 
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # Training dataset
  set.seed(123)
  ps.model.dataset.train <- ps.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6) %>%
    ungroup() %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score training cohort")
    print("################################################")
    print(nrow(ps.model.dataset.train))
    print(table(ps.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "ps.model.train") {
    return(ps.model.dataset.train)
  }
  
  
  # Testing dataset
  ps.model.dataset.test <- subset(ps.model.dataset, !(pated %in% ps.model.dataset.train$pated)) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score testing cohort")
    print("################################################")
    print(nrow(ps.model.dataset.test))
    print(table(ps.model.dataset.test$drugclass))
    
  }
  
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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model")
    print("################################################")
    
  }
  
  hba1c.model.dataset.train <- ps.model.dataset.train %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  hba1c.model.dataset.test <- ps.model.dataset.test %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$multi_drug_start))
    print(table(hba1c.model.dataset.train$multi_drug_start, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$multi_drug_start))
    print(table(hba1c.model.dataset.test$multi_drug_start, hba1c.model.dataset.test$drugclass))
    
  }
  
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
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$timeprevcombo_less61))
    print(table(hba1c.model.dataset.train$timeprevcombo_less61, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$timeprevcombo_less61))
    print(table(hba1c.model.dataset.test$timeprevcombo_less61, hba1c.model.dataset.test$drugclass))
    
  }

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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$hb_extreme_53))
    print(table(hba1c.model.dataset.train$hb_extreme_53, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$hb_extreme_53))
    print(table(hba1c.model.dataset.test$hb_extreme_53, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(hb_extreme_53))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$prehba1c)))
    print(table(is.na(hba1c.model.dataset.train$prehba1c), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$prehba1c)))
    print(table(is.na(hba1c.model.dataset.test$prehba1c), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(prehba1c))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post HbA1c missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal), hba1c.model.dataset.test$drugclass))
    
  }
  
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
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Training cohort")
    print("################################################")
    print("Training Cohort")
    print(nrow(final.hba1c.model.dataset.train))
    print(table(final.hba1c.model.dataset.train$drugclass))
    
  }
  
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
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Testing cohort")
    print("################################################")
    print("Testing Cohort")
    print(nrow(final.hba1c.model.dataset.test))
    print(table(final.hba1c.model.dataset.test$drugclass))
    
  }
  
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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model")
    print("################################################")
    
  }
  
  weight.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(weight.dataset$multi_drug_start))
    print(table(weight.dataset$multi_drug_start, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(weight.dataset$timeprevcombo_less61))
    print(table(weight.dataset$timeprevcombo_less61, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(weight.dataset$hb_extreme_53))
    print(table(weight.dataset$hb_extreme_53, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(weight.dataset$prehba1c)))
    print(table(is.na(weight.dataset$prehba1c), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$preweight)))
    print(table(is.na(weight.dataset$preweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$postweight)))
    print(table(is.na(weight.dataset$postweight), weight.dataset$drugclass))
    
  }
  
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
      agetx, sex, t2dmduration, 
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
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model - final")
    print("################################################")
    print(nrow(final.weight.dataset))
    print(table(final.weight.dataset$drugclass))
    
  }
  
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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model")
    print("################################################")
    
  }
  
  discontinuation.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(discontinuation.dataset$multi_drug_start))
    print(table(discontinuation.dataset$multi_drug_start, discontinuation.dataset$drugclass))
    
  }

  discontinuation.dataset <- discontinuation.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(discontinuation.dataset$timeprevcombo_less61))
    print(table(discontinuation.dataset$timeprevcombo_less61, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(discontinuation.dataset$hb_extreme_53))
    print(table(discontinuation.dataset$hb_extreme_53, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(discontinuation.dataset$prehba1c)))
    print(table(is.na(discontinuation.dataset$prehba1c), discontinuation.dataset$drugclass))
    
  }
  
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
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model - final")
    print("################################################")
    print(nrow(final.discontinuation.dataset))
    print(table(final.discontinuation.dataset$drugclass))
    
  }
  
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
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model")
    print("################################################")
    
  }
  
  egfr.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(egfr.dataset$multi_drug_start))
    print(table(egfr.dataset$multi_drug_start, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(egfr.dataset$timeprevcombo_less61))
    print(table(egfr.dataset$timeprevcombo_less61, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(egfr.dataset$hb_extreme_53))
    print(table(egfr.dataset$hb_extreme_53, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$prehba1c)))
    print(table(is.na(egfr.dataset$prehba1c), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if baseline eGRF is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if baseline eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$preegfr)))
    print(table(is.na(egfr.dataset$preegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(preegfr))
  
  
  ################################################
  ##### Drop if post eGRF is missing
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(postegfr = ifelse(is.na(postegfr12m), postegfr6m, postegfr12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$postegfr)))
    print(table(is.na(egfr.dataset$postegfr), egfr.dataset$drugclass))
    
  }
  
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
      agetx, sex, t2dmduration, 
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
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model - final")
    print("################################################")
    print(nrow(final.egfr.dataset))
    print(table(final.egfr.dataset$drugclass))
    
  }
  
  if (dataset.type == "egfr.dataset") {
    return(final.egfr.dataset)
  }
  
  
}


####################
## Description:
##  - This file includes the framework for patient selection.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


# Loading libraries
library(tidyverse)
library(readr)
library(tableone)
library(stargazer)
library(corrplot)
library(caret)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1/datasets"

## make directory for outputs
dir.create(output_path)

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

load("/slade/CPRD_data/mastermind_2019/sglt2glp1/finalcprd19_analysisreadyv3.Rda")

dataset_raw <- cprd

###############################################################################
###############################################################################
############################### Filter Data ###################################
###############################################################################
###############################################################################

################################################
##### Select only SGLT-2 and GLP-1
################################################

#SGLT2 vs GLP1
cprd <- dataset_raw %>% filter(drugclass=="GLP1" | drugclass=="SGLT2") # nrow = 40708
table(cprd$drugclass)
# table(cprd$drugclass,cprd$yrdrugstart)

#       2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019
# GLP1   169  950 1835 2342 2170 2215 1825 1392 1353  984  977 1015  553
# SGLT2    0    0    0    0    0    0  876 2943 4499 4523 4375 4028 1684

# table(cprd$drugclass,cprd$drugline)
#          1    2    3    4    5    6    7    8    9
# GLP1    51 1134 3850 7084 4388 1115  145   13    0
# SGLT2  106 4168 6064 6482 3734 1776  533   62    3

################################################
##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
################################################

cprd <- cprd %>% mutate(timeprevcom_less61=if_else(timeprevcomb<=61,1,NA_real_))          
table(cprd$timeprevcom_less61) # removing 5455
cprd <- cprd %>% filter(is.na(timeprevcom_less61))
cprd$timeprevcom_less61 <- NULL

################################################
##### Drop if treated with insulin when starting new drug
################################################

table(cprd$INS) # 4488
cprd <- cprd %>% filter(INS==0)      

################################################
##### Drop if first-line treatment
################################################

table(cprd$drugline) # 157
cprd <- cprd %>% filter(drugline != 1)

################################################
##### Drop duplicates (i.e. started treatment on same day)
################################################

cprd$pated <- as.character(cprd$pated)
px <- data.frame(unique(cprd$pated)) #YES

cprddup <- cprd[cprd$pated %in% cprd$pated[duplicated(cprd$pated)],]
cprddup$dummy2 <- 1
cprddup <- cprddup %>% select(patid,dummy2)
cprddup <- cprddup[!duplicated(cprddup$patid),]

#keep only not starting multiple treatments
cprd <- merge(cprd,cprddup,all.x=TRUE,by="patid")
table(cprd$dummy2)
cprd <- cprd %>% filter(is.na(dummy2)) # remove 34

################################################
##### Drop if eGFR<45
################################################

cprd <- cprd %>% mutate(egfr45= if_else(egfr_ckdepi<45,1,0))
cprd <- cprd %>% mutate(egfr45= if_else(is.na(egfr45),0,egfr45))
table(cprd$egfr45)

cprd <- cprd %>% filter(egfr45==0) # remove 543

################################################
##### Drop if HbA1c not in 53 - 120 at baseline
################################################

cprd <- cprd %>% mutate(hb_extreme=if_else(prehba1cmmol<53|prehba1cmmol>120,1,0))
table(cprd$hb_extreme)
cprd <- cprd %>% filter(hb_extreme==0) # remove 1326

################################################
##### Drop if postHbA1c is missing
################################################

summary(cprd$posthba1c_final)
cprd <- cprd %>% filter(!is.na(posthba1c_final)) # remove 10657

################################################
##### Add variable for people that had both drugs
################################################

#Patients with valid data on each therapy    
sglt2glp1resp <- cprd[cprd$patid %in% cprd$patid[duplicated(cprd$patid)],]
# table(sglt2glp1resp$drugclass)   
#save(sglt2glp1resp,file="sglt2glp1resp.Rda")

#define marker in main dataset if data on both drugs available
patids <- sglt2glp1resp %>% 
  dplyr::select(patid) %>% 
  mutate(bothdrugs=1) %>% distinct()
cprd <- merge(cprd,patids,by="patid",all.x=T)
table(cprd$bothdrugs) # 2120

################################################
##### Drop if BMI is less than 5
################################################

cprd$prebmi <- ifelse(cprd$prebmi<5,NA,cprd$prebmi)
summary(cprd$prebmi)
cprd <- cprd %>% filter(!is.na(prebmi)) # remove 325

###############################################################################
###############################################################################
############################# Variable Prep ###################################
###############################################################################
###############################################################################

################################################
##### Code drugline as factor with 5+ group
################################################

cprd$drugline <- ifelse(cprd$drugline>5,5,cprd$drugline)
cprd$drugline <- factor(cprd$drugline)


cprd$druglinena <- ifelse(cprd$druglinena>5,5,cprd$druglinena)
cprd$druglinena <- factor(cprd$druglinena)

################################################
##### Code drugclass as factor of two 
################################################

#Define drugclass as 2 level factor
cprd$drugclass <- as.factor(as.character(cprd$drugclass))
# levels(cprd$drugclass)
# table(cprd$drugclass)

################################################
##### Code ncurrtx as factor of number of current therapies
################################################

cprd$ncurrtx <- factor(cprd$ncurrtx)

################################################
##### Clean up sex variable
################################################

cprd <- cprd %>% mutate(malesex=factor(sex))

################################################
##### Clean up Ethnicity
################################################
cprd <- cprd %>% mutate(ethnicityF=ifelse(cprd$ethnicity5==1,"White",
                                          ifelse(cprd$ethnicity5==2,"Mixed",
                                                 ifelse(cprd$ethnicity5==3,"Asian",
                                                        ifelse(cprd$ethnicity5==4,"Black",
                                                               ifelse(cprd$ethnicity5==5,"Other",
                                                                      ifelse(cprd$ethnicity5=="6"|cprd$ethnicity5=="7","Missing",NA)))))))
cprd$ethnicityF <- factor(cprd$ethnicityF, levels = c("White", "Mixed", "Asian", "Black", "Other", "Missing"))
cprd$ethnicityF <- relevel(cprd$ethnicityF,ref="White")

cprd$nonwhite <- fct_collapse(cprd$ethnicityF,WhiteUnknown=c("White","Missing"),NonWhite=c("Mixed","Asian","Black","Other"))

################################################
##### Code Cholesterol
################################################
cprd$chol <- cprd$pretc

################################################
##### Code SBP
################################################
cprd$sbp <- cprd$presys

################################################
##### Centre scale variables
################################################
cs. <- function(x) scale(x,center=TRUE,scale=TRUE)
cprd$t2dmdurationlog <- log(cprd$t2dmduration)
cprd$t2dmdurationcs <- as.numeric(cs.(cprd$t2dmduration))
cprd$t2dmdurationcslog <- as.numeric(cs.(log(cprd$t2dmduration)))

# hist(cprd$agetx)
cprd$agetxcs <- as.numeric(cs.(cprd$agetx))
cprd$agetxcslog <- as.numeric(cs.(log(cprd$agetx)))

# hist(cprd$precreat)
cprd$precreatcs <- as.numeric(cs.(cprd$precreat))
cprd$precreatcslog <- as.numeric(cs.(log(cprd$precreat)))

# hist(cprd$egfr_ckdepi)
cprd$egfr_ckdepics <- as.numeric(cs.(cprd$egfr_ckdepi))
cprd$egfr_ckdepicslog <- as.numeric(cs.(log(cprd$egfr_ckdepi)))

# hist(cprd$prebmi)
cprd$prebmics <- as.numeric(cs.(cprd$prebmi))
cprd$prebmicslog <- as.numeric(cs.(log(cprd$prebmi)))

cprd$preweightcs <- as.numeric(cs.(cprd$preweight))
cprd$preweightcslog <- as.numeric(cs.(log(cprd$preweight)))

cprd$preldlcs <- as.numeric(cs.(cprd$preldl))
cprd$preldlcslog <- as.numeric(cs.(log(cprd$preldl)))

cprd$prehdlcs <- as.numeric(cs.(cprd$prehdl))
cprd$prehdlcslog <- as.numeric(cs.(log(cprd$prehdl)))

cprd$pretrigcs <- as.numeric(cs.(cprd$pretrig))
cprd$pretriglogcs <- as.numeric(cs.(log(cprd$pretrig)))
cprd$pretriglog <- log(cprd$pretrig)


cprd$prealtcs <- as.numeric(cs.(cprd$prealt))
cprd$prealtlogcs <- as.numeric(cs.(log(cprd$prealt)))
cprd$prealtlog <- log(cprd$prealt)

cprd$preastcs <- as.numeric(cs.(cprd$preast))
cprd$preastlogcs <- as.numeric(cs.(log(cprd$preast)))

cprd$pregluccs <- as.numeric(cs.(cprd$pregluc))
cprd$pregluclogcs <- as.numeric(cs.(log(cprd$pregluc)))

cprd$prealbcs <- as.numeric(cs.(cprd$prealb))
cprd$prealblogcs <- as.numeric(cs.(log(cprd$prealb)))

cprd$prebilcs <- as.numeric(cs.(cprd$prebil))
cprd$prebillogcs <- as.numeric(cs.(log(cprd$prebil)))

cprd$preplateletscs <- as.numeric(cs.(cprd$preplatelets))
cprd$preplateletslogcs <- as.numeric(cs.(log(cprd$preplatelets)))

cprd$cholcs <- as.numeric(cs.(cprd$chol))

cprd$sbpcs <- as.numeric(cs.(cprd$sbp))

################################################
##### Remove all values before 2014
################################################

# cprd <- cprd %>%
#   filter(yrdrugstart > 2014) # removed 4648

###############################################################################
###############################################################################
########################## Add comorbidities ##################################
###############################################################################
###############################################################################

# Add comorbidity
load("/slade/CPRD_data/mastermind_2019/sglt2glp1/cprd.comorbidity.Rda") 

#Merge in comorbidity
cprd.comorbidity$patid <- NULL
cprd <- merge(cprd,cprd.comorbidity,by="pateddrug",all.x=TRUE)

################################################
##### Clean Smoking variable
################################################

#Smoking binary
cprd <- cprd %>% mutate(smok=ifelse(Category=="Active smoker",1,0))
cprd <- cprd %>% mutate(smok=ifelse(is.na(smok),0,smok))

################################################
##### Include Established ASCVD ()
################################################

# Established ASCVD: MI, ischemic stroke, unstable angina with ECG changes, myocardial ischemia on imaging or stress test, revascularisation of coronary, carotid or peripheral arteries.
# High risk ASCVD: Age 55 + Left ventricular hypertrophy or coronary, carotid, lower extremity artery stenosis >50%
# HF especially HFrEF (EF <45%)
# CKD - eGFR 30-60 or UACR>30 (preferably >300) - says 300 in main flow.

cprd <- cprd %>% mutate(ascvd=if_else(predrug.earliest.mi==1|predrug.earliest.stroke==1|predrug.earliest.revasc==1|predrug.earliest.pad==1|predrug.earliest.ihd==1,1,0))
# table(cprd$ascvd)

cprd <- cprd %>% mutate(ascvd.hf=if_else(ascvd==1|predrug.earliest.heartfailure==1,1,0)) 
# table(cprd$ascvd.hf)

cprd <- cprd %>% mutate(ckd=if_else(predrug.earliest.ckd5==1|egfr_ckdepi<60,1,0)) 
cprd <- cprd %>% mutate(ckd=if_else(is.na(ckd),0,ckd))
# table(cprd$ckd)

cprd <- cprd %>% mutate(ascvd.hf.ckd=if_else(ascvd.hf==1|ckd==1,1,0)) 
# table(cprd$ascvd.hf.ckd)


################################################
##### Include risk and non-coronary CVD
################################################

cprd$age <- cprd$agetx

#Step 1:
#Underlying risk  and non-coronary CVD at current age and in 10 years time
#CHD can also be calculated and is additive to estimate total fatal CVD
#UK is low risk
# chol <- 8
# smok <- 1
# sbp <- 180
# age <- 55

alpha_male <- -26.7; alpha_female <- -31.0
p_male <- 5.64; p_female <- 6.62
#CVD, current

cprd <- cprd %>% mutate(S0_current_male = exp(-(exp(alpha_male))*(age-20)^p_male))
cprd <- cprd %>% mutate(S0_current_female = exp(-(exp(alpha_female))*(age-20)^p_female))

#CVD, 10 year      
cprd <- cprd %>% mutate(S0_10yr_male = exp(-(exp(alpha_male))*(age-10)^p_male))
cprd <- cprd %>% mutate(S0_10yr_female = exp(-(exp(alpha_female))*(age-10)^p_female))

#Step 2: linear predictor 

b_chol <- 0.02; b_sbp <- 0.022; b_smok <- 0.63
cprd <- cprd %>% mutate(lp = b_chol*(chol-6) + b_sbp*(sbp-120) + b_smok*smok)

#Step 3: combine

cprd <- cprd %>% mutate(risk_10_male.cvd = 1-((S0_10yr_male^exp(lp))/(S0_current_male^exp(lp))))
cprd <- cprd %>% mutate(risk_10_female.cvd = 1-((S0_10yr_female^exp(lp))/(S0_current_female^exp(lp))))

# summary(cprd$risk_10_female.cvd)
# summary(cprd$risk_10_male.cvd)

cprd <- cprd %>% mutate(score=ifelse(malesex==1,risk_10_male.cvd,risk_10_female.cvd))

cprd <- cprd %>% mutate(score.high=if_else(score<0.05|is.na(score),0,1))
# table(cprd$score.high)

# table(cprd$score.high,cprd$ascvd.hf)

# table(cprd$score.high,cprd$ascvd.hf.ckd)

cprd <- cprd %>% mutate(ascvd.hf.ckd.scorehigh=if_else(ascvd.hf.ckd==1|score.high==1,1,0))

cprd <- cprd %>% mutate(score.excl.mi=ifelse(predrug.earliest.mi==1,NA,score*100))
cprd <- cprd %>% mutate(score.high.excl.mi=ifelse(predrug.earliest.mi==1,NA,score.high))
# table(cprd$score.high.excl.mi)

cprd <- cprd %>% mutate(ascvd.hf.ckd.scorehigh=if_else(ascvd.hf.ckd==1|score.high.excl.mi==1,1,0))

################################################
##### Turn variables into categorical variables
################################################

cprd <- cprd %>%
  mutate(
    sglt2subtype = factor(sglt2subtype),
    glp1subtype = factor(glp1subtype),
    Category = factor(Category),
    predrug.earliest.stroke = factor(predrug.earliest.stroke),
    predrug.5yrrecent.stroke = factor(predrug.5yrrecent.stroke),
    predrug.earliest.tia = factor(predrug.earliest.tia),
    predrug.5yrrecent.tia = factor(predrug.5yrrecent.tia),
    predrug.earliest.revasc = factor(predrug.earliest.revasc),
    predrug.5yrrecent.revasc = factor(predrug.5yrrecent.revasc),
    predrug.earliest.pad = factor(predrug.earliest.pad),
    predrug.5yrrecent.pad = factor(predrug.5yrrecent.pad),
    predrug.earliest.neuropathy = factor(predrug.earliest.neuropathy),
    predrug.5yrrecent.neuropathy = factor(predrug.5yrrecent.neuropathy),
    predrug.earliest.nephropathy = factor(predrug.earliest.nephropathy),
    predrug.5yrrecent.nephropathy = factor(predrug.5yrrecent.nephropathy),
    predrug.earliest.retinopathy = factor(predrug.earliest.retinopathy),
    predrug.5yrrecent.retinopathy = factor(predrug.5yrrecent.retinopathy),
    predrug.earliest.ihd = factor(predrug.earliest.ihd),
    predrug.5yrrecent.ihd = factor(predrug.5yrrecent.ihd),
    predrug.earliest.hypertension = factor(predrug.earliest.hypertension),
    predrug.5yrrecent.hypertension = factor(predrug.5yrrecent.hypertension),
    predrug.earliest.heartfailure = factor(predrug.earliest.heartfailure),
    predrug.5yrrecent.heartfailure = factor(predrug.5yrrecent.heartfailure),
    predrug.earliest.ckd5 = factor(predrug.earliest.ckd5),
    predrug.5yrrecent.ckd5 = factor(predrug.5yrrecent.ckd5),
    predrug.earliest.af = factor(predrug.earliest.af),
    predrug.5yrrecent.af = factor(predrug.5yrrecent.af),
    predrug.earliest.angina = factor(predrug.earliest.angina),
    predrug.5yrrecent.angina = factor(predrug.5yrrecent.angina),
    predrug.earliest.cld = factor(predrug.earliest.cld),
    predrug.5yrrecent.cld = factor(predrug.5yrrecent.cld),
    predrug.earliest.mi = factor(predrug.earliest.mi),
    predrug.5yrrecent.mi = factor(predrug.5yrrecent.mi)
  )






final.all <- cprd

###############################################################################
###############################################################################
##################### Variables chosen for modelling ##########################
###############################################################################
###############################################################################

# ## List variables for the table
# Vars <- c("drugline",       # "Brings to my attention drugline is potentially wrong for anyone who started ANY drug prior
#           "druglinena",     #   to cprd2 data. Make a list of these patients (potentially exclide as drugline == NA)
#           "ncurrtx",        # Number of current therapies
#           "yrdrugstart",    # Year of drug start
#           "agedx",          # Age at diagnosis?
#           "agetx",          # Age at treatment
#           "gender",         # gender var 1/2
#           "malesex",        # sex male-1
#           "glp1subtype",    # subtypes of GLP-1
#           "sglt2subtype",   # subtypes of SGLT-2
#           "t2dmduration",   # duration of treatment
#           "ethnicity5",     # Ethnicity in 5 categories
#           "adherence_t")    # Adherence to therapy
# ## List of variables which are factors
# FactorVars <- c("drugline","druglinena","ncurrtx","gender","malesex","glp1subtype","sglt2subtype","ethnicity5")
# ## BY MATCHED / ID Create the table
# DemographicsTable <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=final.all)
# ## Export as an html document
# # stargazer(print(DemographicsTable, dropEqual = T), out = paste0(output_path, "/DemographicsTableStrat.html"))
# 
# ## List variables for the table
# Vars <- c("precreat","egfr","egfr_ckdepi","egfr_cg",
#           "preldl","prehdl","pretrig","preweight","prebmi",
#           "preast","prealt","preplatelets","prebil","prealb","api")    ## List of variables which are factors
# ## BY MATCHED / ID Create the table
# Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=final.all)
# ## Export as an html document
# # stargazer(print(Table, dropEqual = T))
# 
# ## List variables for the table
# Vars <- c("prehba1cmmol","pregluc","prehba1cchange","iniresp_final","posthba1c_final") #"prehba1cmmol","posthba1c6mmmol","posthba1c12mmmol","iniresp6m","iniresp12m","iniresp6mmmol","iniresp12mmmol"
# ## BY MATCHED / ID Create the table
# Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=final.all)
# ## Export as an html document
# # stargazer(print(Table, dropEqual = T))
# 
# ## List variables for the table
# Vars <- c("fldrug","fltimeon","flt2dmduration","flprehba1cmmol","fliniresp_final","flresid","fladherence_t","timesince1stline")
# FactorVars <- c("fldrug")
# ## BY MATCHED / ID Create the table
# Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=final.all)
# ## Export as an html document
# # stargazer(print(Table, dropEqual = T))
# 
# 
# test <- final.all %>% select("agedx","agetx","precreat","egfr_ckdepi",
#                              "preldl","prehdl","pretrig","preweight","prebmi",
#                              "preast","prealt","preplatelets","prebil","prealb","api","t2dmduration","adherence_t",
#                              "prehba1cmmol","pregluc","prehba1cchange","fladherence_t","fliniresp_final","flt2dmduration","flprehba1cmmol","timesince1stline","flposthba1c_final")
# descrCor <- cor(test, use="na.or.complete")
# # print(descrCor)
# # corrplot(descrCor, order = "FPC", method = "color", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))
# # highlyCorrelated <- findCorrelation(descrCor, cutoff=0.5)
# # highlyCorCol <- colnames(test)[highlyCorrelated]
# # highlyCorCol


#################


final.all <- final.all %>% select(
  patid, pateddrug, bothdrugs,
  ################ Outcome HbA1c
  posthba1c_final,
  # posthba1c6mmmol, posthba1c12mmmol, posthba1c_final6m, hba1cdate6m, hba1cdate12m, respvalid6m,
  # respvalid12m, hba1cdate_final, hba1cdate_final6m,  # removing variants of outcome
  ################ Patient details
  drugclass, ncurrtx, drugline, yrdrugstart, sglt2subtype, glp1subtype, 
  t2dmduration, agetx, malesex, Category, ethnicityF, hba1cmonth,
  ################ Other drugs currently taken
  # MFN, SU, Acarbose, GLP1, Glinide, DPP4, INS, SGLT2, TZD,  # not included betcause we have ncurrtx
  ################ Baseline biomarkers
  # precreat, precreatmgdl, precreatcs, precreatcslog, #not included because eGFR uses creat for the calculation
  # preldl, preldlcs, preldlcslog, preldl, preldlcs, preldlcslog, pretc, chol, cholcs, #calculations in CPRD might not match the equation.
  # pretrig, pretriglog, pretrigcs, pretriglogcs, # should be fasting but CPRD probably wont be
  # pregluc, pregluccs, pregluclogcs, # NEEDS to be fasting but probably wont be. no way to know
  # sbp, sbpcs, # removed because it is used in CSV risk
  preweight, height, prebmi, # preweightcs, preweightcslog, prebmics, prebmicslog, # removing centering and logs
  preast, prealt, # preastcs, preastlogcs, prealtcs, prealtlog, prealtlogcs, # removing centering and logs
  egfr_ckdepi, # egfr, egfr_cg, egfr_cg2, egfr_ckdepi2, egfr2, egfr45, egfr_ckdepics, egfr_ckdepicslog, # removing all variants of egfr_ckdepi
  prealb, prebil, # prealbcs, prealblogcs, prebilcs, prebillogcs, # removing centering and logs
  preplatelets, # preplateletscs, preplateletslogcs, # removing centering and logs
  presys,
  prehdl, # prehdlcs, prehdlcslog, # removing centering and logs
  prehba1cmmol, # prehba1cpc, prehba1cdate, prehba1cchange, iniresp_final,  # removing information about prehba1c
  ################ Derived variables
  score, # risk of CVD
  ################ Diabetes microvascular complications and comorbidities
  predrug.earliest.stroke, predrug.5yrrecent.stroke, predrug.earliest.tia, predrug.5yrrecent.tia,
  predrug.earliest.revasc, predrug.5yrrecent.revasc, predrug.earliest.pad, predrug.5yrrecent.pad,
  predrug.earliest.neuropathy, predrug.5yrrecent.neuropathy, predrug.earliest.nephropathy,
  predrug.5yrrecent.nephropathy, predrug.earliest.retinopathy, predrug.5yrrecent.retinopathy,
  predrug.earliest.ihd, predrug.5yrrecent.ihd, predrug.earliest.hypertension, predrug.5yrrecent.hypertension,
  predrug.earliest.heartfailure, predrug.5yrrecent.heartfailure, predrug.earliest.ckd5, predrug.5yrrecent.ckd5,
  predrug.earliest.af, predrug.5yrrecent.af, predrug.earliest.angina, predrug.5yrrecent.angina,
  predrug.earliest.cld, predrug.5yrrecent.cld, predrug.earliest.mi, predrug.5yrrecent.mi,              
  ################ First line treatment details
  # flt2dmduration, flprehba1cmmol, fliniresp_final, flposthba1c_final, 
  # flstoptime, flstop, # fldate, fldrugstopdate, flresid, # removing derivatives of first-line vars
  ################ Post Therapy side-effects
  # postsys6m, postsys12m, postweight6m, postweight12m, postdrug.first.stroke,
  # postdrug.first.tia, postdrug.first.revasc, postdrug.first.pad, postdrug.first.neuropathy,
  # postdrug.first.nephropathy, postdrug.first.retinopathy, postdrug.first.ihd,
  # postdrug.first.hypertension, postdrug.first.heartfailure, postdrug.first.ckd5,
  # postdrug.first.af, postdrug.first.angina, postdrug.first.cld, postdrug.first.mi,
  ################ Discontinuation
  # stopdrugdate, stopdrug3m_6mFU, stopdrug3m_3mFU, stopdrug6m_6mFU, stopdrug6m_3mFU, stopdrug12m_6mFU, stopdrug12m_3mFU,
  ################ Adherence
  # adherence_t, fladherence_t, flnonadh, nonadh,
  ################ Time difference
  # timediffprecreat, timediffpreldl, timediffprehdl, timediffpretrig, timediffpregluc, timediffpreweight, timediffprebmi, 
  # timediffpreast, timediffprealt, timediffpretc, timediffpresys, timediffpre, timediffpost6m, timediffpost12m, 
  # timediffpost15m, timediffprealb, timediffprebil, timediffpreplatelets, timetolastpx, timetodrugstop, timediffpre91, 
  # timesince1stline, tdiff6m, tdiff12m, tdiff6m61, tdiff12m61,
  
  
  ################ Extra Variables Removed
  # pated, frd2, crd2, druglinena, datedrug, dcorder, drugstopdate, timetochange, timeprevcomb,
  # nextdrugchange, ptdrugclass, nextdrugchangedate, dcstartdate, dcstopdate, startdrugdate, drugcombo, dpp4subtype,
  # susubtype, tzdsubtype, durn_drug, bestdiagdate, gender, eth5, ethcode, ethnicitylong, ethnicity5, nextremdrug,
  # hba1cmmolm_pre1, hba1ctimediff_pre1, hba1cmmolm_pre2, hba1ctimediff_pre2, method, fib4elig, agecat, platecat,
  # yr6mdrugresp, prevdrugstopdate, mindateswitch, switch, dummy, sex, ethblack, api_age, api_pl, api, dummy2,
  # hb_extreme, nonwhite, t2dmdurationlog, t2dmdurationcs, t2dmdurationcslog, agetxcs, agetxcslog,
  # yob, tod, deathdate, indexdate, eth16, smok, ascvd, ascvd.hf, ckd, ascvd.hf.ckd, age, S0_current_male,
  # S0_current_female, S0_10yr_male, S0_10yr_female, lp, risk_10_male.cvd, risk_10_female.cvd, score.high,
  # ascvd.hf.ckd.scorehigh, score.excl.mi, score.high.excl.mi, timeprevcom_less61, agedx,
  
)


###############################################################################
###############################################################################
################# Define Training/Validation Population #######################
###############################################################################
###############################################################################

################################################
##### Total number of individuals: 16443(9964)
################################################

save(final.all,file=paste0(output_path, "/cprd_19_sglt2glp1_allcohort.Rda"))

################################################
##### Define Development dataset (60% on each drug): 9866(5978)
################################################
set.seed(8731)
cprd.dev <- final.all %>% group_by(drugclass) %>% sample_frac(.6)
cprd.dev <- data.frame(cprd.dev)

final.dev <- cprd.dev

save(final.dev,file=paste0(output_path, "/cprd_19_sglt2glp1_devcohort.Rda"))

################################################
##### Define Validation dataset (other 40% on each drug): 6577(3986)
################################################

cprd.val <- subset(final.all, !(pateddrug %in% cprd.dev$pateddrug))

final.val <- cprd.val

save(final.val,file=paste0(output_path, "/cprd_19_sglt2glp1_valcohort.Rda"))

################################################
##### Define Dataset with both drugs taken: 2084(3986)
################################################

#save alternative validation approach
final.bothdrugs <- final.all %>% dplyr::filter(bothdrugs==1)
final.onedrug <- final.all %>% dplyr::filter(is.na(bothdrugs))

save(final.bothdrugs,file=paste0(output_path, "/cprd_19_sglt2glp1_bothdrugs.Rda"))
save(final.onedrug,file=paste0(output_path, "/cprd_19_sglt2glp1_onedrug.Rda"))



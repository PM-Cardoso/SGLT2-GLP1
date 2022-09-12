####################
## Description:
##  - In this file we plot differential treatment effect for all variables from 4.4 model
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")


library(tidyverse)
library(bartMachine)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)


## make directory for outputs
dir.create(paste0(output_path, "/Comparison"))

## make directory for outputs
dir.create(paste0(output_path, "/Comparison/Model_4_4"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
######################### Read Data / Model In ################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))
# name: final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))

###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################

source("0.1.slade_functions.R")

#############################
### Treatment effect for different features
#############################

bart_model_final <- readRDS(paste0(output_path, "/Final_model/With_grf_no_prop/bart_model_final.rds"))


dataset.dev <- final.dev %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

## all variables
if (class(try(
  
  effects_summary_dev <- readRDS(paste0(output_path, "/Comparison/Model_4_4/effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev <- diff_treatment_effect(bart_model_final, dataset.dev, 25)
    
  saveRDS(effects_summary_dev, paste0(output_path, "/Comparison/Model_4_4/effects_summary_dev.rds"))
}


## egfr_ckdepi

# plot treatment effect + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["egfr_ckdepi"]], "egfr_ckdepi", "eGFR")

## prealt

# plot treatment effect + histogram marginal
plot.alt.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prealt"]], "prealt", "ALT")

## prehba1cmmol

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prehba1cmmol"]], "prehba1cmmol", "HbA1c")

## score

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["score"]], "score", "CVD Score")

## agetx

# plot treatment effect + histogram marginal
plot.age.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["agetx"]], "agetx", "Age")

## prehdl

# plot treatment effect + histogram marginal
plot.hdl.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prehdl"]], "prehdl", "HDL")

## prebmi

# plot treatment effect + histogram marginal
plot.bmi.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prebmi"]], "prebmi", "BMI")

## prebil

# plot treatment effect + histogram marginal
plot.bil.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prebil"]], "prebil", "BIL")

## preplatelets

# plot treatment effect + histogram marginal
plot.platelets.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["preplatelets"]], "preplatelets", "Platelets")

## t2dmduration

# plot treatment effect + histogram marginal
plot.t2dmduration.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["t2dmduration"]], "t2dmduration", "t2dmduration")

## prealb

# plot treatment effect + histogram marginal
plot.alb.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prealb"]], "prealb", "ALB")

## presys

# plot treatment effect + histogram marginal
plot.sys.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["presys"]], "presys", "SYS")

## preast

# plot treatment effect + histogram marginal
plot.ast.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["preast"]], "preast", "AST")

## drugline

# plot treatment effect + histogram marginal
plot.drugline.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["drugline"]], "drugline", "Drugline")

## Category Smoker

# plot treatment effect + histogram marginal
plot.category.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["Category"]], "Category", "Smoking")

## ncurrtx

# plot treatment effect + histogram marginal
plot.ncurrtx.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["ncurrtx"]], "ncurrtx", "ncurrtx")

## malesex

# plot treatment effect + histogram marginal
plot.malesex.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["malesex"]], "malesex", "Sex")


#### PDF with all the plots


pdf(file = "Plots/6.4.differential_treatment_effect.pdf")
plot.age.diff.marg
plot.ast.diff.marg
plot.alb.diff.marg
plot.t2dmduration.diff.marg
plot.platelets.diff.marg
plot.bil.diff.marg
plot.bmi.diff.marg
plot.hdl.diff.marg
plot.sys.diff.marg
plot.score.diff.marg
plot.hba1c.diff.marg
plot.alt.diff.marg
plot.egfr.diff.marg
plot.drugline.diff.marg
plot.category.diff.marg
plot.ncurrtx.diff.marg
plot.malesex.diff.marg
dev.off()



########
# Average patient
average.patient <- cbind(
  drugclass = "SGLT2",
  egfr_ckdepi = mean(dataset.dev$egfr_ckdepi, na.rm = TRUE),
  hba1cmonth = 6,
  prealt = mean(dataset.dev$prealt, na.rm = TRUE),
  prehba1cmmol = mean(dataset.dev$prehba1cmmol, na.rm = TRUE),
  score = mean(dataset.dev$score, na.rm = TRUE),
  Category = "Non-smoker",
  drugline = "2",
  ncurrtx = "1",
  yrdrugstart = 2016,
  agetx = mean(dataset.dev$agetx, na.rm = TRUE),
  malesex = "1",
  prehdl = mean(dataset.dev$prehdl, na.rm = TRUE),
  prebmi = mean(dataset.dev$prebmi, na.rm = TRUE),
  prebil = mean(dataset.dev$prebil, na.rm = TRUE),
  preplatelets = mean(dataset.dev$preplatelets, na.rm = TRUE),
  t2dmduration = mean(dataset.dev$t2dmduration, na.rm = TRUE),
  prealb = mean(dataset.dev$prealb, na.rm = TRUE),
  presys = mean(dataset.dev$presys, na.rm = TRUE),
  preast = mean(dataset.dev$preast, na.rm = TRUE)
) %>%
  as.data.frame() %>%
  mutate(drugclass = factor(drugclass, levels = levels(dataset.dev$drugclass)),
         egfr_ckdepi = as.numeric(egfr_ckdepi),
         hba1cmonth = as.numeric(hba1cmonth),
         prealt = as.numeric(prealt),
         prehba1cmmol = as.numeric(prehba1cmmol),
         score = as.numeric(score),
         Category = factor(Category, levels = levels(dataset.dev$Category)),
         drugline = factor(drugline, levels = levels(dataset.dev$drugline)),
         ncurrtx = factor(ncurrtx, levels = levels(dataset.dev$ncurrtx)),
         yrdrugstart = as.numeric(yrdrugstart),
         agetx = as.numeric(agetx),
         malesex = factor(malesex, levels = levels(dataset.dev$malesex)),
         prehdl = as.numeric(prehdl),
         prebmi = as.numeric(prebmi),
         prebil = as.numeric(prebil),
         preplatelets = as.numeric(preplatelets),
         t2dmduration = as.numeric(t2dmduration),
         prealb = as.numeric(prealb),
         presys = as.numeric(presys),
         preast = as.numeric(preast)
  )
       
effects_summary_patient <- diff_treatment_effect(bart_model_final, average.patient, 25)


## egfr_ckdepi

# plot treatment effect + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["egfr_ckdepi"]], "egfr_ckdepi", "eGFR")

## prealt

# plot treatment effect + histogram marginal
plot.alt.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prealt"]], "prealt", "ALT")

## prehba1cmmol

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prehba1cmmol"]], "prehba1cmmol", "HbA1c")

## score

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["score"]], "score", "CVD Score")

## agetx

# plot treatment effect + histogram marginal
plot.age.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["agetx"]], "agetx", "Age")

## prehdl

# plot treatment effect + histogram marginal
plot.hdl.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prehdl"]], "prehdl", "HDL")

## prebmi

# plot treatment effect + histogram marginal
plot.bmi.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prebmi"]], "prebmi", "BMI")

## prebil

# plot treatment effect + histogram marginal
plot.bil.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prebil"]], "prebil", "BIL")

## preplatelets

# plot treatment effect + histogram marginal
plot.platelets.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["preplatelets"]], "preplatelets", "Platelets")

## t2dmduration

# plot treatment effect + histogram marginal
plot.t2dmduration.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["t2dmduration"]], "t2dmduration", "t2dmduration")

## prealb

# plot treatment effect + histogram marginal
plot.alb.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["prealb"]], "prealb", "ALB")

## presys

# plot treatment effect + histogram marginal
plot.sys.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["presys"]], "presys", "SYS")

## preast

# plot treatment effect + histogram marginal
plot.ast.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["preast"]], "preast", "AST")

## drugline

# plot treatment effect + histogram marginal
plot.drugline.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["drugline"]], "drugline", "Drugline")

## Category Smoker

# plot treatment effect + histogram marginal
plot.category.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["Category"]], "Category", "Smoking")

## ncurrtx

# plot treatment effect + histogram marginal
plot.ncurrtx.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["ncurrtx"]], "ncurrtx", "ncurrtx")

## malesex

# plot treatment effect + histogram marginal
plot.malesex.diff.marg <- plot_diff_treatment_effect(effects_summary_patient[["malesex"]], "malesex", "Sex")




#### PDF with all the plots


pdf(file = "Plots/6.4.patient_differential_treatment_effect.pdf")
plot.age.diff.marg
plot.ast.diff.marg
plot.alb.diff.marg
plot.t2dmduration.diff.marg
plot.platelets.diff.marg
plot.bil.diff.marg
plot.bmi.diff.marg
plot.hdl.diff.marg
plot.sys.diff.marg
plot.score.diff.marg
plot.hba1c.diff.marg
plot.alt.diff.marg
plot.egfr.diff.marg
plot.drugline.diff.marg
plot.category.diff.marg
plot.ncurrtx.diff.marg
plot.malesex.diff.marg
dev.off()










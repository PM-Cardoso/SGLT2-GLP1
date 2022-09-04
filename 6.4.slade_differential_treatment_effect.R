####################
## Description:
##  - In this file we plot differential treatment effect for all variables
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
dir.create(paste0(output_path, "/Comparison/Model_7_1"))

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

source("11.slade_functions.R")

#############################
### Treatment effect for different features
#############################

bart_model_final <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/bart_model_final.rds"))


dataset.dev <- final.dev %>%
  filter(yrdrugstart > 2012) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

## egfr_ckdepi

# calculate differential effects
if (class(try(
  
  effects_summary_dev_egfr <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_egfr.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_egfr <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "egfr_ckdepi", 11)
  
  saveRDS(effects_summary_dev_egfr, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_egfr.rds"))
  
}

# plot treatment effect + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_egfr, dataset.dev, "egfr_ckdepi", "eGFR")

## prealt

# calculate differential effects
if (class(try(
  
  effects_summary_dev_alt <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_alt.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_alt <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prealt", 11)
  
  saveRDS(effects_summary_dev_alt, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_alt.rds"))
  
}

# plot treatment effect + histogram marginal
plot.alt.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_alt, dataset.dev, "prealt", "ALT")


## prehba1cmmol

# calculate differential effects
if (class(try(
  
  effects_summary_dev_hba1c <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_hba1c.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_hba1c <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prehba1cmmol", 11)
  
  saveRDS(effects_summary_dev_hba1c, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_hba1c.rds"))
  
}

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_hba1c, dataset.dev, "prehba1cmmol", "HbA1c")


## score

# calculate differential effects
if (class(try(
  
  effects_summary_dev_score <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_score.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_score <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "score", 11)
  
  saveRDS(effects_summary_dev_score, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_score.rds"))
  
}

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_score, dataset.dev, "score", "CVD Score")


## agetx

# calculate differential effects
if (class(try(
  
  effects_summary_dev_age <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_age.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_age <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "agetx", 11)
  
  saveRDS(effects_summary_dev_age, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_age.rds"))
  
}

# plot treatment effect + histogram marginal
plot.age.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_age, dataset.dev, "agetx", "Age")

## prehdl

# calculate differential effects
if (class(try(
  
  effects_summary_dev_hdl <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_hdl.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_hdl <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prehdl", 11)
  
  saveRDS(effects_summary_dev_hdl, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_hdl.rds"))
  
}

# plot treatment effect + histogram marginal
plot.hdl.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_hdl, dataset.dev, "prehdl", "HDL")

## prebmi

# calculate differential effects
if (class(try(
  
  effects_summary_dev_bmi <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_bmi.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_bmi <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prebmi", 11)
  
  saveRDS(effects_summary_dev_bmi, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_bmi.rds"))
  
}

# plot treatment effect + histogram marginal
plot.bmi.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_bmi, dataset.dev, "prebmi", "BMI")

## prebil

# calculate differential effects
if (class(try(
  
  effects_summary_dev_bil <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_bil.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_bil <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prebil", 11)
  
  saveRDS(effects_summary_dev_bil, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_bil.rds"))
  
}

# plot treatment effect + histogram marginal
plot.bil.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_bil, dataset.dev, "prebil", "BIL")

## preplatelets

# calculate differential effects
if (class(try(
  
  effects_summary_dev_platelets <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_platelets.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_platelets <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "preplatelets", 11)
  
  saveRDS(effects_summary_dev_platelets, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_platelets.rds"))
  
}

# plot treatment effect + histogram marginal
plot.platelets.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_platelets, dataset.dev, "preplatelets", "Platelets")


## t2dmduration

# calculate differential effects
if (class(try(
  
  effects_summary_dev_t2dmduration <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_t2dmduration.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_t2dmduration <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "t2dmduration", 11)
  
  saveRDS(effects_summary_dev_t2dmduration, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_t2dmduration.rds"))
  
}

# plot treatment effect + histogram marginal
plot.t2dmduration.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_t2dmduration, dataset.dev, "t2dmduration", "t2dmduration")


## prealb

# calculate differential effects
if (class(try(
  
  effects_summary_dev_alb <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_alb.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_alb <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "prealb", 11)
  
  saveRDS(effects_summary_dev_alb, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_alb.rds"))
  
}

# plot treatment effect + histogram marginal
plot.alb.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_alb, dataset.dev, "prealb", "ALB")


## presys

# calculate differential effects
if (class(try(
  
  effects_summary_dev_sys <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_sys.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_sys <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "presys", 11)
  
  saveRDS(effects_summary_dev_sys, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_sys.rds"))
  
}

# plot treatment effect + histogram marginal
plot.sys.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_sys, dataset.dev, "presys", "SYS")


## preast

# calculate differential effects
if (class(try(
  
  effects_summary_dev_ast <- readRDS(paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_ast.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev_ast <- calc_diff_treatment_effect(bart_model_final, dataset.dev, "preast", 11)
  
  saveRDS(effects_summary_dev_ast, paste0(output_path, "/Comparison/Model_7_1/effects_summary_dev_ast.rds"))
  
}

# plot treatment effect + histogram marginal
plot.ast.diff.marg <- plot_diff_treatment_effect(effects_summary_dev_ast, dataset.dev, "preast", "AST")




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
dev.off()






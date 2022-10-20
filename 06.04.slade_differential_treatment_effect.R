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
dir.create(paste0(output_path, "/Comparison/Model_5"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
######################### Read Data / Model In ################################
###############################################################################
###############################################################################

# name: final.all.extra.vars
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_allcohort.Rda"))
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

bart_model_final <- readRDS(paste0(output_path, "/Final_model/model_7/bart_model_final.rds"))


dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

# ## all variables
# if (class(try(
#   
#   effects_summary_dev <- readRDS(paste0(output_path, "/Comparison/Model_5/effects_summary_dev.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   effects_summary_dev <- diff_treatment_effect(bart_model_final, dataset.dev, 25)
#     
#   saveRDS(effects_summary_dev, paste0(output_path, "/Comparison/Model_5/effects_summary_dev.rds"))
# }
# 
# 
# ## egfr_ckdepi
# 
# # plot treatment effect + histogram marginal
# plot.egfr.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["egfr_ckdepi"]], "egfr_ckdepi", "eGFR", 4)
# 
# ## prealt
# 
# # plot treatment effect + histogram marginal
# plot.alt.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prealt"]], "prealt", "ALT", 4)
# 
# ## prehba1cmmol
# 
# # plot treatment effect + histogram marginal
# plot.hba1c.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prehba1cmmol"]], "prehba1cmmol", "HbA1c", 4)
# 
# ## score.excl.mi
# 
# # plot treatment effect + histogram marginal
# plot.score.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["score.excl.mi"]], "score.excl.mi", "CVD Score", 4)
# 
# ## agetx
# 
# # plot treatment effect + histogram marginal
# plot.age.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["agetx"]], "agetx", "Age", 4)
# 
# ## prehdl
# 
# # plot treatment effect + histogram marginal
# plot.hdl.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prehdl"]], "prehdl", "HDL", 4)
# 
# ## prebmi
# 
# # plot treatment effect + histogram marginal
# plot.bmi.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prebmi"]], "prebmi", "BMI", 4)
# 
# ## prebil
# 
# # plot treatment effect + histogram marginal
# plot.bil.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prebil"]], "prebil", "BIL", 4)
# 
# ## preplatelets
# 
# # plot treatment effect + histogram marginal
# plot.platelets.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["preplatelets"]], "preplatelets", "Platelets", 4)
# 
# ## t2dmduration
# 
# # plot treatment effect + histogram marginal
# plot.t2dmduration.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["t2dmduration"]], "t2dmduration", "t2dmduration", 4)
# 
# ## prealb
# 
# # plot treatment effect + histogram marginal
# plot.alb.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["prealb"]], "prealb", "ALB", 4)
# 
# ## presys
# 
# # plot treatment effect + histogram marginal
# plot.sys.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["presys"]], "presys", "SYS", 4)
# 
# ## preast
# 
# # plot treatment effect + histogram marginal
# plot.ast.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["preast"]], "preast", "AST", 4)
# 
# ## drugline
# 
# # plot treatment effect + histogram marginal
# plot.drugline.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["drugline"]], "drugline", "Drugline")
# 
# ## Category Smoker
# 
# # plot treatment effect + histogram marginal
# plot.category.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["Category"]], "Category", "Smoking")
# 
# ## ncurrtx
# 
# # plot treatment effect + histogram marginal
# plot.ncurrtx.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["ncurrtx"]], "ncurrtx", "ncurrtx")
# 
# ## malesex
# 
# # plot treatment effect + histogram marginal
# plot.malesex.diff.marg <- plot_diff_treatment_effect(effects_summary_dev[["malesex"]], "malesex", "Sex")
# 
# 
# #### PDF with all the plots
# 
# 
# pdf(file = "Plots/6.4.differential_treatment_effect.pdf")
# plot.age.diff.marg
# plot.ast.diff.marg
# plot.alb.diff.marg
# plot.t2dmduration.diff.marg
# plot.platelets.diff.marg
# plot.bil.diff.marg
# plot.bmi.diff.marg
# plot.hdl.diff.marg
# plot.sys.diff.marg
# plot.score.diff.marg
# plot.hba1c.diff.marg
# plot.alt.diff.marg
# plot.egfr.diff.marg
# plot.drugline.diff.marg
# plot.category.diff.marg
# plot.ncurrtx.diff.marg
# plot.malesex.diff.marg
# dev.off()



########
# Average patient
average.patient <- cbind(
  drugclass = "SGLT2",
  egfr_ckdepi = mean(dataset.dev$egfr_ckdepi, na.rm = TRUE),
  hba1cmonth = 12,
  prealt = mean(dataset.dev$prealt, na.rm = TRUE),
  prehba1cmmol = mean(dataset.dev$prehba1cmmol, na.rm = TRUE),
  score.excl.mi = mean(dataset.dev$score.excl.mi, na.rm = TRUE),
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
         score.excl.mi = as.numeric(score.excl.mi),
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


# [1] "drugclass"     "Category"      "drugline"      "egfr_ckdepi"
# [5] "hba1cmonth"    "malesex"       "ncurrtx"       "prehba1cmmol"
# [9] "score.excl.mi"

## egfr_ckdepi

# plot treatment effect + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["egfr_ckdepi"]], 
                                                  variable = "egfr_ckdepi", xtitle = "eGFR", 
                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## prealt

# plot treatment effect + histogram marginal
# plot.alt.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prealt"]], 
#                                                   variable = "prealt", xtitle = "ALT", 
#                                                   thinning = 5, k = 4, ymin = -8, ymax = 8)

## prehba1cmmol

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prehba1cmmol"]], 
                                                 variable = "prehba1cmmol", xtitle = "HbA1c", 
                                                 thinning = 5, k = 4, ymin = -8, ymax = 8)

## score.excl.mi

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["score.excl.mi"]], 
                                                   variable = "score.excl.mi", xtitle = "CVD Score", 
                                                   thinning = 5, k = 4, ymin = -8, ymax = 8)
## agetx

# plot treatment effect + histogram marginal
# plot.age.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["agetx"]], 
#                                                    variable = "agetx", xtitle = "Age", 
#                                                    thinning = 5, k = 4, ymin = -8, ymax = 8)
## prehdl

# plot treatment effect + histogram marginal
# plot.hdl.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prehdl"]], 
#                                                  variable = "prehdl", xtitle = "HDL", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## prebmi

# plot treatment effect + histogram marginal
# plot.bmi.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prebmi"]], 
#                                                  variable = "prebmi", xtitle = "BMI", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## prebil

# plot treatment effect + histogram marginal
# plot.bil.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prebil"]], 
#                                                  variable = "prebil", xtitle = "BIL", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## preplatelets

# plot treatment effect + histogram marginal
# plot.platelets.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["preplatelets"]], 
#                                                  variable = "preplatelets", xtitle = "Platelets", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## t2dmduration

# plot treatment effect + histogram marginal
# plot.t2dmduration.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["t2dmduration"]], 
#                                                        variable = "t2dmduration", xtitle = "t2dmduration", 
#                                                        thinning = 5, k = 4, ymin = -8, ymax = 8)

## prealb

# plot treatment effect + histogram marginal
# plot.alb.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["prealb"]], 
#                                                           variable = "prealb", xtitle = "ALB", 
#                                                           thinning = 5, k = 4, ymin = -8, ymax = 8)

## presys

# plot treatment effect + histogram marginal
# plot.sys.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["presys"]], 
#                                                  variable = "presys", xtitle = "SYS", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## preast

# plot treatment effect + histogram marginal
# plot.ast.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["preast"]], 
#                                                  variable = "preast", xtitle = "AST", 
#                                                  thinning = 5, k = 4, ymin = -8, ymax = 8)

## drugline

# plot treatment effect + histogram marginal
plot.drugline.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["drugline"]], 
                                                 variable = "drugline", xtitle = "Drugline", 
                                                 thinning = 5, ymin = -8, ymax = 8)

## Category Smoker

# plot treatment effect + histogram marginal
plot.category.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["Category"]], 
                                                      variable = "Category", xtitle = "Smoking", 
                                                      thinning = 5, ymin = -8, ymax = 8)

## ncurrtx

# plot treatment effect + histogram marginal
plot.ncurrtx.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["ncurrtx"]], 
                                                      variable = "ncurrtx", xtitle = "ncurrtx", 
                                                      thinning = 5, ymin = -8, ymax = 8)

## malesex

# plot treatment effect + histogram marginal
plot.malesex.diff.marg <- plot_diff_treatment_effect(effects = effects_summary_patient[["malesex"]], 
                                                     variable = "malesex", xtitle = "Sex", 
                                                     thinning = 5, ymin = -8, ymax = 8)



#### PDF with all the plots


pdf(file = "Plots/6.4.patient_differential_treatment_effect.pdf")
# plot.age.diff.marg
# plot.ast.diff.marg
# plot.alb.diff.marg
# plot.t2dmduration.diff.marg
# plot.platelets.diff.marg
# plot.bil.diff.marg
# plot.bmi.diff.marg
# plot.hdl.diff.marg
# plot.sys.diff.marg
plot.score.diff.marg
plot.hba1c.diff.marg
# plot.alt.diff.marg
plot.egfr.diff.marg
plot.drugline.diff.marg
plot.category.diff.marg
plot.ncurrtx.diff.marg
plot.malesex.diff.marg

# grid plot of all plots
cowplot::plot_grid(
  # plot.age.diff.marg,
  # plot.ast.diff.marg,
  # plot.alb.diff.marg,
  # plot.t2dmduration.diff.marg,
  # plot.platelets.diff.marg,
  # plot.bil.diff.marg,
  # plot.bmi.diff.marg,
  # plot.hdl.diff.marg,
  # plot.sys.diff.marg,
  plot.score.diff.marg,
  plot.hba1c.diff.marg,
  # plot.alt.diff.marg,
  plot.egfr.diff.marg,
  plot.drugline.diff.marg,
  plot.category.diff.marg,
  plot.ncurrtx.diff.marg,
  plot.malesex.diff.marg
)
dev.off()



# Patients with complete data
# 9433, 9481, 9517, 9601, 9649, 9241, 9301


# patient 39, MRC_plots.R
# patient final hba1c = dataset.dev[9241,"posthba1c_final"] = 54
specific.patient <- cbind(
  drugclass = "SGLT2",
  egfr_ckdepi = as.numeric(dataset.dev[9433, "egfr_ckdepi"]),
  hba1cmonth = 12,
  # prealt = as.numeric(dataset.dev[9433, "prealt"]),
  prehba1cmmol = as.numeric(dataset.dev[9433, "prehba1cmmol"]),
  score.excl.mi = as.numeric(dataset.dev[9433, "score.excl.mi"]),
  Category = "Non-smoker",
  drugline = "2",
  ncurrtx = "1",
  yrdrugstart = 2016,
  # agetx = as.numeric(dataset.dev[9433, "agetx"]),
  malesex = "0"
  # prehdl = as.numeric(dataset.dev[9433, "prehdl"]),
  # prebmi = as.numeric(dataset.dev[9433, "prebmi"]),
  # prebil = as.numeric(dataset.dev[9433, "prebil"]),
  # preplatelets = as.numeric(dataset.dev[9433, "preplatelets"]),
  # t2dmduration = as.numeric(dataset.dev[9433, "t2dmduration"]),
  # prealb = as.numeric(dataset.dev[9433, "prealb"]),
  # presys = as.numeric(dataset.dev[9433, "presys"]),
  # preast = as.numeric(dataset.dev[9433, "preast"])
) %>%
  as.data.frame() %>%
  mutate(drugclass = factor(drugclass, levels = levels(dataset.dev$drugclass)),
         egfr_ckdepi = as.numeric(egfr_ckdepi),
         hba1cmonth = as.numeric(hba1cmonth),
         # prealt = as.numeric(prealt),
         prehba1cmmol = as.numeric(prehba1cmmol),
         score.excl.mi = as.numeric(score.excl.mi),
         Category = factor(Category, levels = levels(dataset.dev$Category)),
         drugline = factor(drugline, levels = levels(dataset.dev$drugline)),
         ncurrtx = factor(ncurrtx, levels = levels(dataset.dev$ncurrtx)),
         yrdrugstart = as.numeric(yrdrugstart),
         # agetx = as.numeric(agetx),
         malesex = factor(malesex, levels = levels(dataset.dev$malesex)),
         # prehdl = as.numeric(prehdl),
         # prebmi = as.numeric(prebmi),
         # prebil = as.numeric(prebil),
         # preplatelets = as.numeric(preplatelets),
         # t2dmduration = as.numeric(t2dmduration),
         # prealb = as.numeric(prealb),
         # presys = as.numeric(presys),
         # preast = as.numeric(preast)
  )



################
## Different treatment response split between therapies

response_summary_patient <- diff_treatment_response(bart_model_final, specific.patient, 25)


## egfr_ckdepi

# plot treatment response + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["egfr_ckdepi"]], 
                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                  variable = "egfr_ckdepi", xtitle = "eGFR",
                                                  ymin = -25, ymax = -5)

## prealt

# plot treatment effect + histogram marginal
# plot.alt.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prealt"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "prealt", xtitle = "ALT",
#                                                  ymin = -25, ymax = -5)

## prehba1cmmol

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prehba1cmmol"]], 
                                                     pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                   variable = "prehba1cmmol", xtitle = "HbA1c")

## score.excl.mi

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["score.excl.mi"]], 
                                                     pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                   variable = "score.excl.mi", xtitle = "CVD Score",
                                                   ymin = -25, ymax = -5)

## agetx

# plot treatment effect + histogram marginal
# plot.age.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["agetx"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "agetx", xtitle = "Age",
#                                                  ymin = -25, ymax = -5)

## prehdl

# plot treatment effect + histogram marginal
# plot.hdl.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prehdl"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "prehdl", xtitle = "HDL",
#                                                  ymin = -25, ymax = -5)

## prebmi

# plot treatment effect + histogram marginal
# plot.bmi.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prebmi"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "prebmi", xtitle = "BMI",
#                                                  ymin = -25, ymax = -5)

## prebil

# plot treatment effect + histogram marginal
# plot.bil.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prebil"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "prebil", xtitle = "BIL",
#                                                  ymin = -25, ymax = -5)

## preplatelets

# plot treatment effect + histogram marginal
# plot.platelets.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["preplatelets"]],
#                                                          pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]), 
#                                                        variable = "preplatelets", xtitle = "Platelets",
#                                                        ymin = -25, ymax = -5)

## t2dmduration

# plot treatment effect + histogram marginal
# plot.t2dmduration.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["t2dmduration"]], 
#                                                             pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                           variable = "t2dmduration", xtitle = "t2dmduration",
#                                                           ymin = -25, ymax = -5)

## prealb

# plot treatment effect + histogram marginal
# plot.alb.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prealb"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "prealb", xtitle = "ALB",
#                                                  ymin = -25, ymax = -5)

## presys

# plot treatment effect + histogram marginal
# plot.sys.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["presys"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "presys", xtitle = "SYS",
#                                                  ymin = -25, ymax = -5)

## preast

# plot treatment effect + histogram marginal
# plot.ast.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["preast"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                  variable = "preast", xtitle = "AST",
#                                                  ymin = -25, ymax = -5)

## drugline

# plot treatment effect + histogram marginal
plot.drugline.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["drugline"]], 
                                                        pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                      variable = "drugline", xtitle = "Drugline")

## Category Smoker

# plot treatment effect + histogram marginal
plot.category.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["Category"]], 
                                                        pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                      variable = "Category", xtitle = "Smoking")

## ncurrtx

# plot treatment effect + histogram marginal
plot.ncurrtx.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["ncurrtx"]], 
                                                       pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                     variable = "ncurrtx", xtitle = "ncurrtx")

plot.all.female <- patchwork::wrap_plots(
  # plot.age.diff.marg,
  # plot.ast.diff.marg,
  # plot.alb.diff.marg,
  # plot.t2dmduration.diff.marg,
  # plot.platelets.diff.marg,
  # plot.bil.diff.marg,
  # plot.bmi.diff.marg,
  # plot.hdl.diff.marg,
  # plot.sys.diff.marg,
  plot.score.diff.marg,
  # plot.alt.diff.marg,
  plot.egfr.diff.marg,
  plot.hba1c.diff.marg,
  plot.drugline.diff.marg,
  plot.category.diff.marg,
  plot.ncurrtx.diff.marg
)


# patient 39, MRC_plots.R
# patient final hba1c = dataset.dev[9241,"posthba1c_final"] = 54
specific.patient <- cbind(
  drugclass = "SGLT2",
  egfr_ckdepi = as.numeric(dataset.dev[9433, "egfr_ckdepi"]),
  hba1cmonth = 12,
  # prealt = as.numeric(dataset.dev[9433, "prealt"]),
  prehba1cmmol = as.numeric(dataset.dev[9433, "prehba1cmmol"]),
  score.excl.mi = as.numeric(dataset.dev[9433, "score.excl.mi"]),
  Category = "Non-smoker",
  drugline = "2",
  ncurrtx = "1",
  yrdrugstart = 2016,
  # agetx = as.numeric(dataset.dev[9433, "agetx"]),
  malesex = "1"
  # prehdl = as.numeric(dataset.dev[9433, "prehdl"]),
  # prebmi = as.numeric(dataset.dev[9433, "prebmi"]),
  # prebil = as.numeric(dataset.dev[9433, "prebil"]),
  # preplatelets = as.numeric(dataset.dev[9433, "preplatelets"]),
  # t2dmduration = as.numeric(dataset.dev[9433, "t2dmduration"]),
  # prealb = as.numeric(dataset.dev[9433, "prealb"]),
  # presys = as.numeric(dataset.dev[9433, "presys"]),
  # preast = as.numeric(dataset.dev[9433, "preast"])
) %>%
  as.data.frame() %>%
  mutate(drugclass = factor(drugclass, levels = levels(dataset.dev$drugclass)),
         egfr_ckdepi = as.numeric(egfr_ckdepi),
         hba1cmonth = as.numeric(hba1cmonth),
         # prealt = as.numeric(prealt),
         prehba1cmmol = as.numeric(prehba1cmmol),
         score.excl.mi = as.numeric(score.excl.mi),
         Category = factor(Category, levels = levels(dataset.dev$Category)),
         drugline = factor(drugline, levels = levels(dataset.dev$drugline)),
         ncurrtx = factor(ncurrtx, levels = levels(dataset.dev$ncurrtx)),
         yrdrugstart = as.numeric(yrdrugstart),
         # agetx = as.numeric(agetx),
         malesex = factor(malesex, levels = levels(dataset.dev$malesex)),
         # prehdl = as.numeric(prehdl),
         # prebmi = as.numeric(prebmi),
         # prebil = as.numeric(prebil),
         # preplatelets = as.numeric(preplatelets),
         # t2dmduration = as.numeric(t2dmduration),
         # prealb = as.numeric(prealb),
         # presys = as.numeric(presys),
         # preast = as.numeric(preast)
  )



################
## Different treatment response split between therapies

response_summary_patient <- diff_treatment_response(bart_model_final, specific.patient, 25)


## egfr_ckdepi

# plot treatment response + histogram marginal
plot.egfr.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["egfr_ckdepi"]], 
                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                    variable = "egfr_ckdepi", xtitle = "eGFR",
                                                    ymin = -25, ymax = -5)

## prealt

# plot treatment effect + histogram marginal
# plot.alt.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prealt"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "prealt", xtitle = "ALT",
#                                                    ymin = -25, ymax = -5)

## prehba1cmmol

# plot treatment effect + histogram marginal
plot.hba1c.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prehba1cmmol"]], 
                                                     pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                     variable = "prehba1cmmol", xtitle = "HbA1c")

## score.excl.mi

# plot treatment effect + histogram marginal
plot.score.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["score.excl.mi"]], 
                                                     pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                     variable = "score.excl.mi", xtitle = "CVD Score",
                                                     ymin = -25, ymax = -5)

## agetx

# plot treatment effect + histogram marginal
# plot.age.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["agetx"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "agetx", xtitle = "Age",
#                                                    ymin = -25, ymax = -5)

## prehdl

# plot treatment effect + histogram marginal
# plot.hdl.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prehdl"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "prehdl", xtitle = "HDL",
#                                                    ymin = -25, ymax = -5)

## prebmi

# plot treatment effect + histogram marginal
# plot.bmi.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prebmi"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "prebmi", xtitle = "BMI",
#                                                    ymin = -25, ymax = -5)

## prebil

# plot treatment effect + histogram marginal
# plot.bil.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prebil"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "prebil", xtitle = "BIL",
#                                                    ymin = -25, ymax = -5)

## preplatelets

# plot treatment effect + histogram marginal
# plot.platelets.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["preplatelets"]],
#                                                          pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]), 
#                                                          variable = "preplatelets", xtitle = "Platelets",
#                                                          ymin = -25, ymax = -5)

## t2dmduration

# plot treatment effect + histogram marginal
# plot.t2dmduration.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["t2dmduration"]], 
#                                                             pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                             variable = "t2dmduration", xtitle = "t2dmduration",
#                                                             ymin = -25, ymax = -5)

## prealb

# plot treatment effect + histogram marginal
# plot.alb.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["prealb"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "prealb", xtitle = "ALB",
#                                                    ymin = -25, ymax = -5)

## presys

# plot treatment effect + histogram marginal
# plot.sys.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["presys"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "presys", xtitle = "SYS",
#                                                    ymin = -25, ymax = -5)

## preast

# plot treatment effect + histogram marginal
# plot.ast.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["preast"]], 
#                                                    pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
#                                                    variable = "preast", xtitle = "AST",
#                                                    ymin = -25, ymax = -5)

## drugline

# plot treatment effect + histogram marginal
plot.drugline.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["drugline"]], 
                                                        pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                        variable = "drugline", xtitle = "Drugline")

## Category Smoker

# plot treatment effect + histogram marginal
plot.category.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["Category"]], 
                                                        pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                        variable = "Category", xtitle = "Smoking")

## ncurrtx

# plot treatment effect + histogram marginal
plot.ncurrtx.diff.marg <- plot_diff_treatment_response(response = response_summary_patient[["ncurrtx"]], 
                                                       pre_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                       variable = "ncurrtx", xtitle = "ncurrtx")


plot.all.male <- patchwork::wrap_plots(
  # plot.age.diff.marg,
  # plot.ast.diff.marg,
  # plot.alb.diff.marg,
  # plot.t2dmduration.diff.marg,
  # plot.platelets.diff.marg,
  # plot.bil.diff.marg,
  # plot.bmi.diff.marg,
  # plot.hdl.diff.marg,
  # plot.sys.diff.marg,
  plot.score.diff.marg,
  # plot.alt.diff.marg,
  plot.egfr.diff.marg,
  plot.hba1c.diff.marg,
  plot.drugline.diff.marg,
  plot.category.diff.marg,
  plot.ncurrtx.diff.marg
)


pdf(file = "Plots/6.4.patient_differential_treatment_response.pdf")
# plot.age.diff.marg
# plot.ast.diff.marg
# plot.alb.diff.marg
# plot.t2dmduration.diff.marg
# plot.platelets.diff.marg
# plot.bil.diff.marg
# plot.bmi.diff.marg
# plot.hdl.diff.marg
# plot.sys.diff.marg
plot.score.diff.marg
plot.hba1c.diff.marg
# plot.alt.diff.marg
plot.egfr.diff.marg
plot.drugline.diff.marg
plot.category.diff.marg
plot.ncurrtx.diff.marg

# grid plot of all plots
plot.all.male
plot.all.female
dev.off()






####################
## Description:
##  - In this file we exclude all GLP1 patients before 2013
##      in order to investigate how this affects treatment effect.
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
dir.create(paste0(output_path, "/Final_model"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/7.1.Sensitivity"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment"))


## make directory for outputs
dir.create("Plots")



# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))

#### Removal of patients before 2013

final.dev <- final.dev %>%
  filter(yrdrugstart > 2012)

final.val <- final.val %>%
  filter(yrdrugstart > 2012)


###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################

source("0.1.slade_functions.R")


###############################################################################
###############################################################################
############################### Model Fitting #################################
###############################################################################
###############################################################################
##
## Initially we start with 9866 individuals
##

dataset.dev <- final.dev


########
### Building the initial outcome model with var selection
########



## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(

  bart_model <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/bart_model.rds"))

  , silent = TRUE)) == "try-error") {

  bart_model <- bartMachine::bartMachine(X = dataset.dev %>%
                                           select(-patid,
                                                  -pateddrug,
                                                  -bothdrugs,
                                                  -posthba1c_final,
                                                  -sglt2subtype,
                                                  -glp1subtype),
                                         y = dataset.dev[,"posthba1c_final"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)

  saveRDS(bart_model, paste0(output_path, "/Final_model/7.1.Sensitivity/bart_model.rds"))

}


## Perform variable selection from the outcome model

if (class(try(

  bart_var_selection <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/bart_var_selection.rds"))

  , silent = TRUE)) == "try-error") {

  bart_var_selection <- bartMachine::var_selection_by_permute(bart_machine = bart_model,
                                                              num_reps_for_avg = 15,
                                                              num_permute_samples = 150,
                                                              num_trees_for_permute = 50)

  saveRDS(bart_var_selection, paste0(output_path, "/Final_model/7.1.Sensitivity/bart_var_selection.rds"))

}

# Variables selected:
#
# [1] "prehba1cmmol"       "prealt"             "drugclass_SGLT2"
# [4] "hba1cmonth"         "yrdrugstart"        "egfr_ckdepi"
# [7] "drugclass_GLP1"     "score"              "drugline_2"
# [10] "ncurrtx_3"          "drugline_3"         "ncurrtx_2"
# [13] "drugline_4"         "Category_Ex-smoker" "drugline_5"
# [16] "ncurrtx_1"
# male_sex not included this time


if (class(try(

  bart_var_selection_cv <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/bart_var_selection_cv.rds"))

  , silent = TRUE)) == "try-error") {

  bart_var_selection_cv <- bartMachine::var_selection_by_permute_cv(bart_machine = bart_model,
                                                                    k_folds = 15,
                                                                    num_reps_for_avg = 15,
                                                                    num_permute_samples = 150,
                                                                    num_trees_for_permute = 50,
                                                                    num_trees_pred_cv = 50)

  saveRDS(bart_var_selection_cv, paste0(output_path, "/Final_model/7.1.Sensitivity/bart_var_selection_cv.rds"))

}

# Variables selected:
#
# [1] "drugclass_GLP1"  "drugclass_SGLT2" "drugline_2"      "drugline_3"
# [5] "drugline_4"      "drugline_5"      "egfr_ckdepi"     "hba1cmonth"
# [9] "ncurrtx_1"       "ncurrtx_2"       "ncurrtx_3"       "prealt"
# [13] "prebil"          "prehba1cmmol"    "score"           "yrdrugstart"
# Category  not included, prebil extra



# #### GRF variable selection
# current.na.action <- options("na.action")
# options(na.action = "na.pass")
# 
# dataset_model.matrix <- model.matrix(~., final.dev %>%
#                                        select(-patid,
#                                               -pateddrug,
#                                               -bothdrugs)) %>%
#   as.data.frame() %>%
#   select(-`(Intercept)`) %>%
#   mutate(drugclass = drugclassSGLT2) %>%
#   select(-drugclassSGLT2)
# 
# 
# 
# ## Fit a propensity model with all the variables
# if (class(try(
# 
#   prop_model <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/prop_model.rds"))
# 
#   , silent = TRUE)) == "try-error") {
# 
#   # GLP1 is considered target so all the probabilities should be 1-prob
#   prop_model <- bartMachine::bartMachine(X = dataset.dev %>%
#                                            select(-patid,
#                                                   -pateddrug,
#                                                   -bothdrugs,
#                                                   -posthba1c_final,
#                                                   -drugclass,
#                                                   -yrdrugstart,
#                                                   -sglt2subtype,
#                                                   -glp1subtype),
#                                          y = dataset.dev[,"drugclass"],
#                                          use_missing_data = TRUE,
#                                          impute_missingness_with_rf_impute = FALSE,
#                                          impute_missingness_with_x_j_bar_for_lm = TRUE,
#                                          num_trees = 200,
#                                          num_burn_in = 3000,
#                                          num_iterations_after_burn_in = 1000,
#                                          serialize = TRUE)
# 
#   saveRDS(prop_model, paste0(output_path, "/Final_model/7.1.Sensitivity/prop_model.rds"))
# }
# 
# # Perform variables selection of important variables
# 
# if (class(try(
# 
#   prop_var_selection <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/prop_var_selection.rds"))
# 
#   , silent = TRUE)) == "try-error") {
# 
#   prop_var_selection <- bartMachine::var_selection_by_permute(bart_machine = prop_model,
#                                                               num_reps_for_avg = 100,
#                                                               num_permute_samples = 1000,
#                                                               num_trees_for_permute = 50)
# 
#   saveRDS(prop_var_selection, paste0(output_path, "/Final_model/7.1.Sensitivity/prop_var_selection.rds"))
# 
# 
# }
# # Variables selected:
# #
# # [1] "prebmi"              "t2dmduration"        "prehba1cmmol"
# # [4] "preweight"           "drugline_2"          "drugline_5"
# # [7] "ncurrtx_3"           "drugline_4"          "drugline_3"
# # [10] "Category_Non-smoker"
# 
# 
# 
# # Fit final propensity model using selected variables
# 
# if (class(try(
# 
#   prop_model_final <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/prop_model_final.rds"))
# 
#   , silent = TRUE)) == "try-error") {
# 
#   prop_model_final <- bartMachine::bartMachine(X = dataset.dev %>%
#                                                  select(prebmi,
#                                                         t2dmduration,
#                                                         # prealb,
#                                                         # egfr_ckdepi,
#                                                         drugline,
#                                                         prehba1cmmol,
#                                                         ncurrtx,
#                                                         # score,
#                                                         # Category
#                                                         ),
#                                                y = dataset.dev[,"drugclass"],
#                                                use_missing_data = TRUE,
#                                                impute_missingness_with_rf_impute = FALSE,
#                                                impute_missingness_with_x_j_bar_for_lm = TRUE,
#                                                num_trees = 200,
#                                                num_burn_in = 3000,
#                                                num_iterations_after_burn_in = 1000,
#                                                serialize = TRUE)
#   saveRDS(prop_model_final, paste0(output_path, "/Final_model/7.1.Sensitivity/prop_model_final.rds"))
# 
# }
# 
# 
# grf_model <- grf::causal_forest(X = dataset_model.matrix %>%
#                                   select(-posthba1c_final,
#                                          -drugclass,
#                                          -preweight,
#                                          -height),
#                                 Y = dataset_model.matrix[, "posthba1c_final"],
#                                 W = dataset_model.matrix[, "drugclass"],
#                                 W.hat = prop_model_final$p_hat_train)
# 
# 
# grf_var_selection <- grf::variable_importance(grf_model) %>%
#   t()
# 
# colnames(grf_model$X.orig)[which(grf_var_selection > 0.01)]
# 
# # malesex1 <- 0.20329320
# # egfr_ckdepi <- 0.17543374
# # score <- 0.08494679
# # agetx <- 0.07947302
# # prealt <- 0.07819542
# # prehba1cmmol <- 0.05829137
# # prehdl <- 0.04667077
# # preplatelets <- 0.04575637
# # prebmi <- 0.03981382
# # hba1cmonth <- 0.02812800
# # t2dmduration <- 0.02589190
# # presys <- 0.02513101
# # prebil <- 0.02227799
# # yrdrugstart <- 0.01916421
# # prealb <- 0.01903258
# # preast <- 0.01068304
# 
# options(na.action = current.na.action)


## BART model with all the selected variables


if (class(try(
  
  bart_model_final <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/bart_model_final.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model_final <- bartMachine::bartMachine(X = dataset.dev %>%
                                                 select(drugclass,
                                                        # below is vars from both variable selections
                                                        egfr_ckdepi,
                                                        hba1cmonth,
                                                        prealt,
                                                        prehba1cmmol,
                                                        score,
                                                        prebil,
                                                        # below is vars from BART variable selection
                                                        Category, # should be removed, not important for this model
                                                        drugline,
                                                        ncurrtx,
                                                        yrdrugstart,
                                                        # below is the vars from grf variable selection
                                                        agetx,
                                                        malesex,
                                                        prehdl,
                                                        prebmi,
                                                        preplatelets,
                                                        t2dmduration,
                                                        prealb,
                                                        presys,
                                                        preast
                                                 ),
                                               y = dataset.dev[,"posthba1c_final"],
                                               use_missing_data = TRUE,
                                               impute_missingness_with_rf_impute = FALSE,
                                               impute_missingness_with_x_j_bar_for_lm = TRUE,
                                               num_trees = 200,
                                               num_burn_in = 3000,
                                               num_iterations_after_burn_in = 1000,
                                               serialize = TRUE)
  
  saveRDS(bart_model_final, paste0(output_path, "/Final_model/7.1.Sensitivity/bart_model_final.rds"))
  
}



########
### Validation of model
########

# Dev

data_dev <- dataset.dev %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))



## Get posteriors
if (class(try(
  
  posteriors_dev <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/posteriors_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_final, data_dev %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_dev, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/posteriors_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  cred_pred_dev <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_dev <- calc_resid(data_dev, posteriors_dev)
  
  saveRDS(cred_pred_dev, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/cred_pred_dev.rds"))
  
}


# Val
data_val <- final.val %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


### Get posteriors
if (class(try(
  
  posteriors_val <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/posteriors_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_val <- bartMachine::bart_machine_get_posterior(bart_model_final, data_val %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_val, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/posteriors_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  cred_pred_val <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_val <- calc_resid(data_val, posteriors_val)
  
  saveRDS(cred_pred_val, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/cred_pred_val.rds"))
  
}



# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/assessment.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_dev, posteriors_dev)
  
  assessment_values_val <- calc_assessment(data_val, posteriors_val)
  
  assessment <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/assessment.rds"))
  
}


assessment <- assessment %>%
  cbind(Model = "Model") %>%
  as.data.frame() %>%
  mutate(`5%` = as.numeric(`5%`),
         `50%` = as.numeric(`50%`),
         `95%` = as.numeric(`95%`),
         Model = factor(Model)
  )

plot_assessment <- assessment %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(y = Model, xmin = `5%`, xmax = `95%`, colour = Model), width = 0.2) +
  geom_point(aes(x = `50%`, y = Model, shape = Dataset), size = 2, colour = "black") +
  facet_wrap(~statistic, ncol = 1, scales = "free") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(colour = "none")





plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Residuals of Model 4 (Post 2012)")




########
### Validation of effect
########

# calculate effects
if (class(try(
  
  effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev <- calc_effect_summary(bart_model_final, data_dev)
  
  saveRDS(effects_summary_dev, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/effects_summary_dev.rds"))
  
}


# calculate effects
if (class(try(
  
  effects_summary_val <- readRDS(paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_val <- calc_effect_summary(bart_model_final, data_val)
  
  saveRDS(effects_summary_val, paste0(output_path, "/Final_model/7.1.Sensitivity/Assessment/effects_summary_val.rds"))
  
}


## plot effects validation

predicted_observed_dev <- data_dev %>%
  cbind(hba1c_diff = effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

predicted_observed_val <- data_val %>%
  cbind(hba1c_diff = effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_model_final)




## plot histogram of effect

plot_effect_1 <- hist_plot(effects_summary_dev, "", -15, 20)

plot_effect_2 <- hist_plot(effects_summary_val, "", -15, 20)




#### PDF with all the plots


pdf(file = "Plots/7.1.model4_plots.pdf")
plot_residuals
plot_assessment
plot_effects_validation
cowplot::plot_grid(plot_effect_1, plot_effect_2, ncol = 2, nrow = 1, labels = c("A", "B"))
dev.off()

pdf(file = "Plots/7.1.model4_partial_dependence.pdf")
features <- rep(0, length(bart_model_final$training_data_features)) %>%
  as.data.frame() %>%
  t()
colnames(features) <- bart_model_final$training_data_features
features <- features %>%
  as.data.frame() %>%
  select(!contains("_"))
for (i in colnames(features)) {
  pd_plot(bart_model_final, i)
}
dev.off()

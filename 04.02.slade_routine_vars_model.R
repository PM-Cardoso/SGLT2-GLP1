####################
## Description:
##  - In this file we fit a BART model between HbA1c outcome and routine
##      clinical features, using propensity scores as way to deal with confounders
##      We also evaluate goodness of fit + treatment effects validation
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
dir.create(paste0(output_path, "/Final_model/model_2"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/model_2/Assessment"))

## make directory for outputs
dir.create("Plots")



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))

# load all data for range of variable values; name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")


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
### Building the initial prop model with var selection
########

if (class(try(
  
  prop_model <- readRDS(paste0(output_path, "/Final_model/model_2/prop_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # GLP1 is considered target so all the probabilities should be 1-prob
  prop_model <- bartMachine::bartMachine(X = dataset.dev %>%
                                           select(t2dmduration,
                                                  drugline,
                                                  ncurrtx,
                                                  hba1cmonth
                                           ),
                                         y = dataset.dev[,"drugclass"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)
  
  saveRDS(prop_model, paste0(output_path, "/Final_model/model_2/prop_model.rds"))
}

# ggplot() +
#   geom_histogram(aes(x = 1-prop_model$p_hat_train)) +
#   ggtitle("Model 2 propensity scores")
# 


########
### Building the initial outcome model with var selection
########

# Add propensity score as a possible variable
dataset.dev <- dataset.dev %>%
  cbind(prop_score = prop_model$p_hat_train)


## Fit initial model using all the available variables to estimate HbA1c outcome

if (class(try(
  
  bart_model <- readRDS(paste0(output_path, "/Final_model/model_2/bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model <- bartMachine::bartMachine(X = dataset.dev %>%
                                           select(-patid,
                                                  -pateddrug,
                                                  -bothdrugs,
                                                  -posthba1c_final,
                                                  -sglt2subtype,
                                                  -glp1subtype,
                                                  -t2dmduration,
                                                  -drugline,
                                                  -ncurrtx,
                                                  -hba1cmonth,
                                           ),
                                         y = dataset.dev[,"posthba1c_final"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)
  
  saveRDS(bart_model, paste0(output_path, "/Final_model/model_2/bart_model.rds"))
  
}


## Perform variable selection from the outcome model

if (class(try(
  
  bart_var_selection <- readRDS(paste0(output_path, "/Final_model/model_2/bart_var_selection.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_var_selection <- bartMachine::var_selection_by_permute(bart_machine = bart_model)
  
  saveRDS(bart_var_selection, paste0(output_path, "/Final_model/model_2/bart_var_selection.rds"))
  
}


# Variables selected:
# [1] "prehba1cmmol"       "prop_score"         "prealt"
# [4] "egfr_ckdepi"        "yrdrugstart"        "drugclass_SGLT2"
# [7] "drugclass_GLP1"     "malesex_0"          "malesex_1"
# [10] "Category_Ex-smoker"



if (class(try(
  
  bart_var_selection_cv <- readRDS(paste0(output_path, "/Final_model/model_2/bart_var_selection_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_var_selection_cv <- bartMachine::var_selection_by_permute_cv(bart_machine = bart_model)
  
  saveRDS(bart_var_selection_cv, paste0(output_path, "/Final_model/model_2/bart_var_selection_cv.rds"))
  
}


# Variables selected:
# [1] "Category_Ex-smoker" "drugclass_GLP1"     "drugclass_SGLT2"
# [4] "egfr_ckdepi"        "malesex_0"          "malesex_1"
# [7] "prealt"             "prehba1cmmol"       "prop_score"
# [10] "score"



###


## BART model with all the selected variables

if (class(try(

  bart_model_final <- readRDS(paste0(output_path, "/Final_model/model_2/bart_model_final.rds"))

  , silent = TRUE)) == "try-error") {

  bart_model_final <- bartMachine::bartMachine(X = dataset.dev %>%
                                                 select(
                                                   Category,
                                                   drugclass,
                                                   egfr_ckdepi,
                                                   malesex,
                                                   prealt,
                                                   prehba1cmmol,
                                                   prop_score,
                                                   score
                                                 ),
                                               y = dataset.dev[,"posthba1c_final"],
                                               use_missing_data = TRUE,
                                               impute_missingness_with_rf_impute = FALSE,
                                               impute_missingness_with_x_j_bar_for_lm = TRUE,
                                               num_trees = 200,
                                               num_burn_in = 3000,
                                               num_iterations_after_burn_in = 1000,
                                               serialize = TRUE)

  saveRDS(bart_model_final, paste0(output_path, "/Final_model/model_2/bart_model_final.rds"))

}



########
### Validation of model
########

# Dev

data_dev <- dataset.dev %>%
  select(c(patid, pateddrug, posthba1c_final,
           unique(c(colnames(bart_model_final$X),
                    colnames(prop_model$X)))))



## Get posteriors
if (class(try(

  posteriors_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/posteriors_dev.rds"))

  , silent = TRUE)) == "try-error") {

  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_final, data_dev %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_dev, paste0(output_path, "/Final_model/model_2/Assessment/posteriors_dev.rds"))


}

### residuals calculation
if (class(try(

  cred_pred_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/cred_pred_dev.rds"))

  , silent = TRUE)) == "try-error") {

  cred_pred_dev <- calc_resid(data_dev, posteriors_dev, "posthba1c_final")

  saveRDS(cred_pred_dev, paste0(output_path, "/Final_model/model_2/Assessment/cred_pred_dev.rds"))

}


# Val
data_val <- final.val %>%
  select(
    c(patid, pateddrug, posthba1c_final, colnames(prop_model$X))
    )

# calculate prop score
if (class(try(

  prop_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/prop_val.rds"))

  , silent = TRUE)) == "try-error") {

  prop_val <- predict(prop_model, data_val %>%
                        select(
                          colnames(prop_model$X)
                        ))

  saveRDS(prop_val, paste0(output_path, "/Final_model/model_2/Assessment/prop_val.rds"))

}


data_val <- final.val %>%
  cbind(prop_score = prop_val) %>%
  select(c(patid, pateddrug, posthba1c_final,
           unique(c(colnames(bart_model_final$X),
                    colnames(prop_model$X)))))


### Get posteriors
if (class(try(

  posteriors_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/posteriors_val.rds"))

  , silent = TRUE)) == "try-error") {

  posteriors_val <- bartMachine::bart_machine_get_posterior(bart_model_final, data_val %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_val, paste0(output_path, "/Final_model/model_2/Assessment/posteriors_val.rds"))


}

### residuals calculation
if (class(try(

  cred_pred_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/cred_pred_val.rds"))

  , silent = TRUE)) == "try-error") {

  cred_pred_val <- calc_resid(data_val, posteriors_val, "posthba1c_final")

  saveRDS(cred_pred_val, paste0(output_path, "/Final_model/model_2/Assessment/cred_pred_val.rds"))

}



# assessment of R2, RSS, RMSE
if (class(try(

  assessment <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/assessment.rds"))

  , silent = TRUE)) == "try-error") {

  assessment_values_dev <- calc_assessment(data_dev, posteriors_dev, "posthba1c_final")

  assessment_values_val <- calc_assessment(data_val, posteriors_val, "posthba1c_final")

  assessment <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )

  saveRDS(assessment, paste0(output_path, "/Final_model/model_2/Assessment/assessment.rds"))

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





plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Residuals of Model 2")




########
### Validation of effect
########

# calculate effects
if (class(try(

  effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/effects_summary_dev.rds"))

  , silent = TRUE)) == "try-error") {

  effects_summary_dev <- calc_effect_summary(bart_model_final, data_dev)

  saveRDS(effects_summary_dev, paste0(output_path, "/Final_model/model_2/Assessment/effects_summary_dev.rds"))

}


# calculate effects
if (class(try(

  effects_summary_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/effects_summary_val.rds"))

  , silent = TRUE)) == "try-error") {

  effects_summary_val <- calc_effect_summary(bart_model_final, data_val)

  saveRDS(effects_summary_val, paste0(output_path, "/Final_model/model_2/Assessment/effects_summary_val.rds"))

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


# plot_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_model_final)




## plot histogram of effect

plot_effect_1 <- hist_plot(effects_summary_dev, "", -15, 20)

plot_effect_2 <- hist_plot(effects_summary_val, "", -15, 20)

effects_summary_dev_male <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 1)

effects_summary_dev_female <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


effects_summary_val_male <- effects_summary_val %>%
  cbind(malesex = final.val$malesex) %>%
  filter(malesex == 1)

effects_summary_val_female <- effects_summary_val %>%
  cbind(malesex = final.val$malesex) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)




# ##############
# Calculate propensity score for validation

# extracting selected variables for individuals in dataset
data.new <- predicted_observed_dev %>%
  select(patid, pateddrug) %>%
  left_join(final.all.extra.vars %>%
              select(patid,
                     pateddrug,
                     drugclass,
                     yrdrugstart,
                     prebmi,
                     t2dmduration,
                     drugline,
                     prehba1cmmol,
                     egfr_ckdepi,
                     ncurrtx,
                     Category), by = c("patid", "pateddrug"))

# fit propensity model with the variables that influence therapy indication
set.seed(123)
prop_model_dev <- bartMachine::bartMachine(X = data.new %>%
                                             select(yrdrugstart,
                                                    prebmi,
                                                    t2dmduration,
                                                    drugline,
                                                    prehba1cmmol,
                                                    egfr_ckdepi,
                                                    ncurrtx,
                                                    Category),
                                           y = data.new[,"drugclass"],
                                           use_missing_data = TRUE,
                                           impute_missingness_with_rf_impute = FALSE,
                                           impute_missingness_with_x_j_bar_for_lm = TRUE,
                                           num_trees = 200,
                                           num_burn_in = 1000,
                                           num_iterations_after_burn_in = 200,
                                           seed = 123)

# extracting selected variables for individuals in dataset
data.new <- predicted_observed_val %>%
  select(patid, pateddrug) %>%
  left_join(final.all.extra.vars %>%
              select(patid,
                     pateddrug,
                     drugclass,
                     yrdrugstart,
                     prebmi,
                     t2dmduration,
                     drugline,
                     prehba1cmmol,
                     egfr_ckdepi,
                     ncurrtx,
                     Category), by = c("patid", "pateddrug"))

# fit propensity model with the variables that influence therapy indication
set.seed(123)
prop_model_val <- bartMachine::bartMachine(X = data.new %>%
                                             select(yrdrugstart,
                                                    prebmi,
                                                    t2dmduration,
                                                    drugline,
                                                    prehba1cmmol,
                                                    egfr_ckdepi,
                                                    ncurrtx,
                                                    Category),
                                           y = data.new[,"drugclass"],
                                           use_missing_data = TRUE,
                                           impute_missingness_with_rf_impute = FALSE,
                                           impute_missingness_with_x_j_bar_for_lm = TRUE,
                                           num_trees = 200,
                                           num_burn_in = 1000,
                                           num_iterations_after_burn_in = 200,
                                           seed = 123)



# Validating ATE
if (class(try(

  ATE_validation_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_validation_dev <- calc_ATE_validation(predicted_observed_dev, "posthba1c_final", prop_model_dev)

  saveRDS(ATE_validation_dev, paste0(output_path, "/Final_model/model_2/Assessment/ATE_validation_dev.rds"))

}

plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


if (class(try(

  ATE_validation_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_validation_val <- calc_ATE_validation(predicted_observed_val, "posthba1c_final", prop_model_val)

  saveRDS(ATE_validation_val, paste0(output_path, "/Final_model/model_2/Assessment/ATE_validation_val.rds"))

}

plot_ATE_val <- ATE_plot(ATE_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

plot_ATE <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")

  ,

  cowplot::plot_grid(plot_ATE_dev, plot_ATE_val, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(

  ATE_matching_validation_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_matching_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_matching_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_dev, "posthba1c_final", prop_model_dev)

  saveRDS(ATE_matching_validation_dev, paste0(output_path, "/Final_model/model_2/Assessment/ATE_matching_validation_dev.rds"))

}

plot_ATE_dev_prop_score <- ATE_plot(ATE_matching_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

if (class(try(

  ATE_matching_validation_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_matching_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_matching_validation_val <- calc_ATE_validation_prop_matching(predicted_observed_val, "posthba1c_final", prop_model_val)

  saveRDS(ATE_matching_validation_val, paste0(output_path, "/Final_model/model_2/Assessment/ATE_matching_validation_val.rds"))

}

plot_ATE_val_prop_score <- ATE_plot(ATE_matching_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

plot_ATE_prop_score_matching <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")

  ,

  cowplot::plot_grid(plot_ATE_dev_prop_score, plot_ATE_val_prop_score, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(

  ATE_weighting_validation_dev <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_weighting_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "posthba1c_final", prop_model_dev)

  saveRDS(ATE_weighting_validation_dev, paste0(output_path, "/Final_model/model_2/Assessment/ATE_weighting_validation_dev.rds"))

}

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

if (class(try(

  ATE_weighting_validation_val <- readRDS(paste0(output_path, "/Final_model/model_2/Assessment/ATE_weighting_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1c_final", prop_model_val)

  saveRDS(ATE_weighting_validation_val, paste0(output_path, "/Final_model/model_2/Assessment/ATE_weighting_validation_val.rds"))

}

plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

plot_ATE_prop_score_weighting <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")

  ,

  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting, plot_ATE_val_prop_score_weighting, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))




#### PDF with all the plots


pdf(file = "Plots/4.2.model2_plots.pdf")
hist(1-prop_model$p_hat_train)
plot_residuals
plot_assessment
# plot_effects_validation
cowplot::plot_grid(plot_effect_1, plot_effect_2, ncol = 2, nrow = 1, labels = c("A", "B"))

cowplot::plot_grid(

  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)

  ,

  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)

  , nrow = 2, ncol = 1, labels = c("A", "B")
)
plot_ATE
plot_ATE_prop_score_matching
plot_ATE_prop_score_weighting
dev.off()


pdf(file = "Plots/4.2.model2_partial_dependence.pdf")
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




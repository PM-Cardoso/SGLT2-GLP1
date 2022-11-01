####################
## Description:
##  - In this file we:
##    - Fit a BART response model to all variables.
##    - Perform variable selection on BART response model and refit model.
##    - Refit BART response model
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

output_path <- "Samples/SGLT2-GLP1/Aurum"

## make directory for outputs
dir.create(output_path)

## male directory for outputs
dir.create(paste0(output_path, "/response_model"))

## male directory for outputs
dir.create(paste0(output_path, "/response_model/assessment"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")



hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train")


# collect propensity score values
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

hba1c.train <- hba1c.train %>%
  left_join(patient_prop_scores, by = c("patid", "pated"))


## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_model <- readRDS(paste0(output_path, "/response_model/bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  set.seed(123)
  bart_model <- bartMachine::bartMachine(X = hba1c.train %>%
                                           select(-patid,
                                                  -pated,
                                                  -posthba1cfinal,
                                                  #due to missingness
                                                  -preast),
                                         y = hba1c.train[,"posthba1cfinal"] %>%
                                           unlist(),
                                         use_missing_data = TRUE,
                                         num_trees = 50,
                                         num_burn_in = 2000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE,
                                         seed = 123)
  
  saveRDS(bart_model, paste0(output_path, "/response_model/bart_model.rds"))
  
}


## Best number of trees by rmse
if (class(try(
  
  bart_model_rmse_trees <- readRDS(paste0(output_path, "/response_model/bart_model_rmse_trees.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/11.05.response_model_rmse_trees.pdf", width = 18, height = 11)
  set.seed(123)
  bart_model_rmse_trees <- bartMachine::rmse_by_num_trees(bart_model)
  dev.off()
  
  saveRDS(bart_model_rmse_trees, paste0(output_path, "/response_model/bart_model_rmse_trees.rds"))
}


## Fit initial model using all the available variables to estimate HbA1c outcome by CV
if (class(try(
  
  bart_model_cv <- readRDS(paste0(output_path, "/response_model/bart_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  set.seed(123)
  bart_model_cv <- bartMachine::bartMachineCV(X = hba1c.train %>%
                                                select(
                                                  -patid,
                                                  -pated,
                                                  -posthba1cfinal,
                                                  #due to missingness
                                                  -preast
                                                  ),
                                              y = hba1c.train[,"posthba1cfinal"] %>%
                                                unlist(),
                                              use_missing_data = TRUE,
                                              num_burn_in = 2000,
                                              num_iterations_after_burn_in = 1000,
                                              serialize = TRUE,
                                              seed = 123)
  
  saveRDS(bart_model_cv, paste0(output_path, "/response_model/bart_model_cv.rds"))
  
}


########
### Variable selection of response model
########

if (class(try(
  
  vs_bart_model <- readRDS(paste0(output_path, "/response_model/vs_bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/11.05.response_model_vs.pdf", width = 18, height = 11)
  # error with cv
  set.seed(123)
  vs_bart_model <- var_selection_by_permute(bart_model_cv)
  dev.off()
  
  ## Variables selected
  # [1] "prehba1c"   "drugclass"  "preegfr"    "drugline"   "hba1cmonth"
  # [6] "agetx"
  
  
  saveRDS(vs_bart_model, paste0(output_path, "/response_model/vs_bart_model.rds"))
  
}

# Cross-validation selection
if (class(try(
  
  vs_bart_model_cv <- readRDS(paste0(output_path, "/response_model/vs_bart_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  set.seed(123)
  vs_bart_model_cv <- var_selection_by_permute_cv(bart_model_cv)
  
  ## Long version of the var selection
  # [1] "agetx"      "drugclass"  "drugline"   "hba1cmonth" "preegfr"
  # [6] "prehba1c"
  
  
  saveRDS(vs_bart_model_cv, paste0(output_path, "/response_model/vs_bart_model_cv.rds"))
  
}


########
### Refit Response model with selected vars
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(

  bart_model_final <- readRDS(paste0(output_path, "/response_model/bart_model_final.rds"))

  , silent = TRUE)) == "try-error") {

  set.seed(123)
  bart_model_final <- bartMachine::bartMachineCV(X = hba1c.train %>%
                                                   select(
                                                     vs_bart_model_cv$important_vars_cv
                                                     ),
                                                 y = hba1c.train[,"posthba1cfinal"] %>%
                                                   unlist(),
                                                 use_missing_data = TRUE,
                                                 num_burn_in = 2000,
                                                 num_iterations_after_burn_in = 1000,
                                                 serialize = TRUE,
                                                 seed = 123)

  saveRDS(bart_model_final, paste0(output_path, "/response_model/bart_model_final.rds"))

}


########
### Validation of model
########

# Dev

data_dev <- hba1c.train %>%
  select(c(patid, pated, posthba1cfinal,
           colnames(bart_model_final$X)))


## Get posteriors
if (class(try(
  
  posteriors_dev <- readRDS(paste0(output_path, "/response_model/assessment/posteriors_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_final, data_dev %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  
  saveRDS(posteriors_dev, paste0(output_path, "/response_model/assessment/posteriors_dev.rds"))
  
}

### residuals calculation
if (class(try(
  
  cred_pred_dev <- readRDS(paste0(output_path, "/response_model/assessment/cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_dev <- calc_resid(data_dev, posteriors_dev, "posthba1cfinal")
  
  saveRDS(cred_pred_dev, paste0(output_path, "/response_model/assessment/cred_pred_dev.rds"))
  
}

# Val

hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test")

data_val <- hba1c.test %>%
  select(c(patid, pated, posthba1cfinal,
           colnames(bart_model_final$X)))


### Get posteriors
if (class(try(
  
  posteriors_val <- readRDS(paste0(output_path, "/response_model/assessment/posteriors_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_val <- bartMachine::bart_machine_get_posterior(bart_model_final, data_val %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  
  saveRDS(posteriors_val, paste0(output_path, "/response_model/assessment/posteriors_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  cred_pred_val <- readRDS(paste0(output_path, "/response_model/assessment/cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_val <- calc_resid(data_val, posteriors_val, "posthba1cfinal")
  
  saveRDS(cred_pred_val, paste0(output_path, "/response_model/assessment/cred_pred_val.rds"))
  
}


plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Residuals of BART response model")



########
### Validation of effect
########

# calculate effects
if (class(try(
  
  effects_summary_dev <- readRDS(paste0(output_path, "/response_model/assessment/effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev <- calc_effect_summary(bart_model_final, data_dev)
  
  saveRDS(effects_summary_dev, paste0(output_path, "/response_model/assessment/effects_summary_dev.rds"))
  
}


# calculate effects
if (class(try(
  
  effects_summary_val <- readRDS(paste0(output_path, "/response_model/assessment/effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_val <- calc_effect_summary(bart_model_final, data_val)
  
  saveRDS(effects_summary_val, paste0(output_path, "/response_model/assessment/effects_summary_val.rds"))
  
}


## plot histogram of effect

plot_effect_1 <- hist_plot(effects_summary_dev, "", -15, 20)

plot_effect_2 <- hist_plot(effects_summary_val, "", -15, 20)



# Validating ATE

prop_scores_dataset <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

prop_scores_individuals_dev <- data_dev %>%
  select(patid, pated) %>%
  left_join(prop_scores_dataset, by = c("patid", "pated")) %>%
  select(prop.score) %>%
  unlist()

prop_scores_individuals_val <- data_val %>%
  select(patid, pated) %>%
  left_join(prop_scores_dataset, by = c("patid", "pated")) %>%
  select(prop.score) %>%
  unlist()

predicted_observed_dev <- data_dev %>%
  cbind(hba1c_diff = effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

predicted_observed_val <- data_val %>%
  cbind(hba1c_diff = effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


if (class(try(
  
  ATE_matching_validation_dev <- readRDS(paste0(output_path, "/response_model/assessment/ATE_matching_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_dev, "posthba1cfinal", prop_scores_individuals_dev)
  
  saveRDS(ATE_matching_validation_dev, paste0(output_path, "/response_model/assessment/ATE_matching_validation_dev.rds"))
  
}


if (class(try(
  
  ATE_matching_validation_val <- readRDS(paste0(output_path, "/response_model/assessment/ATE_matching_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val <- calc_ATE_validation_prop_matching(predicted_observed_val, "posthba1cfinal", prop_scores_individuals_val)
  
  saveRDS(ATE_matching_validation_val, paste0(output_path, "/response_model/assessment/ATE_matching_validation_val.rds"))
  
}

plot_ATE_dev_prop_score <- ATE_plot(ATE_matching_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


plot_ATE_val_prop_score <- ATE_plot(ATE_matching_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


plot_ATE_prop_score_matching <- patchwork::wrap_plots(list(plot_ATE_dev_prop_score, plot_ATE_val_prop_score), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Effects validation propensity score matching", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5))) # center title of full plot

  


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev <- readRDS(paste0(output_path, "/response_model/assessment/ATE_weighting_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "posthba1cfinal", prop_scores_individuals_dev)
  
  saveRDS(ATE_weighting_validation_dev, paste0(output_path, "/response_model/assessment/ATE_weighting_validation_dev.rds"))
  
}

if (class(try(
  
  ATE_weighting_validation_val <- readRDS(paste0(output_path, "/response_model/assessment/ATE_weighting_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1cfinal", prop_scores_individuals_val)
  
  saveRDS(ATE_weighting_validation_val, paste0(output_path, "/response_model/assessment/ATE_weighting_validation_val.rds"))
  
}

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


plot_ATE_prop_score_weighting <- patchwork::wrap_plots(list(plot_ATE_dev_prop_score_weighting, plot_ATE_val_prop_score_weighting), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Effects validation propensity score inverse weighting", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5))) # center title of full plot



pdf(file = "Plots/11.05.BART_response_model.pdf")
plot_residuals

patchwork::wrap_plots(list(plot_effect_1, plot_effect_2), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Treatment effect heterogeneity", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

plot_ATE_prop_score_matching
plot_ATE_prop_score_weighting
dev.off()




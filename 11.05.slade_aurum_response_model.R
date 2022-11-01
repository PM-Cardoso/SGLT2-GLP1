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
  
  saveRDS(vs_bart_model, paste0(output_path, "/response_model/vs_bart_model.rds"))
  
}

# Cross-validation selection
if (class(try(
  
  vs_bart_model_cv <- readRDS(paste0(output_path, "/response_model/vs_bart_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  set.seed(123)
  vs_bart_model_cv <- var_selection_by_permute_cv(bart_model_cv)
  
  ## Long version of the var selection
  
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










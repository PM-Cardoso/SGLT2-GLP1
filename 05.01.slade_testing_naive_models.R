####################
## Description:
##  - In this file we fit a collection of BART models taking a naive approach
##      utilising routine variables / all variables / propensity scores,
##      alternating between their arrangement.
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
dir.create(paste0(output_path, "/Importance"))

## make directory for outputs
dir.create(paste0(output_path, "/Model_fit"))


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))


###############################################################################
###############################################################################
############################### Model Fitting #################################
###############################################################################
###############################################################################
##
## Initially we start with 9866 individuals, dev: 5978, val: 3888
##


#############################
### Complete model of only routine data, no propensity score (n: 8564(4542))
#############################

data_complete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass, 
    ncurrtx,
    drugline, 
    yrdrugstart, 
    t2dmduration, 
    agetx, 
    malesex, 
    Category, 
    hba1cmonth,
    prebmi,
    prealt,
    egfr_ckdepi,
    prehba1cmmol
  ) %>% 
  drop_na() # removed 1302

## Fit Bart model with all variables.

if (class(try(
  
  bart_comp_routine_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_comp_routine_no_prop <- bartMachine::bartMachine(X = data_complete_routine_dev[,-c(1,2,3)],
                                                        y = data_complete_routine_dev[,3],
                                                        use_missing_data = TRUE,
                                                        impute_missingness_with_rf_impute = FALSE,
                                                        impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                        num_trees = 200,
                                                        num_burn_in = 3000,
                                                        num_iterations_after_burn_in = 1000,
                                                        serialize = TRUE)
  
  saveRDS(bart_comp_routine_no_prop, paste0(output_path, "/Model_fit/bart_comp_routine_no_prop.rds"))
  
}



#############################
### Complete model of only routine data, propensity score (n: 8564(4542))
#############################


data_complete_routine_prop_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass, 
    ncurrtx,
    drugline, 
    yrdrugstart, 
    t2dmduration, 
    agetx, 
    malesex, 
    Category, 
    hba1cmonth,
    prebmi,
    prealt,
    egfr_ckdepi,
    prehba1cmmol
  ) %>% 
  drop_na()


## Fit a propensity model with all the variables
if (class(try(
  
  bart_comp_routine_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # GLP1 is considered target so all the probabilities should be 1-prob
  bart_comp_routine_prop <- bartMachine::bartMachine(X = data_complete_routine_prop_dev[,-c(1,2,3,4,7)],
                                                     y = data_complete_routine_prop_dev[,4],
                                                     use_missing_data = TRUE,
                                                     impute_missingness_with_rf_impute = FALSE,
                                                     impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                     num_trees = 200,
                                                     num_burn_in = 3000,
                                                     num_iterations_after_burn_in = 1000,
                                                     serialize = TRUE)
  
  saveRDS(bart_comp_routine_prop, paste0(output_path, "/Model_fit/bart_comp_routine_prop.rds"))

  }


# Add prop score to dataset
data_complete_routine_prop_dev_prop <- data_complete_routine_prop_dev %>% 
  cbind(prop_score = bart_comp_routine_prop$p_hat_train)


## Fit Bart model with all variables.

if (class(try(
  
  bart_comp_routine_prop_model <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_comp_routine_prop_model <- bartMachine::bartMachine(X = data_complete_routine_prop_dev_prop[,-c(1,2,3)],
                                                           y = data_complete_routine_prop_dev_prop[,3],
                                                           use_missing_data = TRUE,
                                                           impute_missingness_with_rf_impute = FALSE,
                                                           impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                           num_trees = 200,
                                                           num_burn_in = 3000,
                                                           num_iterations_after_burn_in = 1000,
                                                           serialize = TRUE)
  
  saveRDS(bart_comp_routine_prop_model, paste0(output_path, "/Model_fit/bart_comp_routine_prop_model.rds"))
  
}



#############################
### Incomplete model of only routine data, no propensity score (n: 9866(5978))
#############################


data_incomplete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass, 
    ncurrtx,
    drugline, 
    yrdrugstart, 
    t2dmduration, 
    agetx, 
    malesex, 
    Category, 
    hba1cmonth,
    prebmi,
    prealt,
    egfr_ckdepi,
    prehba1cmmol
  ) 

## Fit Bart model with all variables.

if (class(try(
  
  bart_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_routine_no_prop <- bartMachine::bartMachine(X = data_incomplete_routine_dev[,-c(1,2,3)],
                                                          y = data_incomplete_routine_dev[,3],
                                                          use_missing_data = TRUE,
                                                          impute_missingness_with_rf_impute = FALSE,
                                                          impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                          num_trees = 200,
                                                          num_burn_in = 3000,
                                                          num_iterations_after_burn_in = 1000,
                                                          serialize = TRUE)
  
  saveRDS(bart_incomp_routine_no_prop, paste0(output_path, "/Model_fit/bart_incomp_routine_no_prop.rds"))
  
}



#############################
### Incomplete model of all data, no propensity score (n: 9866(5978))
#############################


data_incomplete_dev <- final.dev


if (class(try(
  
  bart_incomp_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_no_prop <- bartMachine::bartMachine(X = data_incomplete_dev[,-c(1,2,3,4,9,10)],
                                                  y = data_incomplete_dev[,4],
                                                  use_missing_data = TRUE,
                                                  impute_missingness_with_rf_impute = FALSE,
                                                  impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                  num_trees = 200,
                                                  num_burn_in = 3000,
                                                  num_iterations_after_burn_in = 1000,
                                                  serialize = TRUE)
  
  saveRDS(bart_incomp_no_prop, paste0(output_path, "/Model_fit/bart_incomp_no_prop.rds"))
  
}


#############################
### VARIABLE SELECTION 1: Incomplete model of all data, no propensity score (n: 9866(5978))
#############################

# Variable selection from Incomplete no propensity model

if (class(try(
  
  vs_incomp_no_prop <- readRDS(paste0(output_path, "/Importance/vs_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  vs_incomp_no_prop <- var_selection_by_permute_cv(bart_incomp_no_prop,
                                                   k_folds = 15,
                                                   num_permute_samples = 100,
                                                   num_trees_pred_cv = 100)
  
  ## Long version of the var selection
  # [1] "drugclass_GLP1" "drugline_2"     "drugline_3"     "drugline_4"
  # [5] "drugline_5"     "egfr_ckdepi"    "hba1cmonth"     "ncurrtx_1"
  # [9] "ncurrtx_2"      "ncurrtx_3"      "prehba1cmmol"   "yrdrugstart"
  
  
  saveRDS(vs_incomp_no_prop, paste0(output_path, "/Importance/vs_incomp_no_prop.rds"))
  
}


data_incomplete_dev_var_select <- final.dev %>%
  select(c(patid,
           pateddrug,
           posthba1c_final,
           # agetx,
           drugclass, 
           drugline, 
           egfr_ckdepi, 
           hba1cmonth, 
           # malesex, 
           ncurrtx, 
           # prealt, 
           prehba1cmmol, 
           # score, 
           yrdrugstart
           # predrug.5yrrecent.pad
           )
  )

# Fit Bart model with variables selected

if (class(try(
  
  bart_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_no_prop_var_select.rds"))
  
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_no_prop_var_select <- bartMachine::bartMachine(X = data_incomplete_dev_var_select[,-c(1,2,3)],
                                                             y = data_incomplete_dev_var_select[,3],
                                                             use_missing_data = TRUE,
                                                             impute_missingness_with_rf_impute = FALSE,
                                                             impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                             num_trees = 200,
                                                             num_burn_in = 3000,
                                                             num_iterations_after_burn_in = 1000,
                                                             serialize = TRUE)
  
  saveRDS(bart_incomp_no_prop_var_select, paste0(output_path, "/Model_fit/bart_incomp_no_prop_var_select.rds"))
  
}




#############################
### Incomplete data, propensity score
#############################

data_incomplete_dev <- final.dev


## Fit a propensity model with all the variables
if (class(try(
  
  bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # GLP1 is considered target so all the probabilities should be 1-prob
  bart_incomp_prop <- bartMachine::bartMachine(X = data_incomplete_dev[,-c(1,2,3,4,5,8,9,10)],
                                               y = data_incomplete_dev[,5],
                                               use_missing_data = TRUE,
                                               impute_missingness_with_rf_impute = FALSE,
                                               impute_missingness_with_x_j_bar_for_lm = TRUE,
                                               num_trees = 200,
                                               num_burn_in = 3000,
                                               num_iterations_after_burn_in = 1000,
                                               serialize = TRUE)
  
  saveRDS(bart_incomp_prop, paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))
  
}


## Add propensity score to dataset

data_incomplete_dev_prop <- final.dev %>% 
  cbind(prop_score = bart_incomp_prop$p_hat_train)


# Fit Bart model with variables selected
if (class(try(
  
  bart_incomp_prop_model <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model.rds"))
  
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_prop_model <- bartMachine::bartMachine(X = data_incomplete_dev_prop[,-c(1,2,3,4,9,10)],
                                                     y = data_incomplete_dev_prop[,4],
                                                     use_missing_data = TRUE,
                                                     impute_missingness_with_rf_impute = FALSE,
                                                     impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                     num_trees = 200,
                                                     num_burn_in = 3000,
                                                     num_iterations_after_burn_in = 1000,
                                                     serialize = TRUE)
  
  saveRDS(bart_incomp_prop_model, paste0(output_path, "/Model_fit/bart_incomp_prop_model.rds"))
  
}


#############################
### VARIABLE SELECTION 1: Incomplete data, propensity score
#############################

# Variable selection from Incomplete no propensity model

vs_incomp_no_prop <- readRDS(paste0(output_path, "/Importance/vs_incomp_no_prop.rds"))
  

bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))


## Add propensity score to dataset
data_incomplete_dev_prop_var_select <- final.dev %>%
  cbind(prop_score = bart_incomp_prop$p_hat_train) %>%
  select(c(patid,
           pateddrug,
           posthba1c_final,
           # agetx,
           drugclass, 
           drugline, 
           egfr_ckdepi, 
           hba1cmonth, 
           # malesex, 
           ncurrtx, 
           # prealt, 
           prehba1cmmol, 
           # score, 
           yrdrugstart,
           # predrug.5yrrecent.pad,
           prop_score)
  )

# Fit Bart model with variables selected
if (class(try(
  
  bart_incomp_prop_model_var_select_1 <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select_1.rds"))
  
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_prop_model_var_select_1 <- bartMachine::bartMachine(X = data_incomplete_dev_prop_var_select[,-c(1,2,3)],
                                                                  y = data_incomplete_dev_prop_var_select[,3],
                                                                  use_missing_data = TRUE,
                                                                  impute_missingness_with_rf_impute = FALSE,
                                                                  impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                                  num_trees = 200,
                                                                  num_burn_in = 3000,
                                                                  num_iterations_after_burn_in = 1000,
                                                                  serialize = TRUE)
  
  saveRDS(bart_incomp_prop_model_var_select_1, paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select_1.rds"))
  
}



#############################
### VARIABLE SELECTION 2: Incomplete data, propensity score
#############################

# Variable selection from Incomplete no propensity model

if (class(try(
  
  vs_incomp_prop_model <- readRDS(paste0(output_path, "/Importance/vs_incomp_prop_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  vs_incomp_prop_model <- var_selection_by_permute_cv(bart_incomp_prop_model,
                                                      k_folds = 15,
                                                      num_permute_samples = 100,
                                                      num_trees_pred_cv = 100)
  
  # [1] "Category_Ex-smoker" "drugclass_GLP1"     "drugclass_SGLT2"
  # [4] "drugline_2"         "drugline_3"         "drugline_4"
  # [7] "drugline_5"         "egfr_ckdepi"        "hba1cmonth"
  # [10] "malesex_0"          "ncurrtx_0"          "ncurrtx_1"
  # [13] "ncurrtx_2"          "ncurrtx_3"          "prealt"
  # [16] "prehba1cmmol"       "score"              "yrdrugstart"
  
  
  saveRDS(vs_incomp_prop_model, paste0(output_path, "/Importance/vs_incomp_prop_model.rds"))
  
}


bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))


## Add propensity score to dataset

data_incomplete_dev_prop_var_select <- final.dev %>%
  cbind(prop_score = bart_incomp_prop$p_hat_train) %>%
  select(c(patid,
           pateddrug,
           posthba1c_final,
           Category,
           drugclass,
           drugline,
           egfr_ckdepi,
           hba1cmonth,
           malesex,
           ncurrtx,
           prehba1cmmol,
           prealt,
           score,
           yrdrugstart,
           prop_score)
  )


# Fit Bart model with variables selected
if (class(try(
  
  bart_incomp_prop_model_var_select <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select.rds"))
  
  
  , silent = TRUE)) == "try-error") {
  
  bart_incomp_prop_model_var_select <- bartMachine::bartMachine(X = data_incomplete_dev_prop_var_select[,-c(1,2,3)],
                                                                y = data_incomplete_dev_prop_var_select[,3],
                                                                use_missing_data = TRUE,
                                                                impute_missingness_with_rf_impute = FALSE,
                                                                impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                                num_trees = 200,
                                                                num_burn_in = 3000,
                                                                num_iterations_after_burn_in = 1000,
                                                                serialize = TRUE)
  
  saveRDS(bart_incomp_prop_model_var_select, paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select.rds"))
  
}







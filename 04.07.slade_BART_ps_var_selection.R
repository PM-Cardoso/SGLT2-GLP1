####################
## Description:
##  - In this file we:
##    - Fit a BART PS model to all variables.
##    - Perform variable selection on PS model and refit model.
##    - Refit BART PS model
##    - Fit a BART HbA1c model with all variables + propensity score
##    - Perform variable selection of BART HbA1c model
##    - Refit BART HbA1c model
##    - Validation
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
dir.create(paste0(output_path, "/Final_model/model_7"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/model_7/Assessment"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.all.extra.vars
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_allcohort.Rda"))

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

# name:final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))


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

# remove old CVD score, add new score
dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(-preweight, -height)


########
### Fit a propensity score model to all variables in dataset
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_ps_model <- readRDS(paste0(output_path, "/Final_model/model_7/bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model <- bartMachine::bartMachine(X = dataset.dev %>%
                                              select(-patid,
                                                     -pateddrug,
                                                     -bothdrugs,
                                                     -posthba1c_final,
                                                     -drugclass,
                                                     -sglt2subtype,
                                                     -glp1subtype),
                                            y = dataset.dev[,"drugclass"],
                                            use_missing_data = TRUE,
                                            impute_missingness_with_rf_impute = FALSE,
                                            impute_missingness_with_x_j_bar_for_lm = TRUE,
                                            num_trees = 200,
                                            num_burn_in = 3000,
                                            num_iterations_after_burn_in = 1000,
                                            serialize = TRUE)
  
  saveRDS(bart_ps_model, paste0(output_path, "/Final_model/model_7/bart_ps_model.rds"))
  
}

########
### Variable selection of propensity score model
########

if (class(try(
  
  vs_bart_ps_model <- readRDS(paste0(output_path, "/Final_model/model_7/vs_bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/4.7.prop_model_vs.pdf", width = 18, height = 11)
  # error with cv
  vs_bart_ps_model <- var_selection_by_permute(bart_ps_model)
  dev.off()
  
  ## Long version of the var selection
  # [1] "yrdrugstart"         "prebmi"              "t2dmduration"
  # [4] "drugline_2"          "drugline_5"          "prehba1cmmol"
  # [7] "egfr_ckdepi"         "ncurrtx_3"           "drugline_3"
  # [10] "Category_Non-smoker"
  
  saveRDS(vs_bart_ps_model, paste0(output_path, "/Final_model/model_7/vs_bart_ps_model.rds"))
  
}


########
### Refit PS model with selected vars
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_ps_model_final <- readRDS(paste0(output_path, "/Final_model/model_7/bart_ps_model_final.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model_final <- bartMachine::bartMachine(X = dataset.dev %>%
                                                        select(
                                                          yrdrugstart,
                                                          prebmi,
                                                          t2dmduration,
                                                          drugline,
                                                          prehba1cmmol,
                                                          egfr_ckdepi,
                                                          ncurrtx,
                                                          Category
                                                        ),
                                                      y = dataset.dev[,"drugclass"],
                                                      use_missing_data = TRUE,
                                                      impute_missingness_with_rf_impute = FALSE,
                                                      impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                      num_trees = 200,
                                                      num_burn_in = 3000,
                                                      num_iterations_after_burn_in = 1000,
                                                      serialize = TRUE)
  
  saveRDS(bart_ps_model_final, paste0(output_path, "/Final_model/model_7/bart_ps_model_final.rds"))
  
}


########
### Fit BART model to all variables + propensity score
########


# Add propensity score as a possible variable
dataset.dev <- dataset.dev %>%
  cbind(prop_score = 1 - prop_model_final$p_hat_train)


## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_model <- readRDS(paste0(output_path, "/Final_model/model_7/bart_model.rds"))
  
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
  
  saveRDS(bart_model, paste0(output_path, "/Final_model/model_7/bart_model.rds"))
  
}

########
### Variable selection of BART model
########

# One-off variable selection
if (class(try(
  
  vs_bart_model <- readRDS(paste0(output_path, "/Final_model/model_7/vs_bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/4.7.response_model_vs.pdf", width = 18, height = 11)
  vs_bart_model <- var_selection_by_permute(bart_model)
  dev.off()
  
  ## Long version of the var selection
  
  saveRDS(vs_bart_model, paste0(output_path, "/Final_model/model_7/vs_bart_model.rds"))
  
}

# Cross-validation selection
if (class(try(
  
  vs_bart_model_cv <- readRDS(paste0(output_path, "/Final_model/model_7/vs_bart_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  vs_bart_model_cv <- var_selection_by_permute_cv(bart_model)
  
  ## Long version of the var selection
  
  saveRDS(vs_bart_model_cv, paste0(output_path, "/Final_model/model_7/vs_bart_model_cv.rds"))
  
}


# ########
# ### Refit BART model with selected variables
# ########
# 
# ## Fit final model using selected variables to estimate HbA1c outcome
# if (class(try(
#   
#   bart_model_final <- readRDS(paste0(output_path, "/Final_model/model_7/bart_model_final.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bart_model_final <- bartMachine::bartMachine(X = dataset.dev %>%
#                                                  select(
#                                                    prehba1cmmol,
#                                                    prebil,
#                                                    drugclass,
#                                                    yrdrugstart,
#                                                    hba1cmonth,
#                                                    score.excl.mi,
#                                                    drugline,
#                                                    malesex,
#                                                    ncurrtx
#                                                  ),
#                                                y = dataset.dev[,"posthba1c_final"],
#                                                use_missing_data = TRUE,
#                                                impute_missingness_with_rf_impute = FALSE,
#                                                impute_missingness_with_x_j_bar_for_lm = TRUE,
#                                                num_trees = 200,
#                                                num_burn_in = 3000,
#                                                num_iterations_after_burn_in = 1000,
#                                                serialize = TRUE)
#   
#   saveRDS(bart_model_final, paste0(output_path, "/Final_model/model_7/bart_model_final.rds"))
#   
# }





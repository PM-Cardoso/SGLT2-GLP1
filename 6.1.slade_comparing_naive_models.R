####################
## Description:
##  - In this file we make a collected of plots comparing assessment
##      measurements of several models from 5.1.
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")


library(tidyverse)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)


## make directory for outputs
dir.create(paste0(output_path, "/Assessment"))



###############################################################################
###############################################################################
######################### Read Data / Model In ################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))
# name: final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))


#############################
### Complete model of only routine data, no propensity score (n: 8564(4542))
#############################

bart_comp_routine_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_no_prop.rds"))


# Dev
data_complete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_no_prop$X)
  ) %>% 
  drop_na()


## Get posteriors
if (class(try(
  
  posteriors_comp_routine_no_prop_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_comp_routine_no_prop_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_comp_routine_no_prop_dev <- bartMachine::bart_machine_get_posterior(bart_comp_routine_no_prop, data_complete_routine_dev %>%
                                                              select(
                                                                colnames(bart_comp_routine_no_prop$X)
                                                              ))
  
  saveRDS(posteriors_comp_routine_no_prop_dev, paste0(output_path, "/Assessment/posteriors_comp_routine_no_prop_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  comp_routine_no_prop_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/comp_routine_no_prop_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_no_prop_cred_pred_dev <- calc_resid(data_complete_routine_dev, posteriors_comp_routine_no_prop_dev)
  
  saveRDS(comp_routine_no_prop_cred_pred_dev, paste0(output_path, "/Assessment/comp_routine_no_prop_cred_pred_dev.rds"))
  
}



# Val
data_complete_routine_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_no_prop$X)
  ) %>% 
  drop_na()


## Get posteriors
if (class(try(
  
  posteriors_comp_routine_no_prop_val <- readRDS(paste0(output_path, "/Assessment/posteriors_comp_routine_no_prop_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_comp_routine_no_prop_val <- bartMachine::bart_machine_get_posterior(bart_comp_routine_no_prop, data_complete_routine_val %>%
                                                                                   select(
                                                                                     colnames(bart_comp_routine_no_prop$X)
                                                                                   ))
  
  saveRDS(posteriors_comp_routine_no_prop_val, paste0(output_path, "/Assessment/posteriors_comp_routine_no_prop_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  comp_routine_no_prop_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/comp_routine_no_prop_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_no_prop_cred_pred_val <- calc_resid(data_complete_routine_val, posteriors_comp_routine_no_prop_val)
  
  saveRDS(comp_routine_no_prop_cred_pred_val, paste0(output_path, "/Assessment/comp_routine_no_prop_cred_pred_val.rds"))
  
}


# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/assessment_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_complete_routine_dev, posteriors_comp_routine_no_prop_dev)
  
  assessment_values_val <- calc_assessment(data_complete_routine_val, posteriors_comp_routine_no_prop_val)
  
  assessment_comp_routine_no_prop <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_comp_routine_no_prop, paste0(output_path, "/Assessment/assessment_comp_routine_no_prop.rds"))
  
}


plot_comp_routine_no_prop <- resid_plot(comp_routine_no_prop_cred_pred_dev,
                                        comp_routine_no_prop_cred_pred_val, 
                                        "Model fitting: Complete data (Routine/No propensity score)")









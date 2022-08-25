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



#############################
### Complete model of only routine data, propensity score (n: 8564(4542))
#############################

bart_comp_routine_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop.rds"))

bart_comp_routine_prop_model <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop_model.rds"))


# Dev
data_complete_routine_prop_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_no_prop$X)
  ) %>% 
  drop_na() %>%
  cbind(prop_score = bart_comp_routine_prop$p_hat_train)


## Get posteriors
if (class(try(
  
  posteriors_complete_routine_prop_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_complete_routine_prop_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_complete_routine_prop_dev <- bartMachine::bart_machine_get_posterior(bart_comp_routine_prop_model, data_complete_routine_prop_dev %>%
                                                                                   select(
                                                                                     colnames(bart_comp_routine_prop_model$X)
                                                                                   ))
  
  saveRDS(posteriors_complete_routine_prop_dev, paste0(output_path, "/Assessment/posteriors_complete_routine_prop_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  comp_routine_prop_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/comp_routine_prop_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_prop_cred_pred_dev <- calc_resid(data_complete_routine_prop_dev, posteriors_complete_routine_prop_dev)
  
  saveRDS(comp_routine_prop_cred_pred_dev, paste0(output_path, "/Assessment/comp_routine_prop_cred_pred_dev.rds"))
  
}


# Val
data_complete_routine_prop_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_no_prop$X)
  ) %>% 
  drop_na()

## calculate propensity score
if (class(try(
  
  prop_score_complete_routine_prop <- readRDS(paste0(output_path,"/Assessment/prop_score_complete_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop_score_complete_routine_prop <- predict(bart_comp_routine_prop, data_complete_routine_prop_val %>%
                                                select(
                                                  colnames(bart_comp_routine_prop$X)
                                                ))
  
  saveRDS(prop_score_complete_routine_prop, paste0(output_path, "/Assessment/prop_score_complete_routine_prop.rds"))
  
}

# Add 
data_complete_routine_prop_val <- data_complete_routine_prop_val %>% 
  cbind(prop_score = prop_score_complete_routine_prop)


## Get posteriors
if (class(try(
  
  posteriors_complete_routine_prop_val <- readRDS(paste0(output_path, "/Assessment/posteriors_complete_routine_prop_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_complete_routine_prop_val <- bartMachine::bart_machine_get_posterior(bart_comp_routine_prop_model, data_complete_routine_prop_val %>%
                                                                                    select(
                                                                                      colnames(bart_comp_routine_prop_model$X)
                                                                                    ))
  
  saveRDS(posteriors_complete_routine_prop_val, paste0(output_path, "/Assessment/posteriors_complete_routine_prop_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  comp_routine_prop_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/comp_routine_prop_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_prop_cred_pred_val <- calc_resid(data_complete_routine_prop_val, posteriors_complete_routine_prop_val)
  
  saveRDS(comp_routine_prop_cred_pred_val, paste0(output_path, "/Assessment/comp_routine_prop_cred_pred_val.rds"))
  
}


# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/assessment_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_complete_routine_prop_dev, posteriors_complete_routine_prop_dev)
  
  assessment_values_val <- calc_assessment(data_complete_routine_prop_val, posteriors_complete_routine_prop_val)
  
  assessment_comp_routine_prop <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_comp_routine_prop, paste0(output_path, "/Assessment/assessment_comp_routine_prop.rds"))
  
}


plot_comp_routine_prop <- resid_plot(comp_routine_prop_cred_pred_dev,
                                     comp_routine_prop_cred_pred_val, 
                                     "Model fitting: Complete data (Routine/Propensity score)")


#############################
### Incomplete model of only routine data, no propensity score (n: 9866(5978))
#############################


bart_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_routine_no_prop.rds"))

# Dev
data_incomplete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_routine_no_prop$X)
  )


## Get posteriors
if (class(try(
  
  posteriors_incomp_routine_no_prop_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_routine_no_prop_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_routine_no_prop_dev <- bartMachine::bart_machine_get_posterior(bart_incomp_routine_no_prop, data_incomplete_routine_dev %>%
                                                                                   select(
                                                                                     colnames(bart_incomp_routine_no_prop$X)
                                                                                   ))
  
  saveRDS(posteriors_incomp_routine_no_prop_dev, paste0(output_path, "/Assessment/posteriors_incomp_routine_no_prop_dev.rds"))
  
  
}


### residuals calculation
if (class(try(
  
  incomp_routine_no_prop_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/incomp_routine_no_prop_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_routine_no_prop_cred_pred_dev <- calc_resid(data_incomplete_routine_dev, posteriors_incomp_routine_no_prop_dev)
  
  saveRDS(incomp_routine_no_prop_cred_pred_dev, paste0(output_path, "/Assessment/incomp_routine_no_prop_cred_pred_dev.rds"))
  
}

# Val
data_incomplete_routine_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_routine_no_prop$X)
  )

## Get posteriors
if (class(try(
  
  posteriors_incomp_routine_no_prop_val <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_routine_no_prop_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_routine_no_prop_val <- bartMachine::bart_machine_get_posterior(bart_incomp_routine_no_prop, data_incomplete_routine_val %>%
                                                                                     select(
                                                                                       colnames(bart_incomp_routine_no_prop$X)
                                                                                     ))
  
  saveRDS(posteriors_incomp_routine_no_prop_val, paste0(output_path, "/Assessment/posteriors_incomp_routine_no_prop_val.rds"))
  
  
}


### residuals calculation
if (class(try(
  
  incomp_routine_no_prop_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/incomp_routine_no_prop_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_routine_no_prop_cred_pred_val <- calc_resid(data_incomplete_routine_val, posteriors_incomp_routine_no_prop_val)
  
  saveRDS(incomp_routine_no_prop_cred_pred_val, paste0(output_path, "/Assessment/incomp_routine_no_prop_cred_pred_val.rds"))
  
}



# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_routine_dev, posteriors_incomp_routine_no_prop_dev)
  
  assessment_values_val <- calc_assessment(data_incomplete_routine_val, posteriors_incomp_routine_no_prop_val)
  
  assessment_incomp_routine_no_prop <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_routine_no_prop, paste0(output_path, "/Assessment/assessment_incomp_routine_no_prop.rds"))
  
}


plot_incomp_routine_no_prop <- resid_plot(incomp_routine_no_prop_cred_pred_dev,
                                          incomp_routine_no_prop_cred_pred_val, 
                                          "Model fitting: Incomplete data (Routine/No propensity score)")



#############################
### Incomplete model of all data, no propensity score (n: 9866(5978))
#############################

bart_incomp_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_no_prop.rds"))

# Dev
data_incomplete_dev <- final.dev

## Get posteriors
if (class(try(
  
  posteriors_incomp_no_prop_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_no_prop_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_no_prop_dev <- bartMachine::bart_machine_get_posterior(bart_incomp_no_prop, data_incomplete_dev %>%
                                                                                     select(
                                                                                       colnames(bart_incomp_no_prop$X)
                                                                                     ))
  
  saveRDS(posteriors_incomp_no_prop_dev, paste0(output_path, "/Assessment/posteriors_incomp_no_prop_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_no_prop_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/incomp_no_prop_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_cred_pred_dev <- calc_resid(data_incomplete_dev, posteriors_incomp_no_prop_dev)
  
  saveRDS(incomp_no_prop_cred_pred_dev, paste0(output_path, "/Assessment/incomp_no_prop_cred_pred_dev.rds"))
  
}


# Val
data_incomplete_val <- final.val

## Get posteriors
if (class(try(
  
  posteriors_incomp_no_prop_val <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_no_prop_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_no_prop_val <- bartMachine::bart_machine_get_posterior(bart_incomp_no_prop, data_incomplete_val %>%
                                                                             select(
                                                                               colnames(bart_incomp_no_prop$X)
                                                                             ))
  
  saveRDS(posteriors_incomp_no_prop_val, paste0(output_path, "/Assessment/posteriors_incomp_no_prop_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_no_prop_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/incomp_no_prop_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_cred_pred_val <- calc_resid(data_incomplete_val, posteriors_incomp_no_prop_val)
  
  saveRDS(incomp_no_prop_cred_pred_val, paste0(output_path, "/Assessment/incomp_no_prop_cred_pred_val.rds"))
  
}


# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_dev, posteriors_incomp_no_prop_dev)
  
  assessment_values_val <- calc_assessment(data_incomplete_val, posteriors_incomp_no_prop_val)
  
  assessment_incomp_no_prop <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_no_prop, paste0(output_path, "/Assessment/assessment_incomp_no_prop.rds"))
  
}


plot_incomp_no_prop <- resid_plot(incomp_no_prop_cred_pred_dev,
                                  incomp_no_prop_cred_pred_val, 
                                  "Model fitting: Incomplete data (No propensity score)")



#############################
### VARIABLE SELECTION: Incomplete model of all data, no propensity score (n: 9866(5978))
#############################

bart_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_no_prop_var_select.rds"))

# Dev
data_incomplete_dev_var_select <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_no_prop_var_select$X)
  )


## Get posteriors
if (class(try(
  
  posteriors_incomp_no_prop_var_select_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_no_prop_var_select_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_no_prop_var_select_dev <- bartMachine::bart_machine_get_posterior(bart_incomp_no_prop_var_select, data_incomplete_dev_var_select %>%
                                                                             select(
                                                                               colnames(bart_incomp_no_prop_var_select$X)
                                                                             ))
  
  saveRDS(posteriors_incomp_no_prop_var_select_dev, paste0(output_path, "/Assessment/posteriors_incomp_no_prop_var_select_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_no_prop_var_select_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/incomp_no_prop_var_select_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_var_select_cred_pred_dev <- calc_resid(data_incomplete_dev_var_select, posteriors_incomp_no_prop_var_select_dev)
  
  saveRDS(incomp_no_prop_var_select_cred_pred_dev, paste0(output_path, "/Assessment/incomp_no_prop_var_select_cred_pred_dev.rds"))
  
}


# Val
data_incomplete_val_var_select <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_no_prop_var_select$X)
  )

## Get posteriors
if (class(try(
  
  posteriors_incomp_no_prop_var_select_val <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_no_prop_var_select_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_no_prop_var_select_val <- bartMachine::bart_machine_get_posterior(bart_incomp_no_prop_var_select, data_incomplete_val_var_select %>%
                                                                                        select(
                                                                                          colnames(bart_incomp_no_prop_var_select$X)
                                                                                        ))
  
  saveRDS(posteriors_incomp_no_prop_var_select_val, paste0(output_path, "/Assessment/posteriors_incomp_no_prop_var_select_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_no_prop_var_select_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/incomp_no_prop_var_select_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_var_select_cred_pred_val <- calc_resid(data_incomplete_val_var_select, posteriors_incomp_no_prop_var_select_val)
  
  saveRDS(incomp_no_prop_var_select_cred_pred_val, paste0(output_path, "/Assessment/incomp_no_prop_var_select_cred_pred_val.rds"))
  
}

# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_no_prop_val_select <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_no_prop_val_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_dev_var_select, posteriors_incomp_no_prop_var_select_dev)
  
  assessment_values_val <- calc_assessment(data_incomplete_val_var_select, posteriors_incomp_no_prop_var_select_val)
  
  assessment_incomp_no_prop_val_select <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_no_prop_val_select, paste0(output_path, "/Assessment/assessment_incomp_no_prop_val_select.rds"))
  
}


plot_incomp_no_prop_var_select <- resid_plot(incomp_no_prop_var_select_cred_pred_dev,
                                             incomp_no_prop_var_select_cred_pred_val, 
                                             "Model fitting: Variable Selection 1, Incomplete data (No propensity score)")



#############################
### Incomplete model of all data, with propensity score
#############################

bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))

bart_incomp_prop_model <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model.rds"))


# Dev
data_incomplete_dev <- final.dev %>%
  cbind(prop_score = bart_incomp_prop$p_hat_train)

## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_dev <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_dev <- bartMachine::bart_machine_get_posterior(bart_comp_routine_prop_model, data_incomplete_dev %>%
                                                                                    select(
                                                                                      colnames(bart_comp_routine_prop_model$X)
                                                                                    ))
  
  saveRDS(posteriors_incomp_prop_dev, paste0(output_path, "/Assessment/posteriors_incomp_prop_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_dev <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_dev <- calc_resid(data_incomplete_dev, posteriors_incomp_prop_dev)
  
  saveRDS(incomp_prop_cred_pred_dev, paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev.rds"))
  
}

# Val
data_incomplete_val <- final.val

## calculate propensity score
if (class(try(
  
  prop_score_incomplete_prop <- readRDS(paste0(output_path,"/Assessment/prop_score_incomplete_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop_score_incomplete_prop <- predict(bart_incomp_prop, data_incomplete_val %>%
                                                select(
                                                  colnames(bart_incomp_prop$X)
                                                ))
  
  saveRDS(prop_score_incomplete_prop, paste0(output_path, "/Assessment/prop_score_incomplete_prop.rds"))
  
}

# Add 
data_incomplete_val <- data_incomplete_val %>%
  cbind(prop_score = prop_score_incomplete_prop)

## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_val <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_val <- bartMachine::bart_machine_get_posterior(bart_comp_routine_prop_model, data_incomplete_val %>%
                                                                          select(
                                                                            colnames(bart_comp_routine_prop_model$X)
                                                                          ))
  
  saveRDS(posteriors_incomp_prop_val, paste0(output_path, "/Assessment/posteriors_incomp_prop_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_val <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_val <- calc_resid(data_incomplete_val, posteriors_incomp_prop_val)
  
  saveRDS(incomp_prop_cred_pred_val, paste0(output_path, "/Assessment/incomp_prop_cred_pred_val.rds"))
  
}

# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_prop <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_dev, posteriors_incomp_prop_dev)
  
  assessment_values_val <- calc_assessment(data_incomplete_val, posteriors_incomp_prop_val)
  
  assessment_incomp_prop <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_prop, paste0(output_path, "/Assessment/assessment_incomp_prop.rds"))
  
}


plot_incomp_prop <- resid_plot(comp_routine_prop_cred_pred_dev,
                               comp_routine_prop_cred_pred_val, 
                               "Model fitting: Incomplete data (Propensity score)")






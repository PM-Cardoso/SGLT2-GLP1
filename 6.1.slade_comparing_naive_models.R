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

# calculate effects
if (class(try(
  
  comp_routine_no_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_no_prop_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_no_prop_effects_summary_dev <- calc_effect_summary(bart_comp_routine_no_prop, data_complete_routine_dev)
  
  saveRDS(comp_routine_no_prop_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/comp_routine_no_prop_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_complete_routine_dev %>%
  cbind(hba1c_diff = comp_routine_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_no_prop_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")



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

# calculate effects
if (class(try(
  
  comp_routine_no_prop_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_no_prop_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_no_prop_effects_summary_val <- calc_effect_summary(bart_comp_routine_no_prop, data_complete_routine_val)
  
  saveRDS(comp_routine_no_prop_effects_summary_val, paste0(output_path, "/Final_model/Assessment/comp_routine_no_prop_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_complete_routine_val %>%
  cbind(hba1c_diff = comp_routine_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_no_prop_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")


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
    colnames(bart_comp_routine_prop_model$X)[which(colnames(bart_comp_routine_prop_model$X) != "prop_score")]
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

# calculate effects
if (class(try(
  
  comp_routine_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_prop_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_prop_effects_summary_dev <- calc_effect_summary(bart_comp_routine_prop_model, data_complete_routine_prop_dev)
  
  saveRDS(comp_routine_prop_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/comp_routine_prop_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_complete_routine_prop_dev %>%
  cbind(hba1c_diff = comp_routine_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_prop_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


# Val
data_complete_routine_prop_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_prop_model$X)[which(colnames(bart_comp_routine_prop_model$X) != "prop_score")]
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

# calculate effects
if (class(try(
  
  comp_routine_prop_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_prop_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_prop_effects_summary_val <- calc_effect_summary(bart_comp_routine_prop_model, data_complete_routine_prop_val)
  
  saveRDS(comp_routine_prop_effects_summary_val, paste0(output_path, "/Final_model/Assessment/comp_routine_prop_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_complete_routine_prop_val %>%
  cbind(hba1c_diff = comp_routine_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_prop_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")


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

# calculate effects
if (class(try(
  
  incomp_routine_no_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_routine_no_prop_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_routine_no_prop_effects_summary_dev <- calc_effect_summary(bart_incomp_routine_no_prop, data_incomplete_routine_dev)
  
  saveRDS(incomp_routine_no_prop_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_routine_no_prop_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_routine_dev %>%
  cbind(hba1c_diff = incomp_routine_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_routine_no_prop_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


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

# calculate effects
if (class(try(
  
  incomp_routine_no_prop_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_routine_no_prop_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_routine_no_prop_effects_summary_val <- calc_effect_summary(bart_incomp_routine_no_prop, data_incomplete_routine_val)
  
  saveRDS(incomp_routine_no_prop_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_routine_no_prop_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_routine_val %>%
  cbind(hba1c_diff = incomp_routine_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_routine_no_prop_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")


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

# calculate effects
if (class(try(
  
  incomp_no_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_no_prop_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_effects_summary_dev <- calc_effect_summary(bart_incomp_no_prop, data_incomplete_dev)
  
  saveRDS(incomp_no_prop_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_no_prop_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_dev %>%
  cbind(hba1c_diff = incomp_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


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

# calculate effects
if (class(try(
  
  incomp_no_prop_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_no_prop_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_effects_summary_val <- calc_effect_summary(bart_incomp_no_prop, data_incomplete_val)
  
  saveRDS(incomp_no_prop_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_no_prop_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_val %>%
  cbind(hba1c_diff = incomp_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")



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

# calculate effects
if (class(try(
  
  incomp_no_prop_var_select_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_no_prop_var_select_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_var_select_effects_summary_dev <- calc_effect_summary(bart_incomp_no_prop_var_select, data_incomplete_dev_var_select)
  
  saveRDS(incomp_no_prop_var_select_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_no_prop_var_select_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_dev_var_select %>%
  cbind(hba1c_diff = incomp_no_prop_var_select_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_var_select_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


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

# calculate effects
if (class(try(
  
  incomp_no_prop_var_select_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_no_prop_var_select_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_no_prop_var_select_effects_summary_val <- calc_effect_summary(bart_incomp_no_prop_var_select, data_incomplete_val_var_select)
  
  saveRDS(incomp_no_prop_var_select_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_no_prop_var_select_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_val_var_select %>%
  cbind(hba1c_diff = incomp_no_prop_var_select_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_var_select_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")



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
  
  posteriors_incomp_prop_dev <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model, data_incomplete_dev %>%
                                                                                    select(
                                                                                      colnames(bart_incomp_prop_model$X)
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

# calculate effects
if (class(try(
  
  incomp_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_effects_summary_dev <- calc_effect_summary(bart_incomp_prop_model, data_incomplete_dev)
  
  saveRDS(incomp_prop_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_prop_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_dev %>%
  cbind(hba1c_diff = incomp_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


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
  
  posteriors_incomp_prop_val <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model, data_incomplete_val %>%
                                                                          select(
                                                                            colnames(bart_incomp_prop_model$X)
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

# calculate effects
if (class(try(
  
  incomp_prop_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_effects_summary_val <- calc_effect_summary(bart_incomp_prop_model, data_incomplete_val)
  
  saveRDS(incomp_prop_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_prop_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_val %>%
  cbind(hba1c_diff = incomp_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")

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


plot_incomp_prop <- resid_plot(incomp_prop_cred_pred_dev,
                               incomp_prop_cred_pred_val, 
                               "Model fitting: Incomplete data (Propensity score)")



#############################
### VARIABLE SELECTION 1: Incomplete model of all data, with propensity score
#############################

bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))

bart_incomp_prop_model_var_select_1 <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select_1.rds"))


# Dev
data_incomplete_dev_var_select_1 <- final.dev %>%
  cbind(prop_score = bart_incomp_prop$p_hat_train) %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_prop_model_var_select_1$X)
  )


## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_dev_var_select_1 <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_dev_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_dev_var_select_1 <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model_var_select_1, data_incomplete_dev_var_select_1 %>%
                                                                          select(
                                                                            colnames(bart_incomp_prop_model_var_select_1$X)
                                                                          ))
  
  saveRDS(posteriors_incomp_prop_dev_var_select_1, paste0(output_path, "/Assessment/posteriors_incomp_prop_dev_var_select_1.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_dev_var_select_1 <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_dev_var_select_1 <- calc_resid(data_incomplete_dev_var_select_1, posteriors_incomp_prop_dev_var_select_1)
  
  saveRDS(incomp_prop_cred_pred_dev_var_select_1, paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev_var_select_1.rds"))
  
}

# calculate effects
if (class(try(
  
  incomp_prop_var_select_1_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_1_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_var_select_1_effects_summary_dev <- calc_effect_summary(bart_incomp_prop_model_var_select_1, data_incomplete_dev_var_select_1)
  
  saveRDS(incomp_prop_var_select_1_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_1_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_dev_var_select_1 %>%
  cbind(hba1c_diff = incomp_prop_var_select_1_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_1_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")


# Val
data_incomplete_val_var_select_1 <- final.val

## calculate propensity score
if (class(try(
  
  prop_score_incomplete_prop <- readRDS(paste0(output_path,"/Assessment/prop_score_incomplete_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop_score_incomplete_prop <- predict(bart_incomp_prop, data_incomplete_val_var_select_1 %>%
                                          select(
                                            colnames(bart_incomp_prop$X)
                                          ))
  
  saveRDS(prop_score_incomplete_prop, paste0(output_path, "/Assessment/prop_score_incomplete_prop.rds"))
  
}

# Add 
data_incomplete_val_var_select_1 <- data_incomplete_val_var_select_1 %>%
  cbind(prop_score = prop_score_incomplete_prop) %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_prop_model_var_select_1$X)
  )

## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_val_var_select_1 <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_val_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_val_var_select_1 <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model_var_select_1, data_incomplete_val_var_select_1 %>%
                                                                                       select(
                                                                                         colnames(bart_incomp_prop_model_var_select_1$X)
                                                                                       ))
  
  saveRDS(posteriors_incomp_prop_val_var_select_1, paste0(output_path, "/Assessment/posteriors_incomp_prop_val_var_select_1.rds"))
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_val_var_select_1 <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_val_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_val_var_select_1 <- calc_resid(data_incomplete_val_var_select_1, posteriors_incomp_prop_val_var_select_1)
  
  saveRDS(incomp_prop_cred_pred_val_var_select_1, paste0(output_path, "/Assessment/incomp_prop_cred_pred_val_var_select_1.rds"))
  
}

# calculate effects
if (class(try(
  
  incomp_prop_var_select_1_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_1_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_var_select_1_effects_summary_val <- calc_effect_summary(bart_incomp_prop_model_var_select_1, data_incomplete_val_var_select_1)
  
  saveRDS(incomp_prop_var_select_1_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_1_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_val_var_select_1 %>%
  cbind(hba1c_diff = incomp_prop_var_select_1_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_1_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")



# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_dev_var_select_1, posteriors_incomp_prop_dev_var_select_1)
  
  assessment_values_val <- calc_assessment(data_incomplete_val_var_select_1, posteriors_incomp_prop_val_var_select_1)
  
  assessment_incomp_prop_var_select_1 <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_prop_var_select_1, paste0(output_path, "/Assessment/assessment_incomp_prop_var_select_1.rds"))
  
}


plot_incomp_prop_var_select_1 <- resid_plot(incomp_prop_cred_pred_dev_var_select_1,
                                            incomp_prop_cred_pred_val_var_select_1, 
                                            "Model fitting: Variable Selection 1, Incomplete data (Propensity score)")



#############################
### VARIABLE SELECTION 2: Incomplete model of all data, with propensity score
#############################

bart_incomp_prop <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop.rds"))

bart_incomp_prop_model_var_select <- readRDS(paste0(output_path, "/Model_fit/bart_incomp_prop_model_var_select.rds"))

# Dev
data_incomplete_dev_var_select <- final.dev %>%
  cbind(prop_score = bart_incomp_prop$p_hat_train) %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_prop_model_var_select$X)
  )

## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_dev_var_select <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_dev_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_dev_var_select <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model_var_select, data_incomplete_dev_var_select %>%
                                                                                       select(
                                                                                         colnames(bart_incomp_prop_model_var_select$X)
                                                                                       ))
  
  saveRDS(posteriors_incomp_prop_dev_var_select, paste0(output_path, "/Assessment/posteriors_incomp_prop_dev_var_select.rds"))
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_dev_var_select <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_dev_var_select <- calc_resid(data_incomplete_dev_var_select, posteriors_incomp_prop_dev_var_select)
  
  saveRDS(incomp_prop_cred_pred_dev_var_select, paste0(output_path, "/Assessment/incomp_prop_cred_pred_dev_var_select.rds"))
  
}

# calculate effects
if (class(try(
  
  incomp_prop_var_select_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_var_select_effects_summary_dev <- calc_effect_summary(bart_incomp_prop_model_var_select, data_incomplete_dev_var_select)
  
  saveRDS(incomp_prop_var_select_effects_summary_dev, paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_effects_summary_dev.rds"))
  
}

## plot effects validation
predicted_observed_dev <- data_incomplete_dev_var_select %>%
  cbind(hba1c_diff = incomp_prop_var_select_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_effects_validation_dev <- plot_full_effects_validation(predicted_observed_dev, dataset = "Dev")



# Val
data_incomplete_val_var_select <- final.val

## calculate propensity score
if (class(try(
  
  prop_score_incomplete_prop <- readRDS(paste0(output_path,"/Assessment/prop_score_incomplete_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  prop_score_incomplete_prop <- predict(bart_incomp_prop, data_incomplete_val_var_select %>%
                                          select(
                                            colnames(bart_incomp_prop$X)
                                          ))
  
  saveRDS(prop_score_incomplete_prop, paste0(output_path, "/Assessment/prop_score_incomplete_prop.rds"))
  
}

# Add 
data_incomplete_val_var_select <- data_incomplete_val_var_select %>%
  cbind(prop_score = prop_score_incomplete_prop) %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_incomp_prop_model_var_select$X)
  )


## Get posteriors
if (class(try(
  
  posteriors_incomp_prop_val_var_select <- readRDS(paste0(output_path, "/Assessment/posteriors_incomp_prop_val_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_incomp_prop_val_var_select <- bartMachine::bart_machine_get_posterior(bart_incomp_prop_model_var_select, data_incomplete_val_var_select %>%
                                                                                     select(
                                                                                       colnames(bart_incomp_prop_model_var_select$X)
                                                                                     ))
  
  saveRDS(posteriors_incomp_prop_val_var_select, paste0(output_path, "/Assessment/posteriors_incomp_prop_val_var_select.rds"))
  
}

### residuals calculation
if (class(try(
  
  incomp_prop_cred_pred_val_var_select <- readRDS(paste0(output_path, "/Assessment/incomp_prop_cred_pred_val_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_cred_pred_val_var_select <- calc_resid(data_incomplete_val_var_select, posteriors_incomp_prop_val_var_select)
  
  saveRDS(incomp_prop_cred_pred_val_var_select, paste0(output_path, "/Assessment/incomp_prop_cred_pred_val_var_select.rds"))
  
}


# calculate effects
if (class(try(
  
  incomp_prop_var_select_effects_summary_val <- readRDS(paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  incomp_prop_var_select_effects_summary_val <- calc_effect_summary(bart_incomp_prop_model_var_select, data_incomplete_val_var_select)
  
  saveRDS(incomp_prop_var_select_effects_summary_val, paste0(output_path, "/Final_model/Assessment/incomp_prop_var_select_effects_summary_val.rds"))
  
}

## plot effects validation
predicted_observed_val <- data_incomplete_val_var_select %>%
  cbind(hba1c_diff = incomp_prop_var_select_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_effects_validation_val <- plot_full_effects_validation(predicted_observed_val, dataset = "Val")


# assessment of R2, RSS, RMSE
if (class(try(
  
  assessment_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/assessment_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  assessment_values_dev <- calc_assessment(data_incomplete_dev_var_select, posteriors_incomp_prop_dev_var_select)
  
  assessment_values_val <- calc_assessment(data_incomplete_val_var_select, posteriors_incomp_prop_val_var_select)
  
  assessment_incomp_prop_var_select <- rbind(
    cbind(t(assessment_values_dev[["r2"]]), Dataset = "Development", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_val[["r2"]]), Dataset = "Validation", statistic = "R2 (bigger is better)"),
    cbind(t(assessment_values_dev[["RSS"]]), Dataset = "Development", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_val[["RSS"]]), Dataset = "Validation", statistic = "RSS (smaller is better)"),
    cbind(t(assessment_values_dev[["RMSE"]]), Dataset = "Development", statistic = "RMSE (smaller is better)"),
    cbind(t(assessment_values_val[["RMSE"]]), Dataset = "Validation", statistic = "RMSE (smaller is better)")
  )
  
  saveRDS(assessment_incomp_prop_var_select, paste0(output_path, "/Assessment/assessment_incomp_prop_var_select.rds"))
  
}


plot_incomp_prop_var_select <- resid_plot(incomp_prop_cred_pred_dev_var_select,
                                            incomp_prop_cred_pred_val_var_select, 
                                            "Model fitting: Variable Selection 2, Incomplete data (Propensity score)")



###############################################################################
###############################################################################
############### Combining Assessment measures for all models ##################
###############################################################################
###############################################################################


assessment <- rbind(
  cbind(assessment_comp_routine_no_prop, Model = "Comp/Routine"),
  cbind(assessment_comp_routine_prop, Model = "Comp/Routine/Prop"),
  cbind(assessment_incomp_routine_no_prop, Model = "Incomp/Routine"),
  cbind(assessment_incomp_no_prop, Model = "Incomp"),
  cbind(assessment_incomp_no_prop_val_select, Model = "Incomp/Var. Selection (1)"),
  cbind(assessment_incomp_prop, Model = "Incomp/Prop"),
  cbind(assessment_incomp_prop_var_select_1, Model = "Incomp/Prop/Var. Selection (1)"),
  cbind(assessment_incomp_prop_var_select, Model = "Incomp/Prop/Var. Selection (2)")
) %>%
  as.data.frame() %>%
  mutate(`5%` = as.numeric(`5%`),
         `50%` = as.numeric(`50%`),
         `95%` = as.numeric(`95%`),
         Model = factor(Model, levels = c("Incomp/Prop/Var. Selection (2)", "Incomp/Prop/Var. Selection (1)", "Incomp/Prop", "Incomp/Var. Selection (1)", "Incomp", "Incomp/Routine", "Comp/Routine/Prop", "Comp/Routine")))


plot_assessment <- assessment %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(y = Model, xmin = `5%`, xmax = `95%`, colour = Model), width = 0.2) +
  geom_point(aes(x = `50%`, y = Model, shape = Dataset), size = 2, colour = "black") +
  facet_wrap(~statistic, ncol = 1, scales = "free") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(colour = "none")



pdf(file = paste0(output_path, "/Assessment/model_residuals.pdf"))
plot_comp_routine_no_prop
plot_comp_routine_no_prop_effects_validation_dev
plot_comp_routine_no_prop_effects_validation_val
plot_comp_routine_prop
plot_comp_routine_prop_effects_validation_dev
plot_comp_routine_prop_effects_validation_val
plot_incomp_routine_no_prop
plot_incomp_routine_no_prop_effects_validation_dev
plot_incomp_routine_no_prop_effects_validation_val
plot_incomp_no_prop
plot_incomp_no_prop_effects_validation_dev
plot_incomp_no_prop_effects_validation_val
plot_incomp_no_prop_var_select
plot_incomp_no_prop_var_select_effects_validation_dev
plot_incomp_no_prop_var_select_effects_validation_val
plot_incomp_prop
plot_incomp_prop_effects_validation_dev
plot_incomp_prop_effects_validation_val
plot_incomp_prop_var_select_1
plot_incomp_prop_var_select_1_effects_validation_dev
plot_incomp_prop_var_select_1_effects_validation_val
plot_incomp_prop_var_select
plot_incomp_prop_var_select_effects_validation_dev
plot_incomp_prop_var_select_effects_validation_val
plot_assessment
dev.off()












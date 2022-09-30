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
library(bartMachine)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)


## make directory for outputs
dir.create(paste0(output_path, "/Assessment"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
######################### Read Data / Model In ################################
###############################################################################
###############################################################################

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

## plot histogram of effect

plot_effect_1 <- hist_plot(comp_routine_no_prop_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- comp_routine_no_prop_effects_summary_dev %>%
  cbind(malesex = data_complete_routine_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
          ) %>%
  filter(malesex == 1)

effects_summary_dev_female <- comp_routine_no_prop_effects_summary_dev %>%
  cbind(malesex = data_complete_routine_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


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

## plot effects validation for development
predicted_observed_dev <- data_complete_routine_dev %>%
  cbind(hba1c_diff = comp_routine_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

## plot effects validation for validation
predicted_observed_val <- data_complete_routine_val %>%
  cbind(hba1c_diff = comp_routine_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_no_prop_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_comp_routine_no_prop)


## plot histogram of effect

plot_effect_2 <- hist_plot(comp_routine_no_prop_effects_summary_val, "", -15, 20)


effects_summary_val_male <- comp_routine_no_prop_effects_summary_val %>%
  cbind(malesex = data_complete_routine_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- comp_routine_no_prop_effects_summary_val %>%
  cbind(malesex = data_complete_routine_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)



plot_comp_routine_no_prop_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Complete data (Routine/No propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)


plot_comp_routine_no_prop_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)


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


##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_comp_routine_no_prop <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_validation_dev_comp_routine_no_prop.rds"))
  
}

plot_ATE_dev_comp_routine_no_prop <- ATE_plot(ATE_validation_dev_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_comp_routine_no_prop <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_validation_val_comp_routine_no_prop.rds"))
  
}

plot_ATE_val_comp_routine_no_prop <- ATE_plot(ATE_validation_val_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_comp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_comp_routine_no_prop, plot_ATE_val_comp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_comp_routine_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_dev_comp_routine_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_comp_routine_no_prop <- ATE_plot(ATE_matching_validation_dev_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_comp_routine_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_val_comp_routine_no_prop.rds"))
  
}

plot_ATE_val_prop_score_comp_routine_no_prop <- ATE_plot(ATE_matching_validation_val_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_comp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_comp_routine_no_prop, plot_ATE_val_prop_score_comp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_comp_routine_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_comp_routine_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_weighting_comp_routine_no_prop  <- ATE_plot(ATE_weighting_validation_dev_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_comp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_comp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_comp_routine_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_comp_routine_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_val_comp_routine_no_prop.rds"))
  
}

plot_ATE_val_prop_score_weighting_comp_routine_no_prop  <- ATE_plot(ATE_weighting_validation_val_comp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_comp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_comp_routine_no_prop, plot_ATE_val_prop_score_weighting_comp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))




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

## plot histogram of effect

plot_effect_1 <- hist_plot(comp_routine_prop_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- comp_routine_prop_effects_summary_dev %>%
  cbind(malesex = data_complete_routine_prop_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- comp_routine_prop_effects_summary_dev %>%
  cbind(malesex = data_complete_routine_prop_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)



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

## plot effects validation for development
predicted_observed_dev <- data_complete_routine_prop_dev %>%
  cbind(hba1c_diff = comp_routine_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

## plot effects validation for validation
predicted_observed_val <- data_complete_routine_prop_val %>%
  cbind(hba1c_diff = comp_routine_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_comp_routine_prop_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_comp_routine_prop_model)


## plot histogram of effect

plot_effect_2 <- hist_plot(comp_routine_prop_effects_summary_val, "", -15, 20)


effects_summary_val_male <- comp_routine_prop_effects_summary_val %>%
  cbind(malesex = data_complete_routine_prop_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- comp_routine_prop_effects_summary_val %>%
  cbind(malesex = data_complete_routine_prop_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)


plot_comp_routine_prop_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Complete data (Routine/Propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)


plot_comp_routine_prop_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)



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



##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_comp_routine_prop <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_comp_routine_prop, paste0(output_path, "/Assessment/ATE_validation_dev_comp_routine_prop.rds"))
  
}

plot_ATE_dev_comp_routine_prop <- ATE_plot(ATE_validation_dev_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_comp_routine_prop <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_comp_routine_prop, paste0(output_path, "/Assessment/ATE_validation_val_comp_routine_prop.rds"))
  
}

plot_ATE_val_comp_routine_prop <- ATE_plot(ATE_validation_val_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_comp_routine_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_comp_routine_prop, plot_ATE_val_comp_routine_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_comp_routine_prop <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_comp_routine_prop, paste0(output_path, "/Assessment/ATE_matching_validation_dev_comp_routine_prop.rds"))
  
}

plot_ATE_dev_prop_score_comp_routine_prop <- ATE_plot(ATE_matching_validation_dev_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_comp_routine_prop <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_comp_routine_prop, paste0(output_path, "/Assessment/ATE_matching_validation_val_comp_routine_prop.rds"))
  
}

plot_ATE_val_prop_score_comp_routine_prop <- ATE_plot(ATE_matching_validation_val_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_comp_routine_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_comp_routine_prop, plot_ATE_val_prop_score_comp_routine_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_comp_routine_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_comp_routine_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_comp_routine_prop.rds"))
  
}

plot_ATE_dev_prop_score_weighting_comp_routine_prop  <- ATE_plot(ATE_weighting_validation_dev_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_comp_routine_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_comp_routine_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_comp_routine_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_comp_routine_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_val_comp_routine_prop.rds"))
  
}

plot_ATE_val_prop_score_weighting_comp_routine_prop  <- ATE_plot(ATE_weighting_validation_val_comp_routine_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_comp_routine_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_comp_routine_prop, plot_ATE_val_prop_score_weighting_comp_routine_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_routine_no_prop_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_routine_no_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_routine_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_routine_no_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_routine_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_routine_dev %>%
  cbind(hba1c_diff = incomp_routine_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


## plot effects validation for validation
predicted_observed_val <- data_incomplete_routine_val %>%
  cbind(hba1c_diff = incomp_routine_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_routine_no_prop_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_routine_no_prop)

## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_routine_no_prop_effects_summary_val, "", -15, 20)


effects_summary_val_male <- incomp_routine_no_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_routine_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_routine_no_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_routine_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)



plot_incomp_routine_no_prop_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Incomplete data (Routine/No propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)

plot_incomp_routine_no_prop_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)


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



##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_routine_no_prop <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_routine_no_prop.rds"))
  
}

plot_ATE_dev_incomp_routine_no_prop <- ATE_plot(ATE_validation_dev_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_routine_no_prop <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_validation_val_incomp_routine_no_prop.rds"))
  
}

plot_ATE_val_incomp_routine_no_prop <- ATE_plot(ATE_validation_val_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_routine_no_prop, plot_ATE_val_incomp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_routine_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_routine_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_incomp_routine_no_prop <- ATE_plot(ATE_matching_validation_dev_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_routine_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_routine_no_prop.rds"))
  
}

plot_ATE_val_prop_score_incomp_routine_no_prop <- ATE_plot(ATE_matching_validation_val_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_routine_no_prop, plot_ATE_val_prop_score_incomp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_routine_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_routine_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_routine_no_prop  <- ATE_plot(ATE_weighting_validation_dev_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_routine_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_routine_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_routine_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_routine_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_routine_no_prop.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_routine_no_prop  <- ATE_plot(ATE_weighting_validation_val_incomp_routine_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_routine_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_routine_no_prop, plot_ATE_val_prop_score_weighting_incomp_routine_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_no_prop_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_no_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_no_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)

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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_dev %>%
  cbind(hba1c_diff = incomp_no_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


## plot effects validation for validation
predicted_observed_val <- data_incomplete_val %>%
  cbind(hba1c_diff = incomp_no_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_no_prop)

## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_no_prop_effects_summary_val, "", -15, 20)


effects_summary_val_male <- incomp_no_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_no_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)




plot_incomp_no_prop_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Incomplete data (No propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)


plot_incomp_no_prop_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)


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


##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_no_prop <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_no_prop, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_no_prop.rds"))
  
}

plot_ATE_dev_incomp_no_prop <- ATE_plot(ATE_validation_dev_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_no_prop <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_no_prop, paste0(output_path, "/Assessment/ATE_validation_val_incomp_no_prop.rds"))
  
}

plot_ATE_val_incomp_no_prop <- ATE_plot(ATE_validation_val_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_no_prop, plot_ATE_val_incomp_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_incomp_no_prop <- ATE_plot(ATE_matching_validation_dev_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_no_prop <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_no_prop, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_no_prop.rds"))
  
}

plot_ATE_val_prop_score_incomp_no_prop <- ATE_plot(ATE_matching_validation_val_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_no_prop, plot_ATE_val_prop_score_incomp_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_no_prop.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_no_prop  <- ATE_plot(ATE_weighting_validation_dev_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_no_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_no_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_no_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_no_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_no_prop.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_no_prop  <- ATE_plot(ATE_weighting_validation_val_incomp_no_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_no_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_no_prop, plot_ATE_val_prop_score_weighting_incomp_no_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))



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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_no_prop_var_select_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_no_prop_var_select_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_no_prop_var_select_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)

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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_dev_var_select %>%
  cbind(hba1c_diff = incomp_no_prop_var_select_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

## plot effects validation for validation
predicted_observed_val <- data_incomplete_val_var_select %>%
  cbind(hba1c_diff = incomp_no_prop_var_select_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_no_prop_var_select_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_no_prop_var_select)


## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_no_prop_var_select_effects_summary_val, "", -15, 20)


effects_summary_val_male <- incomp_no_prop_var_select_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_no_prop_var_select_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)



plot_incomp_no_prop_var_select_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Variable Selection 1, Incomplete data (No propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)

plot_incomp_no_prop_var_select_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)



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

##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_no_prop_var_select <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_dev_incomp_no_prop_var_select <- ATE_plot(ATE_validation_dev_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_no_prop_var_select <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_validation_val_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_val_incomp_no_prop_var_select <- ATE_plot(ATE_validation_val_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_no_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_no_prop_var_select, plot_ATE_val_incomp_no_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_no_prop_var_select <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_dev_prop_score_incomp_no_prop_var_select <- ATE_plot(ATE_matching_validation_dev_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_no_prop_var_select <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_val_prop_score_incomp_no_prop_var_select <- ATE_plot(ATE_matching_validation_val_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_no_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_no_prop_var_select, plot_ATE_val_prop_score_incomp_no_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_no_prop_var_select <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_no_prop_var_select  <- ATE_plot(ATE_weighting_validation_dev_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_no_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_no_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_no_prop_var_select <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_no_prop_var_select, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_no_prop_var_select.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_no_prop_var_select  <- ATE_plot(ATE_weighting_validation_val_incomp_no_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_no_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_no_prop_var_select, plot_ATE_val_prop_score_weighting_incomp_no_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_prop_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_prop_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_dev %>%
  cbind(hba1c_diff = incomp_prop_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


## plot effects validation for validation
predicted_observed_val <- data_incomplete_val %>%
  cbind(hba1c_diff = incomp_prop_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_prop_model)


## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_prop_effects_summary_val, "", -15, 20)


effects_summary_val_male <- incomp_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_prop_effects_summary_val %>%
  cbind(malesex = data_incomplete_val %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)



plot_incomp_prop_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Incomplete data (Propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)

plot_incomp_prop_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)


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


##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_prop <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_prop, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop.rds"))
  
}

plot_ATE_dev_incomp_prop <- ATE_plot(ATE_validation_dev_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_prop <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_prop, paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop.rds"))
  
}

plot_ATE_val_incomp_prop <- ATE_plot(ATE_validation_val_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_prop, plot_ATE_val_incomp_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_prop <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_prop, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop.rds"))
  
}

plot_ATE_dev_prop_score_incomp_prop <- ATE_plot(ATE_matching_validation_dev_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_prop <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_prop, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop.rds"))
  
}

plot_ATE_val_prop_score_incomp_prop <- ATE_plot(ATE_matching_validation_val_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_prop, plot_ATE_val_prop_score_incomp_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_prop  <- ATE_plot(ATE_weighting_validation_dev_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_prop <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_prop <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_prop, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_prop  <- ATE_plot(ATE_weighting_validation_val_incomp_prop[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_prop <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_prop, plot_ATE_val_prop_score_weighting_incomp_prop, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_prop_var_select_1_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_prop_var_select_1_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select_1 %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_prop_var_select_1_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select_1 %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_dev_var_select_1 %>%
  cbind(hba1c_diff = incomp_prop_var_select_1_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


## plot effects validation for validation
predicted_observed_val <- data_incomplete_val_var_select_1 %>%
  cbind(hba1c_diff = incomp_prop_var_select_1_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_1_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_prop_model_var_select_1)

## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_prop_var_select_1_effects_summary_val, "", -15, 20)


effects_summary_val_male <- incomp_prop_var_select_1_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select_1 %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_prop_var_select_1_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select_1 %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)



plot_incomp_prop_var_select_1_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Variable Selection 1, Incomplete data (Propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)


plot_incomp_prop_var_select_1_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)



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


##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_prop_var_select_1 <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_dev_incomp_prop_var_select_1 <- ATE_plot(ATE_validation_dev_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_prop_var_select_1 <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_val_incomp_prop_var_select_1 <- ATE_plot(ATE_validation_val_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_prop_var_select_1 <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_prop_var_select_1, plot_ATE_val_incomp_prop_var_select_1, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_prop_var_select_1 <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_dev_prop_score_incomp_prop_var_select_1 <- ATE_plot(ATE_matching_validation_dev_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_prop_var_select_1 <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_val_prop_score_incomp_prop_var_select_1 <- ATE_plot(ATE_matching_validation_val_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_prop_var_select_1 <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_prop_var_select_1, plot_ATE_val_prop_score_incomp_prop_var_select_1, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_prop_var_select_1 <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_prop_var_select_1  <- ATE_plot(ATE_weighting_validation_dev_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_prop_var_select_1 <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop_var_select_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_prop_var_select_1 <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_prop_var_select_1, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop_var_select_1.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_prop_var_select_1  <- ATE_plot(ATE_weighting_validation_val_incomp_prop_var_select_1[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_prop_var_select_1 <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_prop_var_select_1, plot_ATE_val_prop_score_weighting_incomp_prop_var_select_1, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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

## plot histogram of effect

plot_effect_1 <- hist_plot(incomp_prop_var_select_effects_summary_dev, "", -15, 20)


effects_summary_dev_male <- incomp_prop_var_select_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_dev_female <- incomp_prop_var_select_effects_summary_dev %>%
  cbind(malesex = data_incomplete_dev_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.dev %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)



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

## plot effects validation for development
predicted_observed_dev <- data_incomplete_dev_var_select %>%
  cbind(hba1c_diff = incomp_prop_var_select_effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


## plot effects validation for validation
predicted_observed_val <- data_incomplete_val_var_select %>%
  cbind(hba1c_diff = incomp_prop_var_select_effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


plot_incomp_prop_var_select_effects_validation <- plot_full_effects_validation(predicted_observed_dev, predicted_observed_val, bart_incomp_prop_model_var_select)

## plot histogram of effect

plot_effect_2 <- hist_plot(incomp_prop_var_select_effects_summary_val, "", -15, 20)



effects_summary_val_male <- incomp_prop_var_select_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 1)


effects_summary_val_female <- incomp_prop_var_select_effects_summary_val %>%
  cbind(malesex = data_incomplete_val_var_select %>%
          select(patid, pateddrug) %>%
          left_join(final.val %>%
                      select(patid, pateddrug, malesex), by = c("patid", "pateddrug")) %>%
          select(malesex)
  ) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)




plot_incomp_prop_var_select_effects <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model fitting: Variable Selection 2, Incomplete data (Propensity score)")
  
  ,
  
  cowplot::plot_grid(
    
    plot_effect_1, 
    
    plot_effect_2
    
    , ncol = 2, nrow = 1, labels = c("A", "B")
    
  ), ncol = 1, nrow = 2, rel_heights = c(0.1,1)
  
)


plot_incomp_prop_var_select_effects_genders <- cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)




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



##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev_incomp_prop_var_select <- calc_ATE_validation(predicted_observed_dev)
  
  saveRDS(ATE_validation_dev_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_validation_dev_incomp_prop_var_select.rds"))
  
}

plot_ATE_dev_incomp_prop_var_select <- ATE_plot(ATE_validation_dev_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)


if (class(try(
  
  ATE_validation_val_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val_incomp_prop_var_select <- calc_ATE_validation(predicted_observed_val)
  
  saveRDS(ATE_validation_val_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_validation_val_incomp_prop_var_select.rds"))
  
}

plot_ATE_val_incomp_prop_var_select <- ATE_plot(ATE_validation_val_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12)

plot_ATE_incomp_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_incomp_prop_var_select, plot_ATE_val_incomp_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev_incomp_prop_var_select <- calc_ATE_validation_prop_matching(predicted_observed_dev)
  
  saveRDS(ATE_matching_validation_dev_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_matching_validation_dev_incomp_prop_var_select.rds"))
  
}

plot_ATE_dev_prop_score_incomp_prop_var_select <- ATE_plot(ATE_matching_validation_dev_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_matching_validation_val_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val_incomp_prop_var_select <- calc_ATE_validation_prop_matching(predicted_observed_val)
  
  saveRDS(ATE_matching_validation_val_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_matching_validation_val_incomp_prop_var_select.rds"))
  
}

plot_ATE_val_prop_score_incomp_prop_var_select <- ATE_plot(ATE_matching_validation_val_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching_incomp_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_incomp_prop_var_select, plot_ATE_val_prop_score_incomp_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev_incomp_prop_var_select <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev)
  
  saveRDS(ATE_weighting_validation_dev_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_weighting_validation_dev_incomp_prop_var_select.rds"))
  
}

plot_ATE_dev_prop_score_weighting_incomp_prop_var_select  <- ATE_plot(ATE_weighting_validation_dev_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val_incomp_prop_var_select <- readRDS(paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop_var_select.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val_incomp_prop_var_select <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val)
  
  saveRDS(ATE_weighting_validation_val_incomp_prop_var_select, paste0(output_path, "/Assessment/ATE_weighting_validation_val_incomp_prop_var_select.rds"))
  
}

plot_ATE_val_prop_score_weighting_incomp_prop_var_select  <- ATE_plot(ATE_weighting_validation_val_incomp_prop_var_select[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting_incomp_prop_var_select <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting_incomp_prop_var_select, plot_ATE_val_prop_score_weighting_incomp_prop_var_select, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


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



pdf(file = "Plots/6.1.model_residuals.pdf")
plot_comp_routine_no_prop
plot_comp_routine_prop
plot_incomp_routine_no_prop
plot_incomp_no_prop
plot_incomp_no_prop_var_select
plot_incomp_prop
plot_incomp_prop_var_select_1
plot_incomp_prop_var_select
plot_assessment
dev.off()


pdf(file = "Plots/6.1.model_effects.pdf")
plot_comp_routine_no_prop_effects
plot_comp_routine_no_prop_effects_genders
plot_comp_routine_no_prop_effects_validation
plot_ATE_comp_routine_no_prop
plot_ATE_prop_score_matching_comp_routine_no_prop
plot_ATE_prop_score_weighting_comp_routine_no_prop
plot_comp_routine_prop_effects
plot_comp_routine_prop_effects_genders
plot_comp_routine_prop_effects_validation
plot_ATE_comp_routine_prop
plot_ATE_prop_score_matching_comp_routine_prop
plot_ATE_prop_score_weighting_comp_routine_prop
plot_incomp_routine_no_prop_effects
plot_incomp_routine_no_prop_effects_genders
plot_incomp_routine_no_prop_effects_validation
plot_ATE_incomp_routine_no_prop
plot_ATE_prop_score_matching_incomp_routine_no_prop
plot_ATE_prop_score_weighting_incomp_routine_no_prop
plot_incomp_no_prop_effects
plot_incomp_no_prop_effects_genders
plot_incomp_no_prop_effects_validation
plot_ATE_incomp_no_prop
plot_ATE_prop_score_matching_incomp_no_prop
plot_ATE_prop_score_weighting_incomp_no_prop
plot_incomp_no_prop_var_select_effects
plot_incomp_no_prop_var_select_effects_genders
plot_incomp_no_prop_var_select_effects_validation
plot_ATE_incomp_no_prop_var_select
plot_ATE_prop_score_matching_incomp_no_prop_var_select
plot_ATE_prop_score_weighting_incomp_no_prop_var_select
plot_incomp_prop_effects
plot_incomp_prop_effects_genders
plot_incomp_prop_effects_validation
plot_ATE_incomp_prop
plot_ATE_prop_score_matching_incomp_prop
plot_ATE_prop_score_weighting_incomp_prop
plot_incomp_prop_var_select_1_effects
plot_incomp_prop_var_select_1_effects_genders
plot_incomp_prop_var_select_1_effects_validation
plot_ATE_incomp_prop_var_select_1
plot_ATE_prop_score_matching_incomp_prop_var_select_1
plot_ATE_prop_score_weighting_incomp_prop_var_select_1
plot_incomp_prop_var_select_effects
plot_incomp_prop_var_select_effects_genders
plot_incomp_prop_var_select_effects_validation
plot_ATE_incomp_prop_var_select
plot_ATE_prop_score_matching_incomp_prop_var_select
plot_ATE_prop_score_weighting_incomp_prop_var_select
dev.off()










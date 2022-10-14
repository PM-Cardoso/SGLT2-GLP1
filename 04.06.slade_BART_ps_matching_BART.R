####################
## Description:
##  - In this file we:
##    - Fit a BART PS model to all variables.
##    - Perform variable selection on PS model and refit model.
##    - Match individuals according to their propensity score.
##    - Fit BART model with all variables to matched individuals.
##    - Perform variable selection on BART model
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
dir.create(paste0(output_path, "/Final_model/ps_matching"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/ps_matching/Assessment"))

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
  
  bart_ps_model <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model <- bartMachine::bartMachine(X = dataset.dev %>%
                                              select(-patid,
                                                     -pateddrug,
                                                     -bothdrugs,
                                                     -posthba1c_final,
                                                     -drugclass,
                                                     -yrdrugstart,
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
  
  saveRDS(bart_ps_model, paste0(output_path, "/Final_model/ps_matching/bart_ps_model.rds"))
  
}

########
### Variable selection of propensity score model
########

if (class(try(
  
  vs_bart_ps_model <- readRDS(paste0(output_path, "/Final_model/ps_matching/vs_bart_ps_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/4.6.prop_model_vs.pdf", width = 18, height = 11)
  # error with cv
  vs_bart_ps_model <- var_selection_by_permute(bart_ps_model)
  dev.off()
  
  ## Long version of the var selection
  # [1] "prebmi"                         "t2dmduration"
  # [3] "egfr_ckdepi"                    "drugline_5"
  # [5] "drugline_2"                     "prealt"
  # [7] "prehdl"                         "ncurrtx_3"
  # [9] "score.excl.mi"                  "drugline_4"
  # [11] "Category_Non-smoker"            "predrug.5yrrecent.neuropathy_1"
  # [13] "drugline_3"
  
  saveRDS(vs_bart_ps_model, paste0(output_path, "/Final_model/ps_matching/vs_bart_ps_model.rds"))
  
}


########
### Fit a propensity score model to all variables in dataset (include yrdrugstart)
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_ps_model_yrdrugstart <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_ps_model_yrdrugstart.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model_yrdrugstart <- bartMachine::bartMachine(X = dataset.dev %>%
                                              select(-patid,
                                                     -pateddrug,
                                                     -bothdrugs,
                                                     -posthba1c_final,
                                                     -drugclass,
                                                     # -yrdrugstart,
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
  
  saveRDS(bart_ps_model_yrdrugstart, paste0(output_path, "/Final_model/ps_matching/bart_ps_model_yrdrugstart.rds"))
  
}

########
### Variable selection of propensity score model
########

if (class(try(
  
  vs_bart_ps_model_yrdrugstart <- readRDS(paste0(output_path, "/Final_model/ps_matching/vs_bart_ps_model_yrdrugstart.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/4.6.prop_model_vs_yrdrugstart.pdf", width = 18, height = 11)
  # error with cv
  vs_bart_ps_model_yrdrugstart <- var_selection_by_permute(bart_ps_model_yrdrugstart)
  dev.off()
  
  ## Long version of the var selection
  
  saveRDS(vs_bart_ps_model_yrdrugstart, paste0(output_path, "/Final_model/ps_matching/vs_bart_ps_model_yrdrugstart.rds"))
  
}



########
### Refit PS model with selected vars
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_ps_model_final_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_ps_model_final_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_ps_model_final_dev <- bartMachine::bartMachine(X = dataset.dev %>%
                                                    select(
                                                      ### yrdrugstart included
                                                      # yrdrugstart,
                                                      # prebmi,
                                                      # t2dmduration,
                                                      # drugline,
                                                      # prehba1cmmol,
                                                      # ncurrtx,
                                                      # Category
                                                      ### no yrdrugstart included
                                                      prebmi,
                                                      t2dmduration,
                                                      prealt,
                                                      prehdl,
                                                      egfr_ckdepi,
                                                      drugline,
                                                      ncurrtx,
                                                      score.excl.mi,
                                                      Category,
                                                      predrug.5yrrecent.neuropathy
                                                      ),
                                            y = dataset.dev[,"drugclass"],
                                            use_missing_data = TRUE,
                                            impute_missingness_with_rf_impute = FALSE,
                                            impute_missingness_with_x_j_bar_for_lm = TRUE,
                                            num_trees = 200,
                                            num_burn_in = 3000,
                                            num_iterations_after_burn_in = 1000,
                                            serialize = TRUE)
  
  saveRDS(bart_ps_model_final_dev, paste0(output_path, "/Final_model/ps_matching/bart_ps_model_final_dev.rds"))
  
}


########
### Match patients with propensity scores
##    - Because we have more SGLT2 than GLP1, we match patients to GLP1
########

prop_score = 1 - bart_ps_model_final_dev$p_hat_train

# dataset individuals that had GLP1
rows.glp1 <- which(dataset.dev$drugclass == "GLP1")

# dataset individuals that had SGLT2
rows.sglt2 <- which(dataset.dev$drugclass == "SGLT2")

# list of matched SGLT2 rows to GLP1 (the maximum number of patients is equal to GLP1 number of patients)
matched.sglt2 <- vector(mode = "numeric", length = length(rows.glp1))

# iterate through rows of GLP1
for (l in 1:length(rows.glp1)) {
  # closest SGLT2 row to GLP1
  chosen.row <- which.min(abs(prop_score[rows.sglt2] - prop_score[rows.glp1[l]]))
  
  # check if distance is less than 0.05 (caliper distance)
  if (prop_score[rows.sglt2[chosen.row]] - prop_score[rows.glp1[l]] < 0.05) {
    # if chosen row is within caliper distance
    
    # update list of matched rows
    matched.sglt2[l] <- rows.sglt2[chosen.row]
    
    # remove row from being matched again
    rows.sglt2 <- rows.sglt2[-chosen.row]
    
  } else {
    # if chosen row is outside caliper distance
    
    # update list of matched rows with NA
    matched.sglt2[l] <- NA
    
  }
}

# Combined patients that have been matched
dataset.dev.matched <- rbind(
  # SGLT2 patients that have been matched
  dataset.dev[matched.sglt2[!is.na(matched.sglt2)],],
  # GLP1 patients that have been matched
  dataset.dev[rows.glp1[!is.na(matched.sglt2)],]
)


########
### Compare characteristics of patients between both drugs
########

# seems fine
psych::describeBy(dataset.dev.matched, group = dataset.dev.matched$drugclass, type = 1)


########
### Fit BART model to all variables
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(
  
  bart_model <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model <- bartMachine::bartMachine(X = dataset.dev.matched %>%
                                           select(-patid,
                                                  -pateddrug,
                                                  -bothdrugs,
                                                  -posthba1c_final,
                                                  -sglt2subtype,
                                                  -glp1subtype),
                                         y = dataset.dev.matched[,"posthba1c_final"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)
  
  saveRDS(bart_model, paste0(output_path, "/Final_model/ps_matching/bart_model.rds"))
  
}

########
### Variable selection of BART model
########

# One-off variable selection
if (class(try(
  
  vs_bart_model <- readRDS(paste0(output_path, "/Final_model/ps_matching/vs_bart_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = "Plots/4.6.response_model_vs.pdf", width = 18, height = 11)
  vs_bart_model <- var_selection_by_permute(bart_model)
  dev.off()
  
  ## Long version of the var selection
  
  saveRDS(vs_bart_model, paste0(output_path, "/Final_model/ps_matching/vs_bart_model.rds"))
  
}

# Cross-validation selection
if (class(try(
  
  vs_bart_model_cv <- readRDS(paste0(output_path, "/Final_model/ps_matching/vs_bart_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  vs_bart_model_cv <- var_selection_by_permute_cv(bart_model)
  
  ## Long version of the var selection
  
  saveRDS(vs_bart_model_cv, paste0(output_path, "/Final_model/ps_matching/vs_bart_model_cv.rds"))
  
}


########
### Refit BART model with selected variables
########

## Fit final model using selected variables to estimate HbA1c outcome
if (class(try(

  bart_model_final <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_model_final.rds"))

  , silent = TRUE)) == "try-error") {

  bart_model_final <- bartMachine::bartMachine(X = dataset.dev.matched %>%
                                           select(
                                             prehba1cmmol,
                                             prebil,
                                             drugclass,
                                             yrdrugstart,
                                             hba1cmonth,
                                             score.excl.mi,
                                             drugline,
                                             malesex,
                                             ncurrtx
                                            ),
                                         y = dataset.dev.matched[,"posthba1c_final"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)

  saveRDS(bart_model_final, paste0(output_path, "/Final_model/ps_matching/bart_model_final.rds"))

}



########
### Validation of model
########

# Dev

data_dev <- dataset.dev.matched %>%
  select(c(patid, pateddrug, posthba1c_final,
           colnames(bart_model_final$X)))


## Get posteriors
if (class(try(

  posteriors_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/posteriors_dev.rds"))

  , silent = TRUE)) == "try-error") {

  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_final, data_dev %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/posteriors_dev.rds"))


}

### residuals calculation
if (class(try(

  cred_pred_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/cred_pred_dev.rds"))

  , silent = TRUE)) == "try-error") {

  cred_pred_dev <- calc_resid(data_dev, posteriors_dev, "posthba1c_final")

  saveRDS(cred_pred_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/cred_pred_dev.rds"))

}


# Val

dataset.val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi))


## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(

  bart_ps_model_final_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/bart_ps_model_final_val.rds"))

  , silent = TRUE)) == "try-error") {

  bart_ps_model_final_val <- bartMachine::bartMachine(X = dataset.val %>%
                                                        select(
                                                          ### yrdrugstart included
                                                          yrdrugstart,
                                                          prebmi,
                                                          t2dmduration,
                                                          drugline,
                                                          prehba1cmmol,
                                                          ncurrtx,
                                                          Category
                                                          ### no yrdrugstart included
                                                          # prebmi,
                                                          # t2dmduration,
                                                          # prealt,
                                                          # prehdl,
                                                          # egfr_ckdepi,
                                                          # drugline,
                                                          # ncurrtx,
                                                          # score.excl.mi,
                                                          # Category,
                                                          # predrug.5yrrecent.neuropathy
                                                        ),
                                                      y = dataset.val[,"drugclass"],
                                                      use_missing_data = TRUE,
                                                      impute_missingness_with_rf_impute = FALSE,
                                                      impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                      num_trees = 200,
                                                      num_burn_in = 3000,
                                                      num_iterations_after_burn_in = 1000,
                                                      serialize = TRUE)

  saveRDS(bart_ps_model_final_val, paste0(output_path, "/Final_model/ps_matching/bart_ps_model_final_val.rds"))

}

########
### Match patients with propensity scores
##    - Because we have more SGLT2 than GLP1, we match patients to GLP1
########

prop_score = 1 - bart_ps_model_final_val$p_hat_train

# dataset individuals that had GLP1
rows.glp1 <- which(dataset.val$drugclass == "GLP1")

# dataset individuals that had SGLT2
rows.sglt2 <- which(dataset.val$drugclass == "SGLT2")

# list of matched SGLT2 rows to GLP1 (the maximum number of patients is equal to GLP1 number of patients)
matched.sglt2 <- vector(mode = "numeric", length = length(rows.glp1))

# iterate through rows of GLP1
for (l in 1:length(rows.glp1)) {
  # closest SGLT2 row to GLP1
  chosen.row <- which.min(abs(prop_score[rows.sglt2] - prop_score[rows.glp1[l]]))

  # check if distance is less than 0.05 (caliper distance)
  if (prop_score[rows.sglt2[chosen.row]] - prop_score[rows.glp1[l]] < 0.05) {
    # if chosen row is within caliper distance

    # update list of matched rows
    matched.sglt2[l] <- rows.sglt2[chosen.row]

    # remove row from being matched again
    rows.sglt2 <- rows.sglt2[-chosen.row]

  } else {
    # if chosen row is outside caliper distance

    # update list of matched rows with NA
    matched.sglt2[l] <- NA

  }
}

# Combined patients that have been matched
dataset.val.matched <- rbind(
  # SGLT2 patients that have been matched
  dataset.val[matched.sglt2[!is.na(matched.sglt2)],],
  # GLP1 patients that have been matched
  dataset.val[rows.glp1[!is.na(matched.sglt2)],]
)


data_val <- dataset.val.matched %>%
  select(c(patid, pateddrug, posthba1c_final,
           colnames(bart_model_final$X)))


### Get posteriors
if (class(try(

  posteriors_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/posteriors_val.rds"))

  , silent = TRUE)) == "try-error") {

  posteriors_val <- bartMachine::bart_machine_get_posterior(bart_model_final, data_val %>%
                                                              select(
                                                                colnames(bart_model_final$X)
                                                              ))
  saveRDS(posteriors_val, paste0(output_path, "/Final_model/ps_matching/Assessment/posteriors_val.rds"))


}

### residuals calculation
if (class(try(

  cred_pred_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/cred_pred_val.rds"))

  , silent = TRUE)) == "try-error") {

  cred_pred_val <- calc_resid(data_val, posteriors_val, "posthba1c_final")

  saveRDS(cred_pred_val, paste0(output_path, "/Final_model/ps_matching/Assessment/cred_pred_val.rds"))

}



# assessment of R2, RSS, RMSE
if (class(try(

  assessment <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/assessment.rds"))

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

  saveRDS(assessment, paste0(output_path, "/Final_model/ps_matching/Assessment/assessment.rds"))

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





plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Residuals of Model 6")




########
### Validation of effect
########

# calculate effects
if (class(try(

  effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/effects_summary_dev.rds"))

  , silent = TRUE)) == "try-error") {

  effects_summary_dev <- calc_effect_summary(bart_model_final, data_dev)

  saveRDS(effects_summary_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/effects_summary_dev.rds"))

}


# calculate effects
if (class(try(

  effects_summary_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/effects_summary_val.rds"))

  , silent = TRUE)) == "try-error") {

  effects_summary_val <- calc_effect_summary(bart_model_final, data_val)

  saveRDS(effects_summary_val, paste0(output_path, "/Final_model/ps_matching/Assessment/effects_summary_val.rds"))

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


effects_summary_dev_male <- effects_summary_dev %>%
  cbind(malesex = data_dev$malesex) %>%
  filter(malesex == 1)

effects_summary_dev_female <- effects_summary_dev %>%
  cbind(malesex = data_dev$malesex) %>%
  filter(malesex == 0)


plot_effect_1_male <- hist_plot(effects_summary_dev_male, "Male", -15, 20)

plot_effect_1_female <- hist_plot(effects_summary_dev_female, "Female", -15, 20)


effects_summary_val_male <- effects_summary_val %>%
  cbind(malesex = data_val$malesex) %>%
  filter(malesex == 1)

effects_summary_val_female <- effects_summary_val %>%
  cbind(malesex = data_val$malesex) %>%
  filter(malesex == 0)

plot_effect_2_male <- hist_plot(effects_summary_val_male, "Male", -15, 20)

plot_effect_2_female <- hist_plot(effects_summary_val_female, "Female", -15, 20)


##############
# Validating ATE
if (class(try(

  ATE_validation_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_validation_dev <- calc_ATE_validation(predicted_observed_dev, "posthba1c_final")

  saveRDS(ATE_validation_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_validation_dev.rds"))

}

plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)


if (class(try(

  ATE_validation_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_validation_val <- calc_ATE_validation(predicted_observed_val, "posthba1c_final")

  saveRDS(ATE_validation_val, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_validation_val.rds"))

}

plot_ATE_val <- ATE_plot(ATE_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(hba1c~drugclass+prop_score)")

  ,

  cowplot::plot_grid(plot_ATE_dev, plot_ATE_val, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(

  ATE_matching_validation_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_matching_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_matching_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_dev, "posthba1c_final")

  saveRDS(ATE_matching_validation_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_matching_validation_dev.rds"))

}

plot_ATE_dev_prop_score <- ATE_plot(ATE_matching_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(

  ATE_matching_validation_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_matching_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_matching_validation_val <- calc_ATE_validation_prop_matching(predicted_observed_val, "posthba1c_final")

  saveRDS(ATE_matching_validation_val, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_matching_validation_val.rds"))

}

plot_ATE_val_prop_score <- ATE_plot(ATE_matching_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")

  ,

  cowplot::plot_grid(plot_ATE_dev_prop_score, plot_ATE_val_prop_score, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(

  ATE_weighting_validation_dev <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_weighting_validation_dev.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "posthba1c_final")

  saveRDS(ATE_weighting_validation_dev, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_weighting_validation_dev.rds"))

}

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(

  ATE_weighting_validation_val <- readRDS(paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_weighting_validation_val.rds"))

  , silent = TRUE)) == "try-error") {

  ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1c_final")

  saveRDS(ATE_weighting_validation_val, paste0(output_path, "/Final_model/ps_matching/Assessment/ATE_weighting_validation_val.rds"))

}

plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_weighting <- cowplot::plot_grid(

  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")

  ,

  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting, plot_ATE_val_prop_score_weighting, ncol = 2, nrow = 1, labels = c("A", "B"))

  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))





#### PDF with all the plots


pdf(file = "Plots/4.6.model6_plots.pdf")
plot_residuals
plot_assessment
plot_effects_validation
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

pdf(file = "Plots/4.6.model6_partial_dependence.pdf")
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












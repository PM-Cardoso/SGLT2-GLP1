####################
## Description:
##  - In this file we fit the predicted treatment effects against the 
##      variables used in the model.
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")

library(bartMachine)
library(tidyverse)


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
dir.create(paste0(output_path, "/Final_model/model_7/Treatment_Effects"))


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
              select(patid, pateddrug, score.excl.mi))



# load in original BART model

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/model_7/bart_model_final.rds")

# load in treatment effects

effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/model_7/Assessment/effects_summary_dev.rds")


# Fit new treatment effects model

if (class(try(
  
  bart_model_effects <- readRDS(paste0(output_path, "/Final_model/model_7/Treatment_Effects/bart_model_effects.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model_effects <- bartMachine::bartMachine(X = dataset.dev %>%
                                                 select(colnames(bart_model_final$X)
                                                 ) %>%
                                                   select(-drugclass),
                                               y = effects_summary_dev[,"mean"],
                                               use_missing_data = TRUE,
                                               impute_missingness_with_rf_impute = FALSE,
                                               impute_missingness_with_x_j_bar_for_lm = TRUE,
                                               num_trees = 200,
                                               num_burn_in = 3000,
                                               num_iterations_after_burn_in = 1000,
                                               serialize = TRUE)
  
  saveRDS(bart_model_effects, paste0(output_path, "/Final_model/model_7/Treatment_Effects/bart_model_effects.rds"))
  
}


########
### Validation of model
########

data_dev <- dataset.dev %>%
  cbind(effects = effects_summary_dev$mean) %>%
  select(c(patid, pateddrug, effects, 
           colnames(bart_model_effects$X)))


## Get posteriors
if (class(try(
  
  posteriors_dev <- readRDS(paste0(output_path, "/Final_model/model_7/Treatment_Effects/posteriors_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_effects, data_dev %>%
                                                              select(
                                                                colnames(bart_model_effects$X)
                                                              ))
  saveRDS(posteriors_dev, paste0(output_path, "/Final_model/model_7/Treatment_Effects/posteriors_dev.rds"))
  
  
}


### residuals calculation
if (class(try(
  
  cred_pred_dev <- readRDS(paste0(output_path, "/Final_model/model_7/Treatment_Effects/cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_dev <- calc_resid(data_dev, posteriors_dev, "effects")
  
  saveRDS(cred_pred_dev, paste0(output_path, "/Final_model/model_7/Treatment_Effects/cred_pred_dev.rds"))
  
}



# Plot of predicted vs observed for development dataset
plot_pred_obs <- cred_pred_dev %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
  geom_point(aes(x = orig, y = mean)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
  xlim(min(cred_pred_dev$orig), max(cred_pred_dev$orig)) +
  ylim(min(cred_pred_dev$orig), max(cred_pred_dev$orig)) +
  xlab("Observed HbA1c (mmol/mol)") +
  ylab("Predicted HbA1c (mmol/mol)")

# Plot of standardised residuals for development dataset
plot_resid <- cred_pred_dev %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
  geom_point(aes(x = mean, y = std.resid)) +
  geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
  stat_smooth(aes(x = mean, y = std.resid)) +
  xlim(min(cred_pred_dev$mean), max(cred_pred_dev$mean)) +
  ylim(min(cred_pred_dev$std.resid.low), max(cred_pred_dev$std.resid.high)) +
  xlab("Average Predicted HbA1c (mmol/mol)") +
  ylab("Standardised Residuals")

# variable importance

if (class(try(
  
  variable_importance <- readRDS(paste0(output_path, "/Final_model/model_7/Treatment_Effects/variable_importance.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variable_importance <- bartMachine::investigate_var_importance(bart_model_effects)
  
  saveRDS(variable_importance, paste0(output_path, "/Final_model/model_7/Treatment_Effects/variable_importance.rds"))
  
}



plot_var_importance <- variable_importance$avg_var_props %>%
  t() %>%
  as.data.frame() %>%
  mutate(drugclass = drugclass_SGLT2 + drugclass_GLP1,
         drugline = drugline_2 + drugline_3 + drugline_4 + drugline_5,
         ncurrtx = ncurrtx_0 + ncurrtx_1 + ncurrtx_2 + ncurrtx_3,
         malesex = malesex_0 + malesex_1,
         Category = `Category_Active smoker` + `Category_Ex-smoker` + `Category_Non-smoker`) %>%
  select(-drugclass_SGLT2, -drugclass_GLP1, -drugline_2, -drugline_3, -drugline_4, -drugline_5, -ncurrtx_0, -ncurrtx_1, -ncurrtx_2, -ncurrtx_3, -malesex_0, -malesex_1, -`Category_Active smoker`, -`Category_Ex-smoker`, -`Category_Non-smoker`) %>%
  rename("Therapy" = "drugclass",
         "HbA1c" = "prehba1cmmol",
         # "Year Drug Start" = "yrdrugstart",
         "eGFR" = "egfr_ckdepi",
         "Outcome month" = "hba1cmonth",
         # "Systolic" = "presys",
         # "ALT" = "prealt",
         # "Bilirubin" = "prebil",
         # "AST" = "preast",
         # "HDL" = "prehdl",
         # "Age" = "agetx",
         "CVD score" = "score.excl.mi",
         # "Time to prescription" = "t2dmduration",
         # "BMI" = "prebmi",
         # "Albuminuria" = "prealb",
         # "Platelets" = "preplatelets",
         "Number of past therapies" = "drugline",
         "Number of Current therapies" = "ncurrtx",
         "Sex" = "malesex",
         "Smoking status" = "Category") %>%
  gather() %>%
  mutate(value = value * 100) %>%
  ggplot(aes(x = fct_rev(fct_reorder(key, value)), y = value)) +
  theme_classic() +
  geom_col() +
  ylab("Inclusion proportions (%)") +
  scale_y_continuous(limits = c(0, 30), expand = (c(0,0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate(x = 0, xend=0, y=0, yend=30, colour="black", lwd=0.75, geom="segment")


pdf(file = "Plots/7.3.model_treatment_effect.pdf")
plot_pred_obs
plot_resid
plot_var_importance
dev.off()




















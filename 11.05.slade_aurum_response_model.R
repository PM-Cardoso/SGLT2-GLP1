####################
## Description:
##  - In this file we:
##    - Fit a SparseBCF model to choose the variables to include as control and moderators.
##    - Fit a BCF model with the specific chosen variables used in control and in moderators.
##    - Validate model
####################


library(tidyverse)
library(SparseBCF)
require(bcf)

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
dir.create(paste0(output_path, "/response_model/trees"))

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


# Complete data
hba1c.train.complete <- hba1c.train %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()

if (class(try(
  
  sbcf_variable_selection <- readRDS(paste0(output_path, "/response_model/sbcf_variable_selection.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  sbcf_variable_selection <- SparseBCF(y = hba1c.train.complete$posthba1cfinal,
                                       z = hba1c.train.complete %>%
                                         select(drugclass) %>%
                                         mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>% 
                                         unlist(),
                                       x_control = hba1c.train.complete %>%
                                         select(
                                           -patid,
                                           -pated,
                                           -drugclass,
                                           -posthba1cfinal,
                                           -prop.score) %>%
                                         mutate_all(funs(as.numeric(.))) %>%
                                         as.matrix(), 
                                       pihat = 1-hba1c.train.complete$prop.score, 
                                       OOB = F, 
                                       sparse = T,
                                       update_interval = 250,
                                       ntree_control = 200,
                                       sd_control = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                       base_control = 0.95,
                                       power_control = 2,
                                       ntree_moderate = 200,
                                       sd_moderate = 2 * sd(hba1c.train.complete$posthba1cfinal),
                                       base_moderate = 0.95,
                                       power_moderate = 2,
                                       nburn = 200000,
                                       nsim = 50000
  )
  
  
  saveRDS(sbcf_variable_selection, paste0(output_path, "/response_model/sbcf_variable_selection.rds"))
  
}

### Selecting variables
variables <- hba1c.train.complete %>% select(-patid, -pated, -drugclass, -posthba1cfinal, -prop.score) %>% colnames()


variables_tau_original <- colMeans(sbcf_variable_selection$varprb_tau) %>% t() %>% as.data.frame()
colnames(variables_tau_original) <- variables

if (class(try(
  
  variables_tau <- readRDS(paste0(output_path, "/response_model/variables_tau.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  variables_tau <- colnames(variables_tau_original)[variables_tau_original > (1/(length(variables)))]
  
  saveRDS(variables_tau, file = paste0(output_path, "/response_model/variables_tau.rds"))
  
}


plot_variables_tau <- variables_tau_original %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > (1/(length(variables))), "Above", "Below")) %>%
  ggplot() +
  geom_hline(aes(yintercept = 1/length(variables)), colour = "red") +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = value, colour = colour)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ggtitle("tau model") +
  ylab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "green", "Below" = "red"))




variables_mu_original <- colMeans(sbcf_variable_selection$varprb_mu) %>% t() %>% as.data.frame() 
colnames(variables_mu_original) <- c(variables, "propensity score")


if (class(try(
  
  variables_mu <- readRDS(paste0(output_path, "/response_model/variables_mu.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  dataset.variables.mu <- variables_mu_original %>%
    select(-`propensity score`)
  
  variables_mu <- colnames(dataset.variables.mu)[dataset.variables.mu > (1/(length(variables)))]
  
  saveRDS(variables_mu, file = paste0(output_path, "/response_model/variables_mu.rds"))
  
}

plot_variables_mu <- variables_mu_original %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > (1/(length(variables))), "Above", "Below")) %>%
  ggplot() +
  geom_point(aes(x = forcats::fct_reorder(key, value), y = value, colour = colour)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  geom_hline(aes(yintercept = 1/length(variables)), colour = "red") +
  ggtitle("mu model") +
  ylab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "green", "Below" = "red"))



plot_sigma_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(sbcf_variable_selection$sigma), y = sbcf_variable_selection$sigma)) +
  ggtitle("Sigma") +
  ylab("Sigma") +
  theme(axis.title.x = element_blank())


plot_mu_scale_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(sbcf_variable_selection$mu_scale), y = sbcf_variable_selection$mu_scale)) +
  ggtitle("Mu") +
  ylab("Mu") +
  theme(axis.title.x = element_blank())


plot_tau_scale_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(sbcf_variable_selection$tau_scale), y = sbcf_variable_selection$tau_scale)) +
  ggtitle("Tau") +
  ylab("Tau") +
  theme(axis.title.x = element_blank())


plot_bscale0_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(sbcf_variable_selection$bscale0), y = sbcf_variable_selection$bscale0)) +
  ggtitle("B0") +
  ylab("B0") +
  theme(axis.title.x = element_blank())


plot_bscale1_sparsebcf <- ggplot() +
  geom_path(aes(x = 1:length(sbcf_variable_selection$bscale1), y = sbcf_variable_selection$bscale1)) +
  ggtitle("B1") +
  ylab("B1") +
  theme(axis.title.x = element_blank())


hba1c.train.complete.vs <- hba1c.train.complete %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score)


## Hold out cohort
hba1c.test <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test")


# collect propensity score values
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

hba1c.test <- hba1c.test %>%
  left_join(patient_prop_scores, by = c("patid", "pated"))


# Complete data
hba1c.test.complete.vs <- hba1c.test %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()


if (class(try(
  
  bcf_model <- readRDS(paste0(output_path, "/response_model/bcf_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  bcf_model = bcf::bcf(y = hba1c.train.complete.vs$posthba1cfinal,
                       z = hba1c.train.complete.vs %>%
                         select(drugclass) %>%
                         mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                         unlist(),
                       x_control = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_mu)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       x_moderate = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_tau)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       pihat = 1-hba1c.train.complete.vs$prop.score,
                       nburn = 200000,
                       nsim = 25000,
                       nthin = 1,
                       n_chains = 2,
                       # n_threads was max((RcppParallel::defaultNumThreads() - 2)/n_cores, 1) (this uses all of the server)
                       n_threads = 4,
                       update_interval = 500,
                       ntree_control = 200,
                       sd_control = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_control = 0.95,
                       power_control = 2,
                       ntree_moderate = 200,
                       sd_moderate = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_moderate = 0.95,
                       power_moderate = 2,
                       save_tree_directory = paste0(output_path, "/response_model/trees"))
  
  saveRDS(bcf_model, paste0(output_path, "/response_model/bcf_model.rds"))
  
}


if (class(try(
  
  predictions.hba1c.test <- readRDS(paste0(output_path, "/response_model/predictions.hba1c.test.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.hba1c.test <- predict(object = bcf_model,
                                    x_predict_control = hba1c.test.complete.vs %>%
                                      select(
                                        all_of(variables_mu)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    x_predict_moderate = hba1c.test.complete.vs %>%
                                      select(
                                        all_of(variables_tau)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    pi_pred = 1-hba1c.test.complete.vs$prop.score,
                                    z_pred = hba1c.test.complete.vs %>%
                                      select(drugclass) %>%
                                      mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                      unlist(),
                                    save_tree_directory = paste0(output_path, "/response_model/trees"))
  
  
  saveRDS(predictions.hba1c.test, paste0(output_path, "/response_model/predictions.hba1c.test.rds"))
  
}


#### Model parameters


plot_sigma_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$sigma), y = bcf_model$sigma)) +
  ggtitle("Sigma") +
  ylab("Sigma") +
  theme(axis.title.x = element_blank())


plot_mu_scale_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$mu_scale), y = bcf_model$mu_scale)) +
  ggtitle("Mu") +
  ylab("Mu") +
  theme(axis.title.x = element_blank())


plot_tau_scale_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$tau_scale), y = bcf_model$tau_scale)) +
  ggtitle("Tau") +
  ylab("Tau") +
  theme(axis.title.x = element_blank())


plot_bscale0_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$b0), y = bcf_model$b0)) +
  ggtitle("B0") +
  ylab("B0") +
  theme(axis.title.x = element_blank())


plot_bscale1_bcf <- ggplot() +
  geom_path(aes(x = 1:length(bcf_model$b1), y = bcf_model$b1)) +
  ggtitle("B1") +
  ylab("B1") +
  theme(axis.title.x = element_blank())

#### Residuals

cred_pred_dev <- calc_resid(hba1c.train.complete.vs, bcf_model$mu, "posthba1cfinal")

cred_pred_val <- calc_resid(hba1c.test.complete.vs, predictions.hba1c.test$mu, "posthba1cfinal")

plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Standardised Residuals of BCF model")


#### Effects

data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
  as.data.frame()

plot_effect_1 <- hist_plot(data_dev, "", -15, 25)

data_val <- cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  as.data.frame()

plot_effect_2 <- hist_plot(data_val, "", -15, 25)

####### Strata effect
## male sex
sex_male_dev <- hba1c.train.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(bcf_model$tau)) %>%
  filter(sex == "Male") %>%
  select(-sex)

plot_effect_1_male <- hist_plot(sex_male_dev, "", -15, 25)

## female sex
sex_female_dev <- hba1c.train.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(bcf_model$tau)) %>%
  filter(sex == "Female") %>%
  select(-sex)

plot_effect_1_female <- hist_plot(sex_female_dev, "", -15, 25)


## male sex
sex_male_val <- hba1c.test.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  filter(sex == "Male") %>%
  select(-sex)

plot_effect_2_male <- hist_plot(sex_male_val, "", -15, 25)

## female sex
sex_female_val <- hba1c.test.complete.vs %>%
  select(sex) %>%
  cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  filter(sex == "Female") %>%
  select(-sex)

plot_effect_2_female <- hist_plot(sex_female_val, "", -15, 25)


####### Validation
predicted_observed_dev <- hba1c.train.complete.vs %>%
  cbind(hba1c_diff = colMeans(bcf_model$tau)) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


ATE_matching_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_dev, "posthba1cfinal", hba1c.train.complete.vs$prop.score)

plot_ATE_dev_prop_score <- ATE_plot(ATE_matching_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "posthba1cfinal", hba1c.train.complete.vs$prop.score)

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")



predicted_observed_val <- hba1c.test.complete.vs %>%
  cbind(hba1c_diff = colMeans(predictions.hba1c.test$tau)) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


ATE_matching_validation_val <- calc_ATE_validation_prop_matching(predicted_observed_val, "posthba1cfinal", hba1c.test.complete.vs$prop.score)

plot_ATE_val_prop_score <- ATE_plot(ATE_matching_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1cfinal", hba1c.test.complete.vs$prop.score)

plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")




pdf(width = 10, file = "Plots/11.05.slade_aurum_bcf.pdf")
# convergence of the sparsebcf model analysis
patchwork::wrap_plots(list(plot_sigma_sparsebcf, plot_mu_scale_sparsebcf, plot_tau_scale_sparsebcf, plot_bscale0_sparsebcf, plot_bscale1_sparsebcf), ncol = 1) +
  patchwork::plot_annotation(title = "Parameters for SparseBCF")

# variables for treatment effects
plot_variables_tau
# variables for response
plot_variables_mu

# convergence of the bcf model analysis
patchwork::wrap_plots(list(plot_sigma_bcf, plot_mu_scale_bcf, plot_tau_scale_bcf, plot_bscale0_bcf, plot_bscale1_bcf), ncol = 1) +
  patchwork::plot_annotation(title = "Parameters for BCF")

# Standardised residuals for the bcf model
plot_residuals

patchwork::wrap_plots(list(plot_effect_1, plot_effect_2)) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Treatment effect heterogeneity (A-Developement)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

patchwork::wrap_plots(list(plot_effect_1_male, plot_effect_1_female), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Developement: Treatment effect heterogeneity (A-Male)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot


patchwork::wrap_plots(list(plot_effect_2_male, plot_effect_2_female), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Validation: Treatment effect heterogeneity (A-Male)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

patchwork::wrap_plots(list(plot_ATE_dev_prop_score, plot_ATE_val_prop_score)) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Effects validation propensity score matching (A-Development)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5))) # center title of full plot


patchwork::wrap_plots(list(plot_ATE_dev_prop_score_weighting, plot_ATE_val_prop_score_weighting)) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Effects validation propensity score matching (A-Development)", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5))) # center title of full plot


dev.off()

#:-------------------------------------------------------------------------------------
#

library(rpart)
library(rattle)
library(rpart.plot)

rpart.dataset <- predicted_observed_dev

fit <- rpart(hba1c_diff ~ agetx + sex + drugline + ncurrtx + hba1cmonth + prehba1c + preegfr + preihd + preneuropathy + preretinopathy + preaf, data = rpart.dataset)

rpart.dataset.strata <- predicted_observed_dev %>%
filter(ncurrtx == "2" | ncurrtx == "1") %>%
filter(drugline == "2") %>%
filter(hba1cmonth > 6)

fit2 <- rpart(hba1c_diff ~ agetx + sex + drugline + ncurrtx + hba1cmonth + prehba1c + preegfr + preihd + preneuropathy + preretinopathy + preaf, data = rpart.dataset.strata)


pdf(width = 20, height = 8, file = "Plots/11.05.effect_decision_tree.pdf")

prp(fit, pal.thresh = 0, box.palette="BuGn", extra = "auto", main = "Decision tree for treatment effects using development cohort")

prp(fit2, pal.thresh = 0, box.palette="BuGn", extra = "auto", main = "Strata: drugline = 2, ncurrtx == 1/2, hba1cmonth > 6")

dev.off()



#:-------------------------------------------------------------------------------------
# Prediction of treatment effects for all patients

# patients which already have a calculated effect
patient_effects <- hba1c.train.complete.vs %>%
  select(patid, pated) %>%
  cbind(effects = colMeans(bcf_model$tau)) %>%
  rbind(
    hba1c.test.complete.vs %>%
      select(patid, pated) %>%
      cbind(effects = colMeans(predictions.hba1c.test$tau))
  )


# patients not yet calculated
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

interim.dataset <- full.cohort[!(full.cohort$pated %in% patient_effects$pated),] %>%
  # left_join propensity scores%>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  # select variables to make prediction
  select(patid, pated, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # this results in a lot of missing hba1cmonth, should I set these to 12 months?
  drop_na()


if (class(try(
  
  patient_effects <- readRDS(paste0(output_path, "/response_model/patient_effects.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions.interim <- predict(object = bcf_model,
                                    x_predict_control = interim.dataset %>%
                                      select(
                                        all_of(variables_mu)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    x_predict_moderate = interim.dataset %>%
                                      select(
                                        all_of(variables_tau)
                                      ) %>%
                                      mutate_all(funs(as.numeric(.))) %>%
                                      as.matrix(),
                                    pi_pred = 1-interim.dataset$prop.score,
                                    z_pred = interim.dataset %>%
                                      select(drugclass) %>%
                                      mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                      unlist(),
                                    save_tree_directory = paste0(output_path, "/response_model/trees"))
  
  
  
  patient_effects <- patient_effects %>%
    rbind(
      interim.dataset %>%
        select(patid, pated) %>%
        cbind(effects = colMeans(predictions.interim$tau))
    )
  
  
  saveRDS(patient_effects, paste0(output_path, "/response_model/patient_effects.rds"))
  
}











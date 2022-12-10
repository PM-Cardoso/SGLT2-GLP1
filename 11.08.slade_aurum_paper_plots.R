####################
## Description:
##  - In this file we create the plots used in the paper
####################

library(tidyverse)
library(patchwork)

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")


#:------------------------------------------------------------------------------
# Validation of BCF model (training cohort, trying background colours)

ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")

plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_dev[["effects"]]$N), ")"))

plot_1 <- patchwork::wrap_plots(list(plot_ATE_adjust_validation_dev), ncol = 1) +
  patchwork::plot_annotation(title = "Validation of treatment effects for training cohort") # title of full plot

pdf(width = 7, height = 8, "Plots/Plot_1.pdf")
plot_1
dev.off()


#:------------------------------------------------------------------------------
# Validation of BCF model (histogram + adjusted validation)

# Validation plot dev
ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")

plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_dev[["effects"]]$N), ")"))

# histogram dev
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
  as.data.frame()

plot_effect_dev <- hist_plot(data_dev, "Treatment effects heterogeneity", -15, 20)

# Validation plot val
ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")

plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_val[["effects"]]$N), ")"))

# histogram val
predictions.hba1c.test <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/predictions.hba1c.test.rds")

data_val <- cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
  as.data.frame()

plot_effect_val <- hist_plot(data_val, "Treatment effects heterogeneity", -15, 20)


plot_2.1 <- (plot_effect_dev + plot_ATE_adjust_validation_dev) / (plot_effect_val + plot_ATE_adjust_validation_val)

plot_2.1[[1]] <- plot_2.1[[1]] + plot_layout(tag_level = 'new')
plot_2.1[[2]] <- plot_2.1[[2]] + plot_layout(tag_level = 'new')

plot_2 <- plot_2.1 + 
  patchwork::plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c('A', '1'),
                  title = "Predicted treatment effects in the training and hold-out cohorts",
                  theme = theme(legend.position = "bottom")) # title of full plot

pdf(width = 10, height = 9, "Plots/Plot_2.pdf")
plot_2
dev.off()

#:------------------------------------------------------------------------------
# Other validation plots
# PSM 1:1 posthba1cfinal ~ drugclass
ATE_matching_1_1_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_dev.rds")

plot_ATE_matching_1_1_validation_dev <- ATE_plot(ATE_matching_1_1_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Propensity score matching (n=",sum(ATE_matching_1_1_validation_dev[["effects"]]$n_drug1)*2, ")"))

# PSM 1:1 posthba1cfinal ~ drugclass + adjustment
ATE_matching_1_1_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_dev.rds")

plot_ATE_matching_1_1_adjust_validation_dev <- ATE_plot(ATE_matching_1_1_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Propensity score matching + adjusted (n=",sum(ATE_matching_1_1_adjust_validation_dev[["effects"]]$n_drug1)*2, ")"))

# PSM 1:1 posthba1cfinal ~ drugclass
ATE_matching_1_1_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_val.rds")

plot_ATE_matching_1_1_validation_val <- ATE_plot(ATE_matching_1_1_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Propensity score matching (n=",sum(ATE_matching_1_1_validation_val[["effects"]]$n_drug1)*2, ")"))

# PSM 1:1 posthba1cfinal ~ drugclass + adjustment
ATE_matching_1_1_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_val.rds")

plot_ATE_matching_1_1_adjust_validation_val <- ATE_plot(ATE_matching_1_1_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Propensity score matching + adjusted (n=",sum(ATE_matching_1_1_adjust_validation_val[["effects"]]$n_drug1)*2, ")"))


plot_3.1 <- (plot_ATE_matching_1_1_validation_dev + plot_ATE_matching_1_1_adjust_validation_dev) / (plot_ATE_matching_1_1_validation_val + plot_ATE_matching_1_1_adjust_validation_val)

plot_3.1[[1]] <- plot_3.1[[1]] + plot_layout(tag_level = 'new')
plot_3.1[[2]] <- plot_3.1[[2]] + plot_layout(tag_level = 'new')

plot_3 <- plot_3.1 + 
  patchwork::plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c('A', '1'),
                  title = "Predicted treatment effects in the training and hold-out cohorts",
                  theme = theme(legend.position = "bottom")) # title of full plot

pdf(width = 10, height = 9, "Plots/Plot_3.pdf")
plot_3
dev.off()




#:------------------------------------------------------------------------------
# Residuals 
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
predictions.hba1c.test <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/predictions.hba1c.test.rds")

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

# collect propensity score values
patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

# R2/RMSE
assessment_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/assessment_dev.rds")
assessment_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/assessment_val.rds")


hba1c.train.complete.vs <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score)

hba1c.test.complete.vs <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  # selected variables from SparseBCF
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()


cred_pred_dev <- calc_resid(hba1c.train.complete.vs, bcf_model$mu, "posthba1cfinal")

cred_pred_val <- calc_resid(hba1c.test.complete.vs, predictions.hba1c.test$mu, "posthba1cfinal")

plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Standardised Residuals of BCF model")

# plot histogram of standardised residuals
plot_hist_dev <- cred_pred_dev %>%
  ggplot() +
  geom_density(aes(x = std.resid)) +
  xlab("Standardised Mean Residuals") +
  xlim(min(c(cred_pred_dev$std.resid, cred_pred_val$std.resid)), max(c(cred_pred_dev$std.resid, cred_pred_val$std.resid))) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank())

plot_hist_val <- cred_pred_val %>%
  ggplot() +
  geom_density(aes(x = std.resid)) +
  xlab("Standardised Mean Residuals") +
  xlim(min(c(cred_pred_dev$std.resid, cred_pred_val$std.resid)), max(c(cred_pred_dev$std.resid, cred_pred_val$std.resid))) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank())


plot_4.1 <- (plot_residuals[[1]] +
               ggtitle(paste0("R2: ", round(assessment_dev$r2[2], 2), " (", round(assessment_dev$r2[1], 2), "-", round(assessment_dev$r2[3], 2),") RMSE: ", round(assessment_dev$RMSE[2], 2), " (", round(assessment_dev$RMSE[1], 2),"-", round(assessment_dev$RMSE[3], 2),")")) + plot_hist_dev) / (plot_residuals[[2]] +
                                                                                                                                                                                                                                                                                                       ggtitle(paste0("R2: ", round(assessment_val$r2[2], 2), " (", round(assessment_val$r2[1], 2), "-", round(assessment_val$r2[3], 2),") RMSE: ", round(assessment_val$RMSE[2], 2), " (", round(assessment_val$RMSE[1], 2),"-", round(assessment_val$RMSE[3], 2),")")) + plot_hist_val)

plot_4.1[[1]] <- plot_4.1[[1]] + plot_layout(tag_level = 'new')
plot_4.1[[2]] <- plot_4.1[[2]] + plot_layout(tag_level = 'new')

plot_4 <- plot_4.1 + 
  plot_annotation(tag_levels = c('A', '1'),
                  title = "Standardised residuals in the training and hold-out cohorts",
                  theme = theme(legend.position = "bottom")) # title of full plot


pdf(width = 10, height = 9, "Plots/Plot_4.pdf")
plot_4
dev.off()


#:-----------------------------------------------------------------------------
# Variable importance

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

# load BCF model
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

# load development dataset
hba1c.train.cleaned_up <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau))) %>%
  cbind(effects = colMeans(bcf_model$tau))


m1 <- rms::ols(effects ~ rms::rcs(agetx, 3) + sex + ncurrtx + rms::rcs(prehba1c, 3) + rms::rcs(prebmi, 3) + rms::rcs(preegfr, 3) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = hba1c.train.cleaned_up, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_5 <- as.data.frame(values) %>%
  cbind(variable = c("Number of other current glucose-lowering drugs", "Sex", "eGFR", "Current age", "BMI", "HbA1c", "Retinopathy", "Peripheral arterial disease", "Neuropathy", "Ischaemic heart disease", "Heart failure")) %>%
  mutate(variable = factor(variable),
         values = values * 100) %>%
  ggplot(aes(y = forcats::fct_reorder(variable, values), x = values)) +
  geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(variable, values)), linetype = "dashed") +
  geom_point(size = 2, colour = "black") +
  ggtitle("Relative importance for treatment effect heterogeneity") +
  xlab("Relative Importance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45),
        axis.title.y = element_blank())


pdf("Plots/Plot_5.pdf")
plot_5
dev.off()



#:------------------------------------------------------------------------------
# Heterogeneity in sex

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

# load BCF model
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

# load development dataset
hba1c.train.heterogeneity <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  # drop the variables with the most missingness
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame() %>%
  select(sex) %>%
  cbind(mean = colMeans(bcf_model$tau))

plot_effects_female <- hist_plot(hba1c.train.heterogeneity %>% filter(sex == "Female"), "Females", -15, 18) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal")

plot_effects_male <- hist_plot(hba1c.train.heterogeneity %>% filter(sex == "Male"), "Males", -15, 18) +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal")

plot_6 <- patchwork::wrap_plots(list(
  plot_effects_female,
  plot_effects_male
), ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Heterogeneity of treatment effects",  # title of full plot
                             theme = theme(legend.direction = "horizontal", 
                                           legend.position = "bottom",
                                           legend.box = "horizontal"))


pdf(width = 8, height = 5, "Plots/Plot_6.pdf")
plot_6
dev.off()














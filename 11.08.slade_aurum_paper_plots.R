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


m1 <- rms::ols(effects ~ rms::rcs(agetx, 3) + rms::rcs(hba1cmonth, 3) + sex + ncurrtx + rms::rcs(prehba1c, 3) + rms::rcs(prebmi, 3) + rms::rcs(preegfr, 3) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = hba1c.train.cleaned_up, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_5 <- as.data.frame(values) %>%
  cbind(variable = c("Number of other current glucose-lowering drugs", "Sex", "eGFR", "Current age", "BMI", "HbA1c", "Retinopathy", "Peripheral arterial disease", "Neuropathy", "Ischaemic heart disease", "Heart failure", "Month of outcome")) %>%
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




#:------------------------------------------------------------------------------
# Validation of BCF model (population with/without CVD/HF/CKD)

patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

# Population without CVD/HF/CKD
no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(patient_effects, by = c("patid", "pated"))

predicted_no_co.dataset <- no_co.dataset %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10)) %>%
  drop_na(posthba1cfinal, hba1c_diff)


# Population with CVD/HF/CKD
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

co.dataset <- full.cohort[!(full.cohort$pated %in% no_co.dataset$pated),] %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(patient_effects, by = c("patid", "pated"))

predicted_co.dataset <- co.dataset %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10)) %>%
  drop_na(posthba1cfinal, hba1c_diff)


##:-- No comorbidities
# PSM 1:1
ATE_matching_1_1_no_co.dataset <- calc_ATE_validation_prop_matching(predicted_no_co.dataset, "posthba1cfinal", predicted_no_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

plot_ATE_matching_1_1_no_co.dataset <- ATE_plot(ATE_matching_1_1_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PSM 1:1 + adjusted
ATE_matching_1_1_adjusted_no_co.dataset <- calc_ATE_validation_prop_matching(predicted_no_co.dataset, "posthba1cfinal", predicted_no_co.dataset$prop.score, order = "largest", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)

plot_ATE_matching_1_1_adjusted_no_co.dataset <- ATE_plot(ATE_matching_1_1_adjusted_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# Adjusted
ATE_adjusted_no_co.dataset <- calc_ATE_validation_adjust(predicted_no_co.dataset, "posthba1cfinal", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)

plot_ATE_adjusted_no_co.dataset <- ATE_plot(ATE_adjusted_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

##:-- Comorbidities
# PSM 1:1
ATE_matching_1_1_co.dataset <- calc_ATE_validation_prop_matching(predicted_co.dataset, "posthba1cfinal", predicted_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

plot_ATE_matching_1_1_co.dataset <- ATE_plot(ATE_matching_1_1_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# PSM 1:1 + adjusted
ATE_matching_1_1_adjusted_co.dataset <- calc_ATE_validation_prop_matching(predicted_co.dataset, "posthba1cfinal", predicted_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)

plot_ATE_matching_1_1_adjusted_co.dataset <- ATE_plot(ATE_matching_1_1_adjusted_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")

# Adjusted
ATE_adjusted_co.dataset <- calc_ATE_validation_adjust(predicted_co.dataset, "posthba1cfinal", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)

plot_ATE_adjusted_co.dataset <- ATE_plot(ATE_adjusted_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")


# No comorbidities plots
plot_7.1 <- patchwork::wrap_plots(list(plot_ATE_matching_1_1_no_co.dataset + ggtitle("Match 1:1"),
                                   plot_ATE_matching_1_1_adjusted_no_co.dataset + ggtitle("Match 1:1 adjusted"),
                                   plot_ATE_adjusted_no_co.dataset + ggtitle("Adjusted")), ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = paste0("No comorbidities cohort: Treatment effect validation (n=", predicted_no_co.dataset%>%nrow(),")"), # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

# Comorbidities plots
plot_7.2 <- patchwork::wrap_plots(list(plot_ATE_matching_1_1_co.dataset + ggtitle("Match 1:1"),
                                       plot_ATE_matching_1_1_adjusted_co.dataset + ggtitle("Match 1:1 adjusted"),
                                       plot_ATE_adjusted_co.dataset + ggtitle("Adjusted")), ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A",
                             title = paste0("Comorbidities cohort: Treatment effect validation (n=", predicted_co.dataset%>%nrow(),")"), # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

pdf(width = 12, height = 5, "Plots/Plot_7.pdf")
plot_7.1
plot_7.2
dev.off()


#:------------------------------------------------------------------------------
# HbA1c grouping (population with/without CVD/HF/CKD)

require(forestplot)


patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")

patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")

# load variables used in the BCF model
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

interval_breaks <- c(-5, -3, 0, 3, 5)

# Population without CVD/HF/CKD
no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(patient_effects, by = c("patid", "pated"))

group.hba1c.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals, posthba1cfinal) %>%
  rename("hba1c_diff" = "effects")


# Population with CVD/HF/CKD
full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

co.dataset <- full.cohort[!(full.cohort$pated %in% no_co.dataset$pated),] %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(patient_effects, by = c("patid", "pated"))

group.hba1c.co.dataset <- group_values(data = co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals, posthba1cfinal) %>%
  rename("hba1c_diff" = "effects")


##:-- No comorbidities
# PSM 1:1
ATE_matching_1_1_no_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.no_co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

# PSM 1:1 + adjusted
ATE_matching_1_1_adjusted_no_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.no_co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)

# Adjusted
ATE_adjusted_no_co.dataset <- calc_ATE_validation_adjust(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)

##:-- Comorbidities
# PSM 1:1
ATE_matching_1_1_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)))

# PSM 1:1 + adjusted
ATE_matching_1_1_adjusted_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)

# Adjusted
ATE_adjusted_co.dataset <- calc_ATE_validation_adjust(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)


hba1c_overall_axis_min <- plyr::round_any(floor(min(c(ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_matching_1_1_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_overall_axis_max <- plyr::round_any(ceiling(max(c(ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_matching_1_1_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)


# Propensity score matching
plot_psm_1_1_no_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("No comorbidities population (n=", sum(ATE_matching_1_1_no_co.dataset[["effects"]]$n_drug1)*2, ")"))

plot_psm_1_1_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Comorbidities population (n=", sum(ATE_matching_1_1_co.dataset[["effects"]]$n_drug1)*2, ")"))


# Propensity score matching + adjusted
plot_psm_1_1_adjusted_no_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("No comorbidities population (n=", sum(ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]$n_drug1)*2, ")"))

plot_psm_1_1_adjusted_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Comorbidities population (n=", sum(ATE_matching_1_1_adjusted_co.dataset[["effects"]]$n_drug1)*2, ")"))


# Adjusted
plot_adjusted_no_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(1)%>%unlist(), ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(2)%>%unlist(), ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(3)%>%unlist(), ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(4)%>%unlist(), ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(5)%>%unlist(), ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(6)%>%unlist(), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("No comorbidities population (n=", sum(ATE_adjusted_no_co.dataset[["effects"]]$N), ")"))

plot_adjusted_co_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(1)%>%unlist(), ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(2)%>%unlist(), ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(3)%>%unlist(), ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(4)%>%unlist(), ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(5)%>%unlist(), ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(6)%>%unlist(), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Comorbidities population (n=", sum(ATE_adjusted_co.dataset[["effects"]]$N), ")"))



pdf(width = 15, height = 12, "Plots/Plot_8.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HbA1c treatment effect", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_psm_1_1_no_co_hba1c
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_psm_1_1_co_hba1c
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_psm_1_1_adjusted_no_co_hba1c
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_psm_1_1_adjusted_co_hba1c
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_adjusted_no_co_hba1c
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_adjusted_co_hba1c
upViewport()

dev.off()



#:------------------------------------------------------------------------------
# Additional outcomes: adjusted model
require(forestplot)

output_path <- "Samples/SGLT2-GLP1/Aurum"

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Read in treatment effects
treatment_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

model_variables <- unique(c(variables_mu, variables_tau))[which(unique(c(variables_mu, variables_tau)) != "sex")]


# HbA1c grouping
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.hba1c.dataset <- group_values(data = hba1c.train,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")


ATE_adjust_hba1c <- calc_ATE_validation_adjust(group.hba1c.dataset, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_full <- calc_ATE_validation_adjust(group.hba1c.dataset %>% mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)


hba1c_strata_axis_min <- plyr::round_any(floor(min(c(ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)


plot_adjust_overall_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[1])%>%nrow(), ")"), 
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[3])%>%nrow(), ")"), 
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[4])%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[5])%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[6])%>%nrow(), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "HbA1c",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", group.hba1c.dataset%>%nrow(), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



# Weight change
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the weight adjusted model
predictions_weight_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_overall.rds"))

# predictions for the weight adjusted model full
predictions_weight_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_full.rds"))


weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                        predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


plot_weight_adjusted_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", group.weight.dataset%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


# eGFR change

egfr.dataset <- set_up_data_sglt2_glp1(dataset.type = "egfr.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(egfr.change = postegfr - preegfr)


group.egfr.dataset <- group_values(data = egfr.dataset,
                                   variable = "effects",
                                   breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the egfr adjusted model
predictions_egfr_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_overall.rds"))

# predictions for the egfr adjusted model sex strata
predictions_egfr_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_full.rds"))

egfr_strata_axis_min <- plyr::round_any(floor(min(c(predictions_egfr_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                    predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

egfr_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_egfr_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                      predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


plot_egfr_adjusted_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_egfr_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(intervals==levels(group.egfr.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "eGFR change",
             clip = c(egfr_strata_axis_min, egfr_strata_axis_max),
             xticks = seq(egfr_strata_axis_min, egfr_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Overall population (n=", group.egfr.dataset%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



# Discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)


# predictions for the discontinuation adjusted model
predictions_discontinuation_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds"))

# predictions for the discontinuation model sex strata
predictions_discontinuation_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds"))

discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                  predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                 predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)


plot_discontinuation_adjusted_overall <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_adjusted_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_adjusted_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Discontinuation",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", group.discontinuation.dataset%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


# No comorbidities

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

## CVD outcomes

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds"))

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_cvd_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds"))

plot_no_co_cvd_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "CVD outcomes",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


## HF outcomes

# predictions for the HF outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_overall.rds"))

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_no_co_hf_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_full.rds"))

plot_no_co_hf_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "HF outcomes",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



pdf(width = 20, height = 10, "Plots/Plot_9_a.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 3, heights = unit(c(5, 5), "null"))))
# 1 2 3
# 4 5 6

# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_adjust_overall_hba1c
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 2))
plot_weight_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 3))
plot_egfr_adjusted_overall
upViewport()

# forth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_discontinuation_adjusted_overall
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_no_co_cvd_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 3))
plot_no_co_hf_adjusted_overall
upViewport()

dev.off()

pdf(width = 20, height = 10, "Plots/Plot_9_b.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 3, heights = unit(c(5, 5), "null"))))
# 1 2 3
# 4 5 6

# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_adjust_overall_hba1c
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 2))
plot_no_co_cvd_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 3))
plot_no_co_hf_adjusted_overall
upViewport()

# forth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_discontinuation_adjusted_overall
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_weight_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 3))
plot_egfr_adjusted_overall
upViewport()

dev.off()
















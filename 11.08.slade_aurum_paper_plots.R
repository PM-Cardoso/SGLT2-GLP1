####################
## Description:
##  - In this file we create the plots used in the paper
####################

library(tidyverse)
library(patchwork)
library(scales)

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

# 
# #:------------------------------------------------------------------------------
# # Validation of BCF model (training cohort, trying background colours)
# 
# ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")
# 
# plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = TRUE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_dev[["effects"]]$N), ")"))
# 
# plot_1 <- patchwork::wrap_plots(list(plot_ATE_adjust_validation_dev), ncol = 1) +
#   patchwork::plot_annotation(title = "Validation of treatment effects for training cohort") # title of full plot
# 
# pdf(width = 7, height = 8, "Plots/11.08.plot_1.pdf")
# plot_1
# dev.off()
# 
# 
# #:------------------------------------------------------------------------------
# # Validation of BCF model (histogram + adjusted validation)
# 
# # Validation plot dev
# ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")
# 
# plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_dev[["effects"]]$N), ")"))
# 
# # histogram dev
# bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
# 
# data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
#   as.data.frame()
# 
# plot_effect_dev <- hist_plot(data_dev, "Treatment effects heterogeneity", -15, 20)
# 
# # Validation plot val
# ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")
# 
# plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Adjusted (n=",sum(ATE_adjust_validation_val[["effects"]]$N), ")"))
# 
# # histogram val
# predictions.hba1c.test <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/predictions.hba1c.test.rds")
# 
# data_val <- cbind(mean = colMeans(predictions.hba1c.test$tau)) %>%
#   as.data.frame()
# 
# plot_effect_val <- hist_plot(data_val, "Treatment effects heterogeneity", -15, 20)
# 
# 
# plot_2.1 <- (plot_effect_dev + plot_ATE_adjust_validation_dev) / (plot_effect_val + plot_ATE_adjust_validation_val)
# 
# plot_2.1[[1]] <- plot_2.1[[1]] + plot_layout(tag_level = 'new')
# plot_2.1[[2]] <- plot_2.1[[2]] + plot_layout(tag_level = 'new')
# 
# plot_2 <- plot_2.1 + 
#   patchwork::plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = c('A', '1'),
#                   title = "Predicted treatment effects in the training and hold-out cohorts",
#                   theme = theme(legend.position = "bottom")) # title of full plot
# 
# pdf(width = 10, height = 9, "Plots/11.08.plot_2.pdf")
# plot_2
# dev.off()
# 
# #:------------------------------------------------------------------------------
# # Other validation plots
# 
# # PSM 1:1 posthba1cfinal ~ drugclass
# ATE_matching_1_1_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_dev.rds")
# 
# plot_ATE_matching_1_1_validation_dev <- ATE_plot(ATE_matching_1_1_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Propensity score matching (n=",sum(ATE_matching_1_1_validation_dev[["effects"]]$n_drug1)*2, ")"))
# 
# # PSM 1:1 posthba1cfinal ~ drugclass + adjustment
# ATE_matching_1_1_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_dev.rds")
# 
# plot_ATE_matching_1_1_adjust_validation_dev <- ATE_plot(ATE_matching_1_1_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Propensity score matching + adjusted (n=",sum(ATE_matching_1_1_adjust_validation_dev[["effects"]]$n_drug1)*2, ")"))
# 
# # PSM 1:1 posthba1cfinal ~ drugclass
# ATE_matching_1_1_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_validation_val.rds")
# 
# plot_ATE_matching_1_1_validation_val <- ATE_plot(ATE_matching_1_1_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Propensity score matching (n=",sum(ATE_matching_1_1_validation_val[["effects"]]$n_drug1)*2, ")"))
# 
# # PSM 1:1 posthba1cfinal ~ drugclass + adjustment
# ATE_matching_1_1_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_matching_1_1_adjust_validation_val.rds")
# 
# plot_ATE_matching_1_1_adjust_validation_val <- ATE_plot(ATE_matching_1_1_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Propensity score matching + adjusted (n=",sum(ATE_matching_1_1_adjust_validation_val[["effects"]]$n_drug1)*2, ")"))
# 
# 
# plot_3.1 <- (plot_ATE_matching_1_1_validation_dev + plot_ATE_matching_1_1_adjust_validation_dev) / (plot_ATE_matching_1_1_validation_val + plot_ATE_matching_1_1_adjust_validation_val)
# 
# plot_3.1[[1]] <- plot_3.1[[1]] + plot_layout(tag_level = 'new')
# plot_3.1[[2]] <- plot_3.1[[2]] + plot_layout(tag_level = 'new')
# 
# plot_3 <- plot_3.1 + 
#   patchwork::plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = c('A', '1'),
#                   title = "Predicted treatment effects in the training and hold-out cohorts",
#                   theme = theme(legend.position = "bottom")) # title of full plot
# 
# pdf(width = 10, height = 9, "Plots/11.08.plot_3.pdf")
# plot_3
# dev.off()
# 
#:------------------------------------------------------------------------------
# # Residuals
# 
# bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
# predictions.hba1c.test <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/predictions.hba1c.test.rds")
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# # collect propensity score values
# patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")
# 
# # R2/RMSE
# assessment_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/assessment_dev.rds")
# assessment_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/assessment_val.rds")
# 
# hba1c.train.complete.vs <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   # drop the variables with the most missingness
#   select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
#   # only complete cases
#   drop_na() %>%
#   as.data.frame() %>%
#   select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score)
# 
# hba1c.test.complete.vs <- set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   # selected variables from SparseBCF
#   select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score) %>%
#   # only complete cases
#   drop_na() %>%
#   as.data.frame()
# 
# 
# # ## calculate slope and intercept
# # slope_vector_dev <- c()
# # slope_vector_val <- c()
# # intercept_vector_dev <- c()
# # intercept_vector_val <- c()
# # for (i in 1:nrow(bcf_model$mu)) {
# #   #:------------- Development dataset
# #   # iterate through 
# #   pred <- bcf_model$mu[i, ]
# #   orig <- hba1c.train.complete.vs$posthba1cfinal
# #   # model
# #   mod <- rms::ols(orig ~ pred, x = TRUE, y = TRUE)
# #   # validate
# #   val <- validate(mod, B=10)
# #   slope_vector_dev <- c(slope_vector_dev, val["Slope", 1])
# #   intercept_vector_dev <- c(intercept_vector_dev, val["Intercept", 1])
# #   
# #   #:------------ Validation dataset
# #   # iterate through
# #   pred <- predictions.hba1c.test$mu[i,]
# #   orig <- hba1c.test.complete.vs$posthba1cfinal
# #   # model
# #   mod <- rms::ols(orig ~ pred, x = TRUE, y = TRUE)
# #   # validate
# #   val <- validate(mod, B=10)
# #   slope_vector_val <- c(slope_vector_val, val["Slope", 1])
# #   intercept_vector_val <- c(intercept_vector_val, val["Intercept", 1])
# # }
# 
# 
# cred_pred_dev <- calc_resid(hba1c.train.complete.vs, bcf_model$mu, "posthba1cfinal")
# 
# cred_pred_val <- calc_resid(hba1c.test.complete.vs, predictions.hba1c.test$mu, "posthba1cfinal")
# 
# plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Standardised Residuals of BCF model")
# 
# # plot histogram of standardised residuals
# plot_hist_dev <- cred_pred_dev %>%
#   ggplot() +
#   geom_density(aes(x = std.resid)) +
#   xlab("Standardised Mean Residuals") +
#   xlim(min(c(cred_pred_dev$std.resid, cred_pred_val$std.resid)), max(c(cred_pred_dev$std.resid, cred_pred_val$std.resid))) +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         panel.grid = element_blank())
# 
# plot_hist_val <- cred_pred_val %>%
#   ggplot() +
#   geom_density(aes(x = std.resid)) +
#   xlab("Standardised Mean Residuals") +
#   xlim(min(c(cred_pred_dev$std.resid, cred_pred_val$std.resid)), max(c(cred_pred_dev$std.resid, cred_pred_val$std.resid))) +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         panel.grid = element_blank())
# 
# 
# plot_4.1 <- (plot_residuals[[1]] +
#                ggtitle(paste0("R2: ", round(assessment_dev$r2[2], 2), " (", round(assessment_dev$r2[1], 2), "-", round(assessment_dev$r2[3], 2),") RMSE: ", round(assessment_dev$RMSE[2], 2), " (", round(assessment_dev$RMSE[1], 2),"-", round(assessment_dev$RMSE[3], 2),")")) + plot_hist_dev) / (plot_residuals[[2]] +
#                                                                                                                                                                                                                                                                                                        ggtitle(paste0("R2: ", round(assessment_val$r2[2], 2), " (", round(assessment_val$r2[1], 2), "-", round(assessment_val$r2[3], 2),") RMSE: ", round(assessment_val$RMSE[2], 2), " (", round(assessment_val$RMSE[1], 2),"-", round(assessment_val$RMSE[3], 2),")")) + plot_hist_val)
# 
# plot_4.1[[1]] <- plot_4.1[[1]] + plot_layout(tag_level = 'new')
# plot_4.1[[2]] <- plot_4.1[[2]] + plot_layout(tag_level = 'new')
# 
# plot_4 <- plot_4.1 +
#   plot_annotation(tag_levels = c('A', '1'),
#                   title = "Standardised residuals in the training and hold-out cohorts",
#                   theme = theme(legend.position = "bottom")) # title of full plot
# 
# 
# pdf(width = 10, height = 9, "Plots/11.08.plot_4.pdf")
# plot_4
# dev.off()
# 
#
# #:-----------------------------------------------------------------------------
# # Variable importance
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# # load BCF model
# bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
# 
# # load development dataset
# hba1c.train.cleaned_up <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
#   # drop the variables with the most missingness
#   select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
#   # only complete cases
#   drop_na() %>%
#   as.data.frame() %>%
#   select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau))) %>%
#   cbind(effects = colMeans(bcf_model$tau))
# 
# 
# m1 <- rms::ols(effects ~ rms::rcs(agetx, 3) + sex + ncurrtx + rms::rcs(prehba1c, 3) + rms::rcs(prebmi, 3) + rms::rcs(preegfr, 3) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = hba1c.train.cleaned_up, x = TRUE, y = TRUE)
# 
# values <- plot(anova(m1), what = 'proportion R2')
# 
# plot_5 <- as.data.frame(values) %>%
#   cbind(variable = c("Number of other current glucose-lowering drugs", "Sex", "eGFR", "Current age", "BMI", "HbA1c", "Retinopathy", "Peripheral arterial disease", "Neuropathy", "Ischaemic heart disease", "Heart failure")) %>%
#   mutate(variable = factor(variable),
#          values = values * 100) %>%
#   ggplot(aes(y = forcats::fct_reorder(variable, values), x = values)) +
#   geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(variable, values)), linetype = "dashed") +
#   geom_point(size = 2, colour = "black") +
#   ggtitle("Relative importance for treatment effect heterogeneity") +
#   xlab("Relative Importance (%)") +
#   theme_bw() +
#   theme(axis.text.y = element_text(angle = 45),
#         axis.title.y = element_blank())
# 
# 
# pdf("Plots/11.08.plot_5.pdf")
# plot_5
# dev.off()
# 
# 
# #:------------------------------------------------------------------------------
# # Heterogeneity in sex
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# # load BCF model
# bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
# 
# # load development dataset
# hba1c.train.heterogeneity <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
#   # drop the variables with the most missingness
#   select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
#   # only complete cases
#   drop_na() %>%
#   as.data.frame() %>%
#   select(sex) %>%
#   cbind(mean = colMeans(bcf_model$tau))
# 
# plot_effects_female <- hist_plot(hba1c.train.heterogeneity %>% filter(sex == "Female"), "Females", -15, 18) +
#   theme(legend.direction = "horizontal", 
#         legend.position = "bottom",
#         legend.box = "horizontal")
# 
# plot_effects_male <- hist_plot(hba1c.train.heterogeneity %>% filter(sex == "Male"), "Males", -15, 18) +
#   theme(legend.direction = "horizontal", 
#         legend.position = "bottom",
#         legend.box = "horizontal")
# 
# plot_6 <- patchwork::wrap_plots(list(
#   plot_effects_female,
#   plot_effects_male
# ), ncol = 2) +
#   patchwork::plot_layout(guides = "collect") +
#   patchwork::plot_annotation(tag_levels = "A",
#                              title = "Heterogeneity of treatment effects",  # title of full plot
#                              theme = theme(legend.direction = "horizontal", 
#                                            legend.position = "bottom",
#                                            legend.box = "horizontal"))
# 
# 
# pdf(width = 8, height = 5, "Plots/11.08.plot_6.pdf")
# plot_6
# dev.off()
# 
# 
# #:------------------------------------------------------------------------------
# # Validation of BCF model (population with/without CVD/HF/CKD)
# 
# patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")
# 
# patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# # Population without CVD/HF/CKD
# no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   left_join(patient_effects, by = c("patid", "pated"))
# 
# predicted_no_co.dataset <- no_co.dataset %>%
#   rename("hba1c_diff" = "effects") %>%
#   mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
#          hba1c_diff.q = ntile(hba1c_diff, 10)) %>%
#   drop_na(posthba1cfinal, hba1c_diff)
# 
# 
# # Population with CVD/HF/CKD
# full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")
# 
# co.dataset <- full.cohort[!(full.cohort$pated %in% no_co.dataset$pated),] %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   left_join(patient_effects, by = c("patid", "pated"))
# 
# predicted_co.dataset <- co.dataset %>%
#   rename("hba1c_diff" = "effects") %>%
#   mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
#          hba1c_diff.q = ntile(hba1c_diff, 10)) %>%
#   drop_na(posthba1cfinal, hba1c_diff)
# 
# 
# ##:-- No comorbidities
# # PSM 1:1
# ATE_matching_1_1_no_co.dataset <- calc_ATE_validation_prop_matching(predicted_no_co.dataset, "posthba1cfinal", predicted_no_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
# 
# plot_ATE_matching_1_1_no_co.dataset <- ATE_plot(ATE_matching_1_1_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # PSM 1:1 + adjusted
# ATE_matching_1_1_adjusted_no_co.dataset <- calc_ATE_validation_prop_matching(predicted_no_co.dataset, "posthba1cfinal", predicted_no_co.dataset$prop.score, order = "largest", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)
# 
# plot_ATE_matching_1_1_adjusted_no_co.dataset <- ATE_plot(ATE_matching_1_1_adjusted_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # Adjusted
# ATE_adjusted_no_co.dataset <- calc_ATE_validation_adjust(predicted_no_co.dataset, "posthba1cfinal", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)
# 
# plot_ATE_adjusted_no_co.dataset <- ATE_plot(ATE_adjusted_no_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# ##:-- Comorbidities
# # PSM 1:1
# ATE_matching_1_1_co.dataset <- calc_ATE_validation_prop_matching(predicted_co.dataset, "posthba1cfinal", predicted_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
# 
# plot_ATE_matching_1_1_co.dataset <- ATE_plot(ATE_matching_1_1_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # PSM 1:1 + adjusted
# ATE_matching_1_1_adjusted_co.dataset <- calc_ATE_validation_prop_matching(predicted_co.dataset, "posthba1cfinal", predicted_co.dataset$prop.score, order = "largest", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)
# 
# plot_ATE_matching_1_1_adjusted_co.dataset <- ATE_plot(ATE_matching_1_1_adjusted_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# # Adjusted
# ATE_adjusted_co.dataset <- calc_ATE_validation_adjust(predicted_co.dataset, "posthba1cfinal", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)
# 
# plot_ATE_adjusted_co.dataset <- ATE_plot(ATE_adjusted_co.dataset[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci")
# 
# 
# # No comorbidities plots
# plot_7.1 <- patchwork::wrap_plots(list(plot_ATE_matching_1_1_no_co.dataset + ggtitle("Match 1:1"),
#                                        plot_ATE_matching_1_1_adjusted_no_co.dataset + ggtitle("Match 1:1 adjusted"),
#                                        plot_ATE_adjusted_no_co.dataset + ggtitle("Adjusted")), ncol = 3) +
#   patchwork::plot_layout(guides = "collect") +
#   patchwork::plot_annotation(tag_levels = "A",
#                              title = paste0("No comorbidities cohort: Treatment effect validation (n=", predicted_no_co.dataset%>%nrow(),")"), # title of full plot
#                              theme = theme(plot.title = element_text(hjust = 0.5),
#                                            legend.position = "bottom")) # center title of full plot
# 
# # Comorbidities plots
# plot_7.2 <- patchwork::wrap_plots(list(plot_ATE_matching_1_1_co.dataset + ggtitle("Match 1:1"),
#                                        plot_ATE_matching_1_1_adjusted_co.dataset + ggtitle("Match 1:1 adjusted"),
#                                        plot_ATE_adjusted_co.dataset + ggtitle("Adjusted")), ncol = 3) +
#   patchwork::plot_layout(guides = "collect") +
#   patchwork::plot_annotation(tag_levels = "A",
#                              title = paste0("Comorbidities cohort: Treatment effect validation (n=", predicted_co.dataset%>%nrow(),")"), # title of full plot
#                              theme = theme(plot.title = element_text(hjust = 0.5),
#                                            legend.position = "bottom")) # center title of full plot
# 
# pdf(width = 12, height = 5, "Plots/11.08.plot_7.pdf")
# plot_7.1
# plot_7.2
# dev.off()
# 
# 
# #:------------------------------------------------------------------------------
# # HbA1c grouping (population with/without CVD/HF/CKD)
# 
# require(forestplot)
# 
# patient_prop_scores <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/patient_prop_scores.rds")
# 
# patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# interval_breaks <- c(-5, -3, 0, 3, 5)
# 
# # Population without CVD/HF/CKD
# no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   left_join(patient_effects, by = c("patid", "pated"))
# 
# group.hba1c.no_co.dataset <- group_values(data = no_co.dataset,
#                                           variable = "effects",
#                                           breaks = interval_breaks) %>%
#   drop_na(intervals, posthba1cfinal) %>%
#   rename("hba1c_diff" = "effects")
# 
# 
# # Population with CVD/HF/CKD
# full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")
# 
# co.dataset <- full.cohort[!(full.cohort$pated %in% no_co.dataset$pated),] %>%
#   left_join(patient_prop_scores, by = c("patid", "pated")) %>%
#   left_join(patient_effects, by = c("patid", "pated"))
# 
# group.hba1c.co.dataset <- group_values(data = co.dataset,
#                                        variable = "effects",
#                                        breaks = interval_breaks) %>%
#   drop_na(intervals, posthba1cfinal) %>%
#   rename("hba1c_diff" = "effects")
# 
# 
# ##:-- No comorbidities
# # PSM 1:1
# ATE_matching_1_1_no_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.no_co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
# 
# # PSM 1:1 + adjusted
# ATE_matching_1_1_adjusted_no_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.no_co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)
# 
# # Adjusted
# ATE_adjusted_no_co.dataset <- calc_ATE_validation_adjust(group.hba1c.no_co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preneuropathy", "preretinopathy", "t2dmduration", "drugline", "hba1cmonth", "prealt"), adjust = TRUE)
# 
# ##:-- Comorbidities
# # PSM 1:1
# ATE_matching_1_1_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)))
# 
# # PSM 1:1 + adjusted
# ATE_matching_1_1_adjusted_co.dataset <- calc_ATE_validation_prop_matching(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.co.dataset$prop.score, quantile_var = "intervals", order = "largest", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)
# 
# # Adjusted
# ATE_adjusted_co.dataset <- calc_ATE_validation_adjust(group.hba1c.co.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_tau, variables_mu)), adjust = TRUE)
# 
# 
# hba1c_overall_axis_min <- plyr::round_any(floor(min(c(ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
#                                                       ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
#                                                       ATE_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
#                                                       ATE_matching_1_1_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
#                                                       ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
#                                                       ATE_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)
# 
# hba1c_overall_axis_max <- plyr::round_any(ceiling(max(c(ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
#                                                         ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
#                                                         ATE_adjusted_no_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
#                                                         ATE_matching_1_1_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
#                                                         ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
#                                                         ATE_adjusted_co.dataset[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)
# 
# 
# # Propensity score matching
# plot_psm_1_1_no_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("No comorbidities population (n=", sum(ATE_matching_1_1_no_co.dataset[["effects"]]$n_drug1)*2, ")"))
# 
# plot_psm_1_1_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("Comorbidities population (n=", sum(ATE_matching_1_1_co.dataset[["effects"]]$n_drug1)*2, ")"))
# 
# 
# # Propensity score matching + adjusted
# plot_psm_1_1_adjusted_no_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching + adjusted",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("No comorbidities population (n=", sum(ATE_matching_1_1_adjusted_no_co.dataset[["effects"]]$n_drug1)*2, ")"))
# 
# plot_psm_1_1_adjusted_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_matching_1_1_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_matching_1_1_adjusted_co.dataset[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching + adjusted",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("Comorbidities population (n=", sum(ATE_matching_1_1_adjusted_co.dataset[["effects"]]$n_drug1)*2, ")"))
# 
# 
# # Adjusted
# plot_adjusted_no_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_adjusted_no_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(1)%>%unlist(), ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(2)%>%unlist(), ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(3)%>%unlist(), ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(4)%>%unlist(), ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(5)%>%unlist(), ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_adjusted_no_co.dataset[["effects"]]%>%select(N)%>%slice(6)%>%unlist(), ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching + adjusted",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("No comorbidities population (n=", sum(ATE_adjusted_no_co.dataset[["effects"]]$N), ")"))
# 
# plot_adjusted_co_hba1c <- rbind(
#   cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
#   ATE_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
#   cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
#   ATE_adjusted_co.dataset[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6)
# ) %>%
#   as.data.frame() %>%
#   mutate(obs = as.numeric(obs),
#          lci = as.numeric(lci),
#          uci = as.numeric(uci),
#          intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(1)%>%unlist(), ")"), 
#                             ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(2)%>%unlist(), ")"), 
#                                    ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(3)%>%unlist(), ")"), 
#                                           ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(4)%>%unlist(), ")"), 
#                                                  ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(5)%>%unlist(), ")"),
#                                                         ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_adjusted_co.dataset[["effects"]]%>%select(N)%>%slice(6)%>%unlist(), ")"), intervals))))))) %>%
#   rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
#   forestplot(labeltext = intervals,
#              ci.vertices = TRUE,
#              title = "Propensity score matching + adjusted",
#              xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
#              ci.vertices.height = 0.05,
#              boxsize = .1,
#              txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
#              xlab = "Average treatment effect (mmol/mol)") %>%
#   fp_add_header(paste0("Comorbidities population (n=", sum(ATE_adjusted_co.dataset[["effects"]]$N), ")"))
# 
# 
# 
# pdf(width = 15, height = 12, "Plots/11.08.plot_8.pdf")
# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(nrow = 4,
#                                            ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# # title
# grid.text("HbA1c treatment effect", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
# 
# # first plot
# pushViewport(viewport(layout.pos.row = 2,
#                       layout.pos.col = 1))
# plot_psm_1_1_no_co_hba1c
# upViewport()
# # second plot
# pushViewport(viewport(layout.pos.row = 2,
#                       layout.pos.col = 2))
# plot_psm_1_1_co_hba1c
# upViewport()
# # third plot
# pushViewport(viewport(layout.pos.row = 3,
#                       layout.pos.col = 1))
# plot_psm_1_1_adjusted_no_co_hba1c
# upViewport()
# # forth plot
# pushViewport(viewport(layout.pos.row = 3,
#                       layout.pos.col = 2))
# plot_psm_1_1_adjusted_co_hba1c
# upViewport()
# # fifth plot
# pushViewport(viewport(layout.pos.row = 4,
#                       layout.pos.col = 1))
# plot_adjusted_no_co_hba1c
# upViewport()
# # sixth plot
# pushViewport(viewport(layout.pos.row = 4,
#                       layout.pos.col = 2))
# plot_adjusted_co_hba1c
# upViewport()
# 
# dev.off()
# 
# 
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
hba1c <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  rbind(
    set_up_data_sglt2_glp1(dataset.type = "hba1c.test") %>%
      left_join(patient_prop_scores, by = c("patid", "pated")) %>%
      left_join(treatment_effects, by = c("patid", "pated"))
  )

group.hba1c.dataset <- group_values(data = hba1c,
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
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.hba1c.dataset%>%filter(intervals == levels(group.hba1c.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
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
  fp_add_header(paste0("Overall population (n=", format(group.hba1c.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
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
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 1, f = floor)

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
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.weight.dataset%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Weight change",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 1),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", format(group.weight.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
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
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.discontinuation.dataset%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals)))))),
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
  fp_add_header(paste0("Overall population (n=", format(group.discontinuation.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


# No microvascular complications

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

micro_comp.dataset <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.micro_comp.dataset <- group_values(data = micro_comp.dataset,
                                         variable = "effects",
                                         breaks = interval_breaks) %>%
  drop_na(intervals)

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_overall.rds"))

# predictions for the CVD outcomes in the population with no CVD/HF/CKD
predictions_micro_comp_micro_comp_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_micro_comp_micro_comp_stan_adjusted_full.rds"))


# Adjusted
plot_micro_comp_micro_comp_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_micro_comp_micro_comp_stan_adjusted_overall %>% slice(4:6),
  predictions_micro_comp_micro_comp_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.micro_comp.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[1])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.micro_comp.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[2])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.micro_comp.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[3])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.micro_comp.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[4])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.micro_comp.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[5])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.micro_comp.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%nrow(), ", event = ", group.micro_comp.dataset%>%filter(intervals==levels(group.micro_comp.dataset$intervals)[6])%>%filter(postdrug_micro_comp_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Microvascular complications",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.micro_comp.dataset), ")")) %>%
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
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
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
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
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
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", event = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
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
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



pdf(width = 20, height = 10, "Plots/11.08.plot_9.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 3, heights = unit(c(0.2, 5, 0.2, 5), "null"))))
# 1 2 3
# 4 5 6
grid.text(expression(bol("A.1")), vp = viewpost(layout.pos.row = 1, layout.pos.col = 1),just = "right")
grid.text(expression(bol("A.2")), vp = viewpost(layout.pos.row = 1, layout.pos.col = 2),just = "right")
grid.text(expression(bol("A.3")), vp = viewpost(layout.pos.row = 1, layout.pos.col = 3),just = "right")
# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_adjust_overall_hba1c
upViewport()

# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_weight_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 3))
plot_discontinuation_adjusted_overall
upViewport()

grid.text(expression(bol("B.1")), vp = viewpost(layout.pos.row = 3, layout.pos.col = 1),just = "right")
grid.text(expression(bol("B.2")), vp = viewpost(layout.pos.row = 3, layout.pos.col = 2),just = "right")
grid.text(expression(bol("B.3")), vp = viewpost(layout.pos.row = 3, layout.pos.col = 3),just = "right")
# forth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_micro_comp_micro_comp_adjusted_overall
upViewport()

# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_no_co_cvd_adjusted_overall
upViewport()

# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 3))
plot_no_co_hf_adjusted_overall
upViewport()

dev.off()





#:------------------------------------------------------------------------------
# Summary of BCF model: a) histogram of treatment effects in development,
#   b) decile calibration of treatment effects (development)
#   c) decile calibration of treatment effects

# Validation plot dev
ATE_adjust_validation_dev <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_dev.rds")

plot_ATE_adjust_validation_dev <- ATE_plot(ATE_adjust_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste0("Development cohort (n=", format(sum(ATE_adjust_validation_dev[["effects"]]$N),big.mark=",",scientific=FALSE), ")")) +
  scale_y_continuous(label=comma) +
  xlim(-12, 12) +
  ylim(-12, 12)

# histogram dev
bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")

data_dev <- cbind(mean = colMeans(bcf_model$tau)) %>%
  as.data.frame()

plot_effect_dev <- hist_plot(data_dev, paste0("Development cohort (n=", format(sum(ATE_adjust_validation_dev[["effects"]]$N),big.mark=",",scientific=FALSE), ")"), -15, 20) +
  scale_y_continuous(label=comma)


# Validation plot val
ATE_adjust_validation_val <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/ATE_adjust_validation_val.rds")

plot_ATE_adjust_validation_val <- ATE_plot(ATE_adjust_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -12, 12, colour_background = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste0("Validation cohort (n=", format(sum(ATE_adjust_validation_val[["effects"]]$N),big.mark=",",scientific=FALSE), ")")) +
  scale_y_continuous(label=comma) +
  xlim(-12, 12) +
  ylim(-12, 12)


plot_10.1 <- plot_effect_dev | plot_ATE_adjust_validation_dev | plot_ATE_adjust_validation_val

plot_10 <- plot_10.1 +
  plot_annotation(tag_levels = "A",
                  title = "Model development and validation",
                  theme = theme(legend.position = "bottom")) # title of full plot

pdf(width = 15, height = 6, "Plots/11.08.plot_10.pdf")
plot_10
dev.off()

# 
# 
# #:------------------------------------------------------------------------------
# # Model explanability: a) relative importance for treatment effect,
# #   b) distribution of treatment effects by sex
# #   c) differential treatment effects
# 
# ### a)
# 
# # load variables used in the BCF model
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")
# 
# # load BCF model
# bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/bcf_model.rds")
# 
# # load development dataset
# hba1c.train.cleaned_up <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
#   # drop the variables with the most missingness
#   select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
#   # only complete cases
#   drop_na() %>%
#   as.data.frame() %>%
#   select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau))) %>%
#   cbind(effects = colMeans(bcf_model$tau))
# 
# 
# m1 <- rms::ols(effects ~ rms::rcs(agetx, 3) + sex + ncurrtx + rms::rcs(prehba1c, 3) + rms::rcs(prebmi, 3) + rms::rcs(preegfr, 3) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = hba1c.train.cleaned_up, x = TRUE, y = TRUE)
# 
# values <- plot(anova(m1), what = 'proportion R2')
# 
# plot_10.1 <- as.data.frame(values) %>%
#   cbind(variable = c("Number of other current glucose-lowering drugs", "Sex", "eGFR", "Current age", "BMI", "HbA1c", "Retinopathy", "Peripheral arterial disease", "Neuropathy", "Ischaemic heart disease", "Heart failure")) %>%
#   mutate(variable = factor(variable),
#          values = values * 100) %>%
#   ggplot(aes(y = forcats::fct_reorder(variable, values), x = values)) +
#   geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(variable, values)), linetype = "dashed") +
#   geom_point(size = 2, colour = "black") +
#   ggtitle("Relative importance for treatment effect heterogeneity") +
#   xlab("Relative Importance (%)") +
#   theme_bw() +
#   theme(axis.text.y = element_text(angle = 45),
#         axis.title.y = element_blank())
# 
# ### b)
# 
# 
# hba1c.heterogeneity <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
#   left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
#   select(sex, effects) %>%
#   drop_na() %>%
#   rename("mean" = "effects")
# 
# 
# percent_var = round(((hba1c.heterogeneity %>% filter(sex == "Female" & mean < 0) %>% nrow()) / hba1c.heterogeneity %>% filter(sex == "Female") %>% nrow())*100)
# 
# plot_10.2 <- hist_plot(hba1c.heterogeneity %>% filter(sex == "Female"), "Females", -15, 18) +
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom",
#         legend.box = "horizontal") +
#   labs(subtitle = paste0("SGLT2i = ",
#                          percent_var,
#                          "%, GLP1-RA = ",
#                          100-percent_var,
#                          "%"))
# 
# percent_var = round(((hba1c.heterogeneity %>% filter(sex == "Male" & mean < 0) %>% nrow()) / hba1c.heterogeneity %>% filter(sex == "Male") %>% nrow())*100)
# 
# plot_10.3 <- hist_plot(hba1c.heterogeneity %>% filter(sex == "Male"), "Males", -15, 18) +
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom",
#         legend.box = "horizontal") +
#   labs(subtitle = paste0("SGLT2i = ",
#                          percent_var,
#                          "%, GLP1-RA = ",
#                          100-percent_var,
#                          "%"))
# 
# ### c)
# 
# # treatment effects
# patient_effects <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds")
# 
# # Full cohort for average values
# hba1c.train <- set_up_data_sglt2_glp1(dataset.type="hba1c.train") %>%
#   left_join(patient_effects, by = c("patid", "pated"))
# 
# levels(hba1c.train$sex) <- c("Females", "Males")
# levels(hba1c.train$ncurrtx) <- c("0", "1", "2", "3", "4+")
# 
# #:------------------------------------------------------------------------
# # Stratify data by variables
# # sex
# plot_sex_strata <- hba1c.train %>%
#   select(sex, effects) %>%
#   group_by(sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = sex, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.15) +
#   ylim(-10, 10) +
#   ggtitle("Sex") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5))
# 
# # ncurrtx
# plot_ncurrtx_strata <- hba1c.train %>%
#   select(ncurrtx, effects, sex) %>%
#   group_by(ncurrtx, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = ncurrtx, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.5) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Number of other current\nglucose-lowering drugs") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # prepad
# plot_prepad_strata <- hba1c.train %>%
#   select(prepad, effects, sex) %>%
#   group_by(prepad, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = prepad, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Peripheral arterial disease") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # preihd
# plot_preihd_strata <- hba1c.train %>%
#   select(preihd, effects, sex) %>%
#   group_by(preihd, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = preihd, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Ischaemic heart disease") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # preneuropathy
# plot_preneuropathy_strata <- hba1c.train %>%
#   select(preneuropathy, effects, sex) %>%
#   group_by(preneuropathy, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = preneuropathy, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Neuropathy") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # preretinopathy
# plot_preretinopathy_strata <- hba1c.train %>%
#   select(preretinopathy, effects, sex) %>%
#   group_by(preretinopathy, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = preretinopathy, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Retinopathy") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # preheartfailure
# plot_preheartfailure_strata <- hba1c.train %>%
#   select(preheartfailure, effects, sex) %>%
#   group_by(preheartfailure, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   ggplot(aes(x = preheartfailure, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Heart failure") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_rect(fill="white"))
# 
# # prehba1c
# breaks_hba1c <- quantile(hba1c.train$prehba1c, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
# 
# plot_prehba1c_strata <- group_values(data = hba1c.train,
#                                      variable = "prehba1c",
#                                      breaks = breaks_hba1c) %>%
#   select(intervals, effects, sex) %>%
#   group_by(intervals, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   drop_na() %>%
#   ggplot(aes(x = intervals, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("HbA1c") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust=1),
#         strip.background = element_rect(fill="white"))
# 
# # preegfr
# breaks_egfr <- quantile(hba1c.train$preegfr, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
# 
# plot_preegfr_strata <- group_values(data = hba1c.train,
#                                     variable = "preegfr",
#                                     breaks = breaks_egfr) %>%
#   select(intervals, effects, sex) %>%
#   group_by(intervals, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   drop_na() %>%
#   ggplot(aes(x = intervals, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("eGFR") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust=1),
#         strip.background = element_rect(fill="white"))
# 
# # agetx
# breaks_agetx <- quantile(hba1c.train$agetx, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
# 
# plot_agetx_strata <- group_values(data = hba1c.train,
#                                   variable = "agetx",
#                                   breaks = breaks_agetx) %>%
#   select(intervals, effects, sex) %>%
#   group_by(intervals, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   drop_na() %>%
#   ggplot(aes(x = intervals, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("Current age") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust=1),
#         strip.background = element_rect(fill="white"))
# 
# # prebmi
# breaks_prebmi <- quantile(hba1c.train$prebmi, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
# 
# plot_prebmi_strata <- group_values(data = hba1c.train,
#                                    variable = "prebmi",
#                                    breaks = breaks_prebmi) %>%
#   select(intervals, effects, sex) %>%
#   group_by(intervals, sex) %>%
#   mutate(mean = mean(effects,  na.rm = TRUE),
#          lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
#          uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
#   select(-effects) %>%
#   ungroup() %>%
#   unique() %>%
#   drop_na() %>%
#   ggplot(aes(x = intervals, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
#   geom_point() +
#   facet_wrap(~sex) +
#   ylim(-10, 10) +
#   ggtitle("BMI") +
#   ylab("Predicted treatment effects (mmol/mol)") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.title = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust=1),
#         strip.background = element_rect(fill="white"))
# 
# ### combine
# 
# pdf(width = 17, height = 12, "Plots/11.08.plot_11.pdf")
# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(nrow = 8,
#                                            ncol = 12, heights = unit(c(0.15, 0.2, 1.5, 0.1, 0.1, 1, 1, 1), "null"))))
# # title
# grid.text("Model explainability", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:12))
# 
# # legend
# grid.text(expression(bold("A")), vp = viewport(layout.pos.row = 2, layout.pos.col = 1), just = "right")
# grid.text(expression(bold("B.1")), vp = viewport(layout.pos.row = 2, layout.pos.col = 7), just = "right")
# grid.text(expression(bold("B.2")), vp = viewport(layout.pos.row = 2, layout.pos.col = 10), just = "right")
# 
# # first plot
# pushViewport(viewport(layout.pos.row = 3,
#                       layout.pos.col = 1:6))
# grid.draw(ggplotGrob(plot_10.1))
# upViewport()
# 
# # second plot
# pushViewport(viewport(layout.pos.row = 3,
#                       layout.pos.col = 7:9))
# grid.draw(ggplotGrob(plot_10.2))
# upViewport()
# 
# # third plot
# pushViewport(viewport(layout.pos.row = 3,
#                       layout.pos.col = 10:12))
# grid.draw(ggplotGrob(plot_10.3))
# upViewport()
# 
# # legend
# grid.text(expression(bold("C")), vp = viewport(layout.pos.row = 5, layout.pos.col = 1), just = "right")
# 
# 
# # forth plot
# pushViewport(viewport(layout.pos.row = 6,
#                       layout.pos.col = 1:3))
# grid.draw(ggplotGrob(plot_ncurrtx_strata +
#                        ylab("") +
#                        theme(axis.title.y = element_text())))
# upViewport()
# 
# # fifth plot
# pushViewport(viewport(layout.pos.row = 6,
#                       layout.pos.col = 4:6))
# grid.draw(ggplotGrob(plot_sex_strata))
# upViewport()
# 
# # sixth plot
# pushViewport(viewport(layout.pos.row = 6,
#                       layout.pos.col = 7:9))
# grid.draw(ggplotGrob(plot_preegfr_strata))
# upViewport()
# 
# # seventh plot
# pushViewport(viewport(layout.pos.row = 6,
#                       layout.pos.col = 10:12))
# grid.draw(ggplotGrob(plot_prehba1c_strata))
# upViewport()
# 
# # nineth plot
# pushViewport(viewport(layout.pos.row = 7,
#                       layout.pos.col = 4:6))
# grid.draw(ggplotGrob(plot_prehba1c_strata))
# upViewport()
# 
# # tenth plot
# pushViewport(viewport(layout.pos.row = 7,
#                       layout.pos.col = 7:9))
# grid.draw(ggplotGrob(plot_preretinopathy_strata))
# upViewport()
# 
# # eleventh plot
# pushViewport(viewport(layout.pos.row = 7,
#                       layout.pos.col = 10:12))
# grid.draw(ggplotGrob(plot_prepad_strata))
# upViewport()
# 
# # twelfth plot
# pushViewport(viewport(layout.pos.row = 8,
#                       layout.pos.col = 1:3))
# grid.draw(ggplotGrob(plot_preneuropathy_strata +
#                        ylab("") +
#                        theme(axis.title.y = element_text())))
# upViewport()
# 
# # thirteenth plot
# pushViewport(viewport(layout.pos.row = 8,
#                       layout.pos.col = 4:6))
# grid.draw(ggplotGrob(plot_preihd_strata))
# upViewport()
# 
# # fourteenth plot
# pushViewport(viewport(layout.pos.row = 8,
#                       layout.pos.col = 7:9))
# grid.draw(ggplotGrob(plot_preheartfailure_strata))
# upViewport()
# 
# # eighth plot
# pushViewport(viewport(layout.pos.row = 7,
#                       layout.pos.col = 1:3))
# grid.draw(ggplotGrob(plot_prebmi_strata +
#                        theme(axis.title.y = element_text(size = 12))))
# upViewport()
# 
# 
# dev.off()
# 
# 
# 
# 
# #:------------------------------------------------------------------------------
# ## Variable selection of sparcebcf
# variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_tau_original_chain_1.rds")
# 
# plot_variables_tau <- variables_tau %>%
#   rename("Current age" = "agetx",
#          "Sex" = "sex",
#          "Duration of diabetes" = "t2dmduration",
#          "Number of glucose-lowering drug classes ever prescribed" = "drugline",
#          "Number of other current glucose-lowering drugs" = "ncurrtx",
#          "Month of Hba1C measure" = "hba1cmonth",
#          "HbA1c" = "prehba1c",
#          "BMI" = "prebmi",
#          "eGFR" = "preegfr",
#          "Albumin" = "prealbuminblood",
#          "Alanine transaminase" = "prealt",
#          "Bilirubin" = "prebilirubin",
#          "HDL" = "prehdl",
#          "Mean arterial blood pressure" = "premap",
#          "Total cholesterol" = "pretotalcholesterol",
#          "Angina" = "preangina",
#          "Chronic liver disease" = "precld",
#          "Nephropathy" = "prediabeticnephropathy",
#          "Heart failure" = "preheartfailure",
#          "Hypertension" = "prehypertension",
#          "Ischaemic heart disease" = "preihd",
#          "Myocardial infarction" = "premyocardialinfarction",
#          "Neuropathy" = "preneuropathy",
#          "Peripheral arterial disease" = "prepad",
#          "Retinopathy" = "preretinopathy",
#          "Atherosclerotic cardiovascular disease" = "prerevasc",
#          "Stroke" = "prestroke",
#          "Transient ischaemic attack" = "pretia",
#          "Atrial fibrillation" = "preaf") %>% 
#   as.data.frame() %>%
#   gather(key, value) %>%
#   arrange(desc(value)) %>%
#   mutate(key = factor(key),
#          colour = ifelse(value > 0.023, "Above", "Below")) %>%
#   ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
#   geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
#   geom_point(size = 2) +
#   geom_vline(aes(xintercept = 0.023), colour = "red") +
#   ggtitle("Moderator model") +
#   xlab("Posterior splitting probabilities") +
#   scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
#   theme_light() +
#   theme(axis.text.y = element_text(angle = 30),
#         axis.title.y = element_blank(),
#         legend.position = "none")
# 
# 
# variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_mu_original_chain_1.rds")
# 
# plot_variables_mu <- variables_mu %>%
#   rename("Current age" = "agetx",
#          "Sex" = "sex",
#          "Duration of diabetes" = "t2dmduration",
#          "Number of glucose-lowering drug classes ever prescribed" = "drugline",
#          "Number of other current glucose-lowering drugs" = "ncurrtx",
#          "Month of Hba1C measure" = "hba1cmonth",
#          "HbA1c" = "prehba1c",
#          "BMI" = "prebmi",
#          "eGFR" = "preegfr",
#          "Albumin" = "prealbuminblood",
#          "Alanine transaminase" = "prealt",
#          "Bilirubin" = "prebilirubin",
#          "HDL" = "prehdl",
#          "Mean arterial blood pressure" = "premap",
#          "Total cholesterol" = "pretotalcholesterol",
#          "Angina" = "preangina",
#          "Chronic liver disease" = "precld",
#          "Nephropathy" = "prediabeticnephropathy",
#          "Heart failure" = "preheartfailure",
#          "Hypertension" = "prehypertension",
#          "Ischaemic heart disease" = "preihd",
#          "Myocardial infarction" = "premyocardialinfarction",
#          "Neuropathy" = "preneuropathy",
#          "Peripheral arterial disease" = "prepad",
#          "Retinopathy" = "preretinopathy",
#          "Atherosclerotic cardiovascular disease" = "prerevasc",
#          "Stroke" = "prestroke",
#          "Transient ischaemic attack" = "pretia",
#          "Atrial fibrillation" = "preaf",
#          "Propensity score" = "propensity score") %>% 
#   as.data.frame() %>%
#   gather(key, value) %>%
#   arrange(desc(value)) %>%
#   mutate(key = factor(key),
#          colour = ifelse(value > 0.023, "Above", "Below")) %>%
#   ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
#   geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
#   geom_point(size = 2) +
#   geom_vline(aes(xintercept = 0.023), colour = "red") +
#   ggtitle("Control model") +
#   xlab("Posterior splitting probabilities") +
#   scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
#   theme_light() +
#   theme(axis.text.y = element_text(angle = 30),
#         axis.title.y = element_blank(),
#         legend.position = "none")
# 
# 
# pdf(width = 14, height = 6, "Plots/11.08.plot_12.pdf")
# patchwork::wrap_plots(list(plot_variables_mu, plot_variables_tau), ncol = 2) +
#   patchwork::plot_annotation(title = "BART model variable selection") # title of full plot
# dev.off()
# 
# 
# #:------------------------------------------------------------------------------
# ## Trace plots for sparsebcf and bcf
# 
# # Sparse BCF
# 
# sigma_sparsebcf <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/sigma_sparsebcf_vector.rds")
#   
# plot_sigma_sparsebcf <- ggplot() +
#   geom_path(aes(x = 1:length(sigma_sparsebcf), y = sigma_sparsebcf)) +
#   ggtitle("Sparse BCF") +
#   ylab("Sigma") +
#   theme_light() +
#   theme(axis.title.x = element_blank()) +
#   scale_x_continuous(label=comma)
# 
# 
# # BCF - propensity score
# 
# sigma_bcf_model_prop <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/sigma_bcf_model_prop_vector.rds")
# 
# plot_sigma_prop_bcf <- ggplot() +
#   geom_path(aes(x = 1:length(sigma_bcf_model_prop), y = sigma_bcf_model_prop)) +
#   ggtitle("BCF with propensity score in Control model") +
#   ylab("Sigma") +
#   theme_light() +
#   theme(axis.title.x = element_blank()) +
#   scale_x_continuous(label=comma)
# 
# 
# # BCF - no propensity score
# 
# sigma_bcf_model <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/sigma_bcf_model_vector.rds")
# 
# plot_sigma_bcf <- ggplot() +
#   geom_path(aes(x = 1:length(sigma_bcf_model), y = sigma_bcf_model)) +
#   ggtitle("BCF without propensity score in Control model") +
#   ylab("Sigma") +
#   theme_light() +
#   theme(axis.title.x = element_blank()) +
#   scale_x_continuous(label=comma)
# 
# 
# 
# pdf(width = 14, height = 6, "Plots/11.08.plot_13.pdf")
# patchwork::wrap_plots(list(plot_sigma_sparsebcf, plot_sigma_prop_bcf, plot_sigma_bcf), ncol = 1, row = 3)
# dev.off()
# 




















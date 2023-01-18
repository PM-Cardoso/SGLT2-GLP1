####################
## Description:
##  - In this file we calculate the treatment effect / response for 
##      the average patient for a range of values at each covariate.
####################

# library(bcf)
library(tidyverse)

## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1/Aurum"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path, "/differential_response"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

# variables chosen
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))
variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

# treatment effects
patient_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# Full cohort for average values
hba1c.train <- set_up_data_sglt2_glp1(dataset.type="hba1c.train") %>%
  left_join(patient_effects, by = c("patid", "pated"))

levels(hba1c.train$sex) <- c("Females", "Males")
levels(hba1c.train$ncurrtx) <- c("0", "1", "2", "3", "4+")



#:------------------------------------------------------------------------
# Stratify data by variables
# sex
plot_sex_strata <- hba1c.train %>%
  select(sex, effects) %>%
  group_by(sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = sex, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_point() +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.15) +
  ylim(-10, 10) +
  ggtitle("Sex") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# ncurrtx
plot_ncurrtx_strata <- hba1c.train %>%
  select(ncurrtx, effects, sex) %>%
  group_by(ncurrtx, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = ncurrtx, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.5) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Number of other current\nglucose-lowering drugs") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# prepad
plot_prepad_strata <- hba1c.train %>%
  select(prepad, effects, sex) %>%
  group_by(prepad, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = prepad, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Peripheral arterial disease") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# preihd
plot_preihd_strata <- hba1c.train %>%
  select(preihd, effects, sex) %>%
  group_by(preihd, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = preihd, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Ischaemic heart disease") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# preneuropathy
plot_preneuropathy_strata <- hba1c.train %>%
  select(preneuropathy, effects, sex) %>%
  group_by(preneuropathy, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = preneuropathy, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Neuropathy") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# preretinopathy
plot_preretinopathy_strata <- hba1c.train %>%
  select(preretinopathy, effects, sex) %>%
  group_by(preretinopathy, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = preretinopathy, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Retinopathy") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# preheartfailure
plot_preheartfailure_strata <- hba1c.train %>%
  select(preheartfailure, effects, sex) %>%
  group_by(preheartfailure, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = preheartfailure, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.3) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Heart failure") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill="white"))

# prehba1c
breaks_hba1c <- quantile(hba1c.train$prehba1c, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)

plot_prehba1c_strata <- group_values(data = hba1c.train,
                              variable = "prehba1c",
                              breaks = breaks_hba1c) %>%
  select(intervals, effects, sex) %>%
  group_by(intervals, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("HbA1c") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_rect(fill="white"))
  
# preegfr
breaks_egfr <- quantile(hba1c.train$preegfr, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)

plot_preegfr_strata <- group_values(data = hba1c.train,
                              variable = "preegfr",
                              breaks = breaks_egfr) %>%
  select(intervals, effects, sex) %>%
  group_by(intervals, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("eGFR") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_rect(fill="white"))

# agetx
breaks_agetx <- quantile(hba1c.train$agetx, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)

plot_agetx_strata <- group_values(data = hba1c.train,
                             variable = "agetx",
                             breaks = breaks_agetx) %>%
  select(intervals, effects, sex) %>%
  group_by(intervals, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Current age") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_rect(fill="white"))

# prebmi
breaks_prebmi <- quantile(hba1c.train$prebmi, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)

plot_prebmi_strata <- group_values(data = hba1c.train,
                                  variable = "prebmi",
                                  breaks = breaks_prebmi) %>%
  select(intervals, effects, sex) %>%
  group_by(intervals, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("BMI") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_rect(fill="white"))

# hba1cmonth
breaks_hba1cmonth <- quantile(hba1c.train$hba1cmonth, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)

plot_hba1cmonth_strata <- group_values(data = hba1c.train,
                                       variable = "hba1cmonth",
                                       breaks = breaks_hba1cmonth) %>%
  select(intervals, effects, sex) %>%
  group_by(intervals, sex) %>%
  mutate(mean = mean(effects,  na.rm = TRUE),
         lci = quantile(effects, probs = c(0.25), na.rm = TRUE),
         uci = quantile(effects, probs = c(0.75),  na.rm = TRUE)) %>%
  select(-effects) %>%
  ungroup() %>%
  unique() %>%
  drop_na() %>%
  ggplot(aes(x = intervals, y = mean)) +
  geom_hline(aes(yintercept = 0), colour = "red") +
  geom_errorbar(aes(x = intervals, y = mean, ymax = uci, ymin = lci)) +
  geom_point() +
  facet_wrap(~sex) +
  ylim(-10, 10) +
  ggtitle("Month of outcome") +
  ylab("Predicted treatment effects (mmol/mol)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_rect(fill="white"))



# Combine all
plot_strata <- patchwork::wrap_plots(
  list(
    plot_ncurrtx_strata,
    plot_sex_strata,
    plot_preegfr_strata,
    plot_agetx_strata,
    plot_prebmi_strata +
      theme(axis.title.y = element_text(size = 12)),
    plot_prehba1c_strata,
    plot_preretinopathy_strata,
    plot_prepad_strata,
    plot_preneuropathy_strata,
    plot_preihd_strata,
    plot_preheartfailure_strata
)) +
  patchwork::plot_annotation(
    title = "IQR of treatment effects for covariate strata"
  )

plot_strata <- recordPlot(plot_strata)

saveRDS(plot_strata, paste0(output_path, "/differential_response/plot_strata.rds"))

pdf(width = 17, height = 9, "Plots/11.07.diff_treat_effect.pdf")
plot_strata
dev.off()



# #:------------------------------------------------------------------------
# # Predictions average patient
# 
# 
# # Make average patient with the seleted variables
# average.patient <- cbind(
#   drugline = 1,
#   ncurrtx = 2,
#   hba1cmonth = 12,
#   prehba1c = mean(hba1c.train$prehba1c, na.rm = TRUE),
#   preegfr = mean(hba1c.train$preegfr, na.rm = TRUE),
#   prepad = 1,
#   agetx = mean(hba1c.train$agetx, na.rm = TRUE),
#   preihd = 1,
#   preneuropathy = 1,
#   preretinopathy = 1,
#   preaf = 1
#   ) %>%
#   as.data.frame() %>%
#   mutate(
#     drugline = factor(drugline, levels = 1:length(levels(hba1c.train$drugline)), labels = levels(hba1c.train$drugline)),
#     ncurrtx = factor(ncurrtx, levels = 1:length(levels(hba1c.train$ncurrtx)), labels = levels(hba1c.train$ncurrtx)),
#     prepad = factor(prepad, levels = 1:length(levels(hba1c.train$prepad)), labels = levels(hba1c.train$prepad)),
#     preihd = factor(preihd, levels = 1:length(levels(hba1c.train$preihd)), labels = levels(hba1c.train$preihd)),
#     preneuropathy = factor(preneuropathy, levels = 1:length(levels(hba1c.train$preneuropathy)), labels = levels(hba1c.train$preneuropathy)),
#     preretinopathy = factor(preretinopathy, levels = 1:length(levels(hba1c.train$preretinopathy)), labels = levels(hba1c.train$preretinopathy)),
#     preaf = factor(preaf, levels = 1:length(levels(hba1c.train$preaf)), labels = levels(hba1c.train$preaf))
#   )
# 
# # drugline
# 
# current.average.patient.drugline <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$drugline)))) %>%
#   mutate(drugline = rep(levels(hba1c.train$drugline), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$drugline)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          drugline = factor(drugline, levels = levels(hba1c.train$drugline)))
# 
# predictions.drugline <- predict(object = bcf_model,
#                                x_predict_control = current.average.patient.drugline %>%
#                                  select(
#                                    all_of(variables_mu)
#                                  ) %>%
#                                  mutate_all(funs(as.numeric(.))) %>%
#                                  as.matrix(),
#                                x_predict_moderate = current.average.patient.drugline %>%
#                                  select(
#                                    all_of(variables_tau)
#                                  ) %>%
#                                  mutate_all(funs(as.numeric(.))) %>%
#                                  as.matrix(),
#                                pi_pred = rep(0.5, nrow(current.average.patient.drugline)),
#                                z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$drugline))), 2),
#                                save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_drugline <- cbind(
#   values = current.average.patient.drugline$drugline,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$drugline))), 2),
#   mean = apply(predictions.drugline$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.drugline$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.drugline$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient.drugline$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$drugline))), labels = levels(hba1c.train$drugline)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Drugline") +
#   ylim(-8, 12)
# 
# 
# 
# # ncurrtx
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$ncurrtx)))) %>%
#   mutate(ncurrtx = rep(levels(hba1c.train$ncurrtx), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$ncurrtx)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          ncurrtx = factor(ncurrtx, levels = levels(hba1c.train$ncurrtx)))
# 
# predictions.ncurrtx <- predict(object = bcf_model,
#                               x_predict_control = current.average.patient %>%
#                                 select(
#                                   all_of(variables_mu)
#                                 ) %>%
#                                 mutate_all(funs(as.numeric(.))) %>%
#                                 as.matrix(),
#                               x_predict_moderate = current.average.patient %>%
#                                 select(
#                                   all_of(variables_tau)
#                                 ) %>%
#                                 mutate_all(funs(as.numeric(.))) %>%
#                                 as.matrix(),
#                               pi_pred = rep(0.5, nrow(current.average.patient)),
#                               z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$ncurrtx))), 2),
#                               save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_ncurrtx <- cbind(
#   values = current.average.patient$ncurrtx,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$ncurrtx))), 2),
#   mean = apply(predictions.ncurrtx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.ncurrtx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.ncurrtx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$ncurrtx))), labels = levels(hba1c.train$ncurrtx)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Ncurrtx") +
#   ylim(-8, 12)
# 
# 
# # prepad
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$prepad)))) %>%
#   mutate(prepad = rep(levels(hba1c.train$prepad), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$prepad)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          prepad = factor(prepad, levels = levels(hba1c.train$prepad)))
# 
# predictions.prepad <- predict(object = bcf_model,
#                               x_predict_control = current.average.patient %>%
#                                 select(
#                                   all_of(variables_mu)
#                                 ) %>%
#                                 mutate_all(funs(as.numeric(.))) %>%
#                                 as.matrix(),
#                               x_predict_moderate = current.average.patient %>%
#                                 select(
#                                   all_of(variables_tau)
#                                 ) %>%
#                                 mutate_all(funs(as.numeric(.))) %>%
#                                 as.matrix(),
#                               pi_pred = rep(0.5, nrow(current.average.patient)),
#                               z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$prepad))), 2),
#                               save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_prepad <- cbind(
#   values = current.average.patient$prepad,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$prepad))), 2),
#   mean = apply(predictions.prepad$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.prepad$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.prepad$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$prepad))), labels = levels(hba1c.train$prepad)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Prepad") +
#   ylim(-8, 12)
# 
# 
# 
# # preihd
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$preihd)))) %>%
#   mutate(preihd = rep(levels(hba1c.train$preihd), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$preihd)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          preihd = factor(preihd, levels = levels(hba1c.train$preihd)))
# 
# predictions.preihd <- predict(object = bcf_model,
#                                      x_predict_control = current.average.patient %>%
#                                        select(
#                                          all_of(variables_mu)
#                                        ) %>%
#                                        mutate_all(funs(as.numeric(.))) %>%
#                                        as.matrix(),
#                                      x_predict_moderate = current.average.patient %>%
#                                        select(
#                                          all_of(variables_tau)
#                                        ) %>%
#                                        mutate_all(funs(as.numeric(.))) %>%
#                                        as.matrix(),
#                                      pi_pred = rep(0.5, nrow(current.average.patient)),
#                                      z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$preihd))), 2),
#                                      save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_preihd <- cbind(
#   values = current.average.patient$preihd,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$preihd))), 2),
#   mean = apply(predictions.preihd$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.preihd$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.preihd$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$preihd))), labels = levels(hba1c.train$preihd)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Preihd") +
#   ylim(-8, 12)
# 
# 
# 
# # preneuropathy
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$preneuropathy)))) %>%
#   mutate(preneuropathy = rep(levels(hba1c.train$preneuropathy), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$preneuropathy)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          preneuropathy = factor(preneuropathy, levels = levels(hba1c.train$preneuropathy)))
# 
# predictions.preneuropathy <- predict(object = bcf_model,
#                                       x_predict_control = current.average.patient %>%
#                                         select(
#                                           all_of(variables_mu)
#                                         ) %>%
#                                         mutate_all(funs(as.numeric(.))) %>%
#                                         as.matrix(),
#                                       x_predict_moderate = current.average.patient %>%
#                                         select(
#                                           all_of(variables_tau)
#                                         ) %>%
#                                         mutate_all(funs(as.numeric(.))) %>%
#                                         as.matrix(),
#                                       pi_pred = rep(0.5, nrow(current.average.patient)),
#                                       z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$preneuropathy))), 2),
#                                       save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_preneuropathy <- cbind(
#   values = current.average.patient$preneuropathy,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$preneuropathy))), 2),
#   mean = apply(predictions.preneuropathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.preneuropathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.preneuropathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$preneuropathy))), labels = levels(hba1c.train$preneuropathy)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Preneuropathy") +
#   ylim(-8, 12)
# 
# 
# # preretinopathy
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$preretinopathy)))) %>%
#   mutate(preretinopathy = rep(levels(hba1c.train$preretinopathy), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$preretinopathy)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          preretinopathy = factor(preretinopathy, levels = levels(hba1c.train$preretinopathy)))
# 
# predictions.preretinopathy <- predict(object = bcf_model,
#                              x_predict_control = current.average.patient %>%
#                                select(
#                                  all_of(variables_mu)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              x_predict_moderate = current.average.patient %>%
#                                select(
#                                  all_of(variables_tau)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              pi_pred = rep(0.5, nrow(current.average.patient)),
#                              z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$preretinopathy))), 2),
#                              save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_preretinopathy <- cbind(
#   values = current.average.patient$preretinopathy,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$preretinopathy))), 2),
#   mean = apply(predictions.preretinopathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.preretinopathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.preretinopathy$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$preretinopathy))), labels = levels(hba1c.train$preretinopathy)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Preretinopathy") +
#   ylim(-8, 12)
# 
# 
# 
# 
# # preaf
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(levels(hba1c.train$preaf)))) %>%
#   mutate(preaf = rep(levels(hba1c.train$preaf), 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(levels(hba1c.train$preaf)))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)),
#          preaf = factor(preaf, levels = levels(hba1c.train$preaf)))
# 
# predictions.preaf <- predict(object = bcf_model,
#                              x_predict_control = current.average.patient %>%
#                                select(
#                                  all_of(variables_mu)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              x_predict_moderate = current.average.patient %>%
#                                select(
#                                  all_of(variables_tau)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              pi_pred = rep(0.5, nrow(current.average.patient)),
#                              z_pred = rep(rep(c(0,1), each = length(levels(hba1c.train$preaf))), 2),
#                              save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_preaf <-  cbind(
#   values = current.average.patient$preaf,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(levels(hba1c.train$preaf))), 2),
#   mean = apply(predictions.preaf$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.preaf$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.preaf$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = c(1:length(levels(hba1c.train$preaf))), labels = levels(hba1c.train$preaf)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Preaf") +
#   ylim(-8, 12)
# 
# 
# # prehba1c
# # breaks_hba1c
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(breaks_hba1c))) %>%
#   mutate(prehba1c = rep(breaks_hba1c, 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(breaks_hba1c))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)))
# 
# predictions.prehba1c <- predict(object = bcf_model,
#                                x_predict_control = current.average.patient %>%
#                                  select(
#                                    all_of(variables_mu)
#                                  ) %>%
#                                  mutate_all(funs(as.numeric(.))) %>%
#                                  as.matrix(),
#                                x_predict_moderate = current.average.patient %>%
#                                  select(
#                                    all_of(variables_tau)
#                                  ) %>%
#                                  mutate_all(funs(as.numeric(.))) %>%
#                                  as.matrix(),
#                                pi_pred = rep(0.5, nrow(current.average.patient)),
#                                z_pred = rep(rep(c(0,1), each = length(breaks_hba1c)), 2),
#                                save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_prehba1c <- cbind(
#   values = round(current.average.patient$prehba1c),
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(breaks_hba1c)), 2),
#   mean = apply(predictions.prehba1c$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.prehba1c$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.prehba1c$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = round(breaks_hba1c)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Prehba1c") +
#   ylim(-8, 12)
# 
# # preegfr
# # breaks_egfr
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(breaks_egfr))) %>%
#   mutate(preegfr = rep(breaks_egfr, 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(breaks_egfr))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)))
# 
# predictions.preegfr <- predict(object = bcf_model,
#                              x_predict_control = current.average.patient %>%
#                                select(
#                                  all_of(variables_mu)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              x_predict_moderate = current.average.patient %>%
#                                select(
#                                  all_of(variables_tau)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              pi_pred = rep(0.5, nrow(current.average.patient)),
#                              z_pred = rep(rep(c(0,1), each = length(breaks_egfr)), 2),
#                              save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_preegfr <- cbind(
#   values = round(current.average.patient$preegfr),
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(breaks_egfr)), 2),
#   mean = apply(predictions.preegfr$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.preegfr$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.preegfr$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = round(breaks_egfr)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Preegfr") +
#   ylim(-8, 12)
# 
# 
# 
# # agetx
# # breaks_agetx
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(breaks_agetx))) %>%
#   mutate(agetx = rep(breaks_agetx, 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(breaks_agetx))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)))
# 
# predictions.agetx <- predict(object = bcf_model,
#                              x_predict_control = current.average.patient %>%
#                                select(
#                                  all_of(variables_mu)
#                                  ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              x_predict_moderate = current.average.patient %>%
#                                select(
#                                  all_of(variables_tau)
#                                  ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              pi_pred = rep(0.5, nrow(current.average.patient)),
#                              z_pred = rep(rep(c(0,1), each = length(breaks_agetx)), 2),
#                              save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_agetx <- cbind(
#   values = round(current.average.patient$agetx),
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(breaks_agetx)), 2),
#   mean = apply(predictions.agetx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.agetx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.agetx$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
#   ) %>%
#   as.data.frame() %>%
#   mutate(values = factor(values, levels = round(breaks_agetx)),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Agetx") +
#   ylim(-8, 12)
# 
# # hba1cmonth
# # breaks_hba1cmonth
# 
# current.average.patient <- average.patient %>%
#   slice(rep(1, 2*length(breaks_hba1cmonth))) %>%
#   mutate(hba1cmonth = rep(breaks_hba1cmonth, 2)) %>%
#   slice(rep(1:n(), 2)) %>%
#   cbind(sex = rep(c(2,1), each = 2*length(breaks_hba1cmonth))) %>%
#   mutate(sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex)))
# 
# predictions.hba1cmonth <- predict(object = bcf_model,
#                              x_predict_control = current.average.patient %>%
#                                select(
#                                  all_of(variables_mu)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              x_predict_moderate = current.average.patient %>%
#                                select(
#                                  all_of(variables_tau)
#                                ) %>%
#                                mutate_all(funs(as.numeric(.))) %>%
#                                as.matrix(),
#                              pi_pred = rep(0.5, nrow(current.average.patient)),
#                              z_pred = rep(rep(c(0,1), each = length(breaks_hba1cmonth)), 2),
#                              save_tree_directory = paste0(output_path, "/response_model/trees"))
# 
# plot_hba1cmonth <- cbind(
#   values = current.average.patient$hba1cmonth,
#   drugclass = rep(rep(c("GLP1", "SGLT2"), each = length(breaks_hba1cmonth)), 2),
#   mean = apply(predictions.hba1cmonth$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.5))),
#   uci = apply(predictions.hba1cmonth$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.95))),
#   lci = apply(predictions.hba1cmonth$tau, MARGIN = 2, function(x) quantile(c(x), probs = c(0.05))),
#   sex = current.average.patient$sex
# ) %>%
#   as.data.frame() %>%
#   mutate(values = as.numeric(values),
#          drugclass = factor(drugclass, levels = levels(hba1c.train$drugclass)),
#          mean = as.numeric(mean),
#          uci = as.numeric(uci),
#          lci = as.numeric(lci),
#          sex = factor(sex, levels = 1:length(levels(hba1c.train$sex)), labels = levels(hba1c.train$sex))
#   ) %>%
#   ggplot(aes(x = values, y = mean)) +
#   geom_hline(aes(yintercept = 0), colour = "red") +
#   geom_point() +
#   geom_errorbar(aes(x = values, y = mean, ymax = uci, ymin = lci)) +
#   facet_wrap(~sex) +
#   xlab("Hba1cmonth") +
#   ylim(-8, 12)
# 
# 
# # Combine all
# plot_predictions <- patchwork::wrap_plots(
#   list(
#     plot_drugline,
#     plot_ncurrtx,
#     plot_prepad,
#     plot_preihd,
#     plot_preneuropathy,
#     plot_preretinopathy,
#     plot_preaf,
#     plot_prehba1c,
#     plot_preegfr,
#     plot_agetx,
#     plot_hba1cmonth
#   )) +
#   patchwork::plot_annotation(
#     title = "Prediction of treatment effect for the average patient"
#   )

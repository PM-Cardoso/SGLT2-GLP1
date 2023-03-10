####################
## Description:
##  - In this file we create the plots used in the supplementary material of the paper
####################


## Load libraries
library(tidyverse)
library(patchwork)
library(scales)
library(forestplot)
library(rms)

## Load functions required
source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

## make directory for outputs
dir.create("Plots/Paper")

## make directory for outputs
dir.create("Plots/Paper/Sup_Mat")


#### Supplementary manuscript plots



#:--------------------------------------------------------------------------------
### Variable selection from the Sparse BCF model


## Variable selection of sparcebcf for the moderator effect
variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_tau_original_chain_1.rds")

plot_variables_tau <- variables_tau %>%
  rename("Current age" = "agetx",
         "Sex" = "sex",
         "Duration of diabetes" = "t2dmduration",
         "Number of glucose-lowering drug classes ever prescribed" = "drugline",
         "Number of other current glucose-lowering drugs" = "ncurrtx",
         "Month of Hba1C measure" = "hba1cmonth",
         "HbA1c" = "prehba1c",
         "BMI" = "prebmi",
         "eGFR" = "preegfr",
         "Albumin" = "prealbuminblood",
         "Alanine transaminase" = "prealt",
         "Bilirubin" = "prebilirubin",
         "HDL" = "prehdl",
         "Mean arterial blood pressure" = "premap",
         "Total cholesterol" = "pretotalcholesterol",
         "Angina" = "preangina",
         "Chronic liver disease" = "precld",
         "Nephropathy" = "prediabeticnephropathy",
         "Heart failure" = "preheartfailure",
         "Hypertension" = "prehypertension",
         "Ischaemic heart disease" = "preihd",
         "Myocardial infarction" = "premyocardialinfarction",
         "Neuropathy" = "preneuropathy",
         "Peripheral arterial disease" = "prepad",
         "Retinopathy" = "preretinopathy",
         "Atherosclerotic cardiovascular disease" = "prerevasc",
         "Stroke" = "prestroke",
         "Transient ischaemic attack" = "pretia",
         "Atrial fibrillation" = "preaf") %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below")) %>%
  ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
  geom_vline(aes(xintercept = 0.023), colour = "red") +
  geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
  geom_point(size = 2) +
  ggtitle("Moderator model") +
  xlab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
  theme_light() +
  theme(axis.text.y = element_text(angle = 30, face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none")


## Variable selection of sparcebcf for the prognostic effect
variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/assessment/variables_mu_original_chain_1.rds")

plot_variables_mu <- variables_mu %>%
  rename("Current age" = "agetx",
         "Sex" = "sex",
         "Duration of diabetes" = "t2dmduration",
         "Number of glucose-lowering drug classes ever prescribed" = "drugline",
         "Number of other current glucose-lowering drugs" = "ncurrtx",
         "Month of Hba1C measure" = "hba1cmonth",
         "HbA1c" = "prehba1c",
         "BMI" = "prebmi",
         "eGFR" = "preegfr",
         "Albumin" = "prealbuminblood",
         "Alanine transaminase" = "prealt",
         "Bilirubin" = "prebilirubin",
         "HDL" = "prehdl",
         "Mean arterial blood pressure" = "premap",
         "Total cholesterol" = "pretotalcholesterol",
         "Angina" = "preangina",
         "Chronic liver disease" = "precld",
         "Nephropathy" = "prediabeticnephropathy",
         "Heart failure" = "preheartfailure",
         "Hypertension" = "prehypertension",
         "Ischaemic heart disease" = "preihd",
         "Myocardial infarction" = "premyocardialinfarction",
         "Neuropathy" = "preneuropathy",
         "Peripheral arterial disease" = "prepad",
         "Retinopathy" = "preretinopathy",
         "Atherosclerotic cardiovascular disease" = "prerevasc",
         "Stroke" = "prestroke",
         "Transient ischaemic attack" = "pretia",
         "Atrial fibrillation" = "preaf",
         "Propensity score" = "propensity score") %>%
  as.data.frame() %>%
  gather(key, value) %>%
  arrange(desc(value)) %>%
  mutate(key = factor(key),
         colour = ifelse(value > 0.023, "Above", "Below")) %>%
  ggplot(aes(y = forcats::fct_reorder(key, value), x = value, colour = colour)) +
  geom_vline(aes(xintercept = 0.023), colour = "red") +
  geom_segment(aes(x = 0, xend = value, yend = forcats::fct_reorder(key, value)), linetype = "dashed") +
  geom_point(size = 2) +
  ggtitle("Control model") +
  xlab("Posterior splitting probabilities") +
  scale_colour_manual(values = c("Above" = "black", "Below" = "grey")) +
  theme_light() +
  theme(axis.text.y = element_text(angle = 30, face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none")

# combining plots
plot_sup_1 <- patchwork::wrap_plots(list(plot_variables_mu, plot_variables_tau), ncol = 2) +
  patchwork::plot_annotation(tag_levels = c('A'),
                             title = "BART model variable selection") # title of full plot

# pdf of plot
pdf(width = 14, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_sup_1.pdf")
plot_sup_1
dev.off()


#:--------------------------------------------------------------------------------
### Sex difference in trials

# plot_sup_2


#:--------------------------------------------------------------------------------
### Figure of treatment effects validation for Ethnicity

interval_breaks <- c(-5, -3, 0, 3, 5)

variables_mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_mu.rds")

variables_tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/variables_tau.rds")

# Create dataset with all the variables needed
ethnicity.dataset <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated")) %>%
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "effects"))) %>%
  drop_na() %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other", "Mixed"), labels = c("White", "South Asian", "Black", "Other/Mixed", "Other/Mixed")))

# Create variable corresponding to the clinical subgroups
ethnicity.dataset.grouped <- group_values(data = ethnicity.dataset,
                                          variable = "effects",
                                          breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

# Subsetting dataset ethnicity - White
predicted_observed_white <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "White")

# Subsetting dataset ethnicity - South Asian
predicted_observed_asian <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "South Asian")

# Subsetting dataset ethnicity - Black
predicted_observed_black <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Black")

# Subsetting dataset ethnicity - Other/Mixed
predicted_observed_mixed <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Other/Mixed")



# White
ATE_adjust_white <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# South Asian
ATE_adjust_asian <- calc_ATE(data = predicted_observed_asian%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Black
ATE_adjust_black <- calc_ATE(data = predicted_observed_black%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Mixed / Other
ATE_adjust_mixed <- calc_ATE(data = predicted_observed_mixed%>%mutate(intervals=as.numeric(intervals)), validation_type = "Adjust",
                             variable = "posthba1cfinal", quantile_var = "intervals",
                             order = "largest", breakdown = unique(c(variables_mu, variables_tau)))


### Overall analysis (1 group)
# White
ATE_adjust_white_overall <- calc_ATE(data = predicted_observed_white%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# South Asian
ATE_adjust_asian_overall <- calc_ATE(data = predicted_observed_asian%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Black
ATE_adjust_black_overall <- calc_ATE(data = predicted_observed_black%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Mixed / Other
ATE_adjust_mixed_overall <- calc_ATE(data = predicted_observed_mixed%>%mutate(intervals=as.numeric(1)), validation_type = "Adjust",
                                     variable = "posthba1cfinal", quantile_var = "intervals",
                                     order = "largest", breakdown = unique(c(variables_mu, variables_tau)))

# Set up axis limits for the forest plots
hba1c_strata_axis_min <- plyr::round_any(floor(min(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_asian_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_black_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  ATE_adjust_mixed_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% min()
))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max()
))), 2, f = ceiling)


# plot of forest plot
plot_sup_3 <- rbind(
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "South Asian", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Black", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Other/Mixed", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "White"),
  cbind(ATE_adjust_asian[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "South Asian"),
  cbind(ATE_adjust_black[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "Black"),
  cbind(ATE_adjust_mixed[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3), group = "Other/Mixed"),
  cbind(obs = 50, lci = 50, uci = 50, group = "White", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "South Asian", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Black", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(obs = 50, lci = 50, uci = 50, group = "Other/Mixed", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(ATE_adjust_white[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "White"),
  cbind(ATE_adjust_asian[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "South Asian"),
  cbind(ATE_adjust_black[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "Black"),
  cbind(ATE_adjust_mixed[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6), group = "Other/Mixed"),
  cbind(ATE_adjust_white_overall[["effects"]] %>% select(obs, lci, uci), group = "White", intervals = "Average treatment effect"),
  cbind(ATE_adjust_asian_overall[["effects"]] %>% select(obs, lci, uci), group = "South Asian", intervals = "Average treatment effect"),
  cbind(ATE_adjust_black_overall[["effects"]] %>% select(obs, lci, uci), group = "Black", intervals = "Average treatment effect"),
  cbind(ATE_adjust_mixed_overall[["effects"]] %>% select(obs, lci, uci), group = "Other/Mixed", intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, ">5 mmol/mol",
                            ifelse(intervals == 2, "3-5 mmol/mol", 
                                   ifelse(intervals == 3, "0-3 mmol/mol",
                                          ifelse(intervals == 4, "0-3 mmol/mol", 
                                                 ifelse(intervals == 5, "3-5 mmol/mol",
                                                        ifelse(intervals == 6, ">5 mmol/mol", intervals)))))),
         group = factor(group)) %>%
  rename("lower" = "lci", "upper" = "uci", "mean" = "obs", "labeltext" = "intervals") %>%
  mutate(n_white = c(rep(NA, 4), rep(format(ATE_adjust_white[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_white[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_white[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_asian = c(rep(NA, 4), rep(format(ATE_adjust_asian[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_asian[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_asian[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_black = c(rep(NA, 4), rep(format(ATE_adjust_black[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_black[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_black[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         n_mixed = c(rep(NA, 4), rep(format(ATE_adjust_mixed[["effects"]]$N[1:3],big.mark=",",scientific=FALSE), 4),
                     rep(NA, 4), rep(format(ATE_adjust_mixed[["effects"]]$N[4:6],big.mark=",",scientific=FALSE), 4),
                     rep(format(sum(ATE_adjust_mixed[["effects"]]$N),big.mark=",",scientific=FALSE), 4)),
         mean = as.numeric(mean)) %>%
  group_by(group) %>%
  forestplot(labeltext = c("labeltext", "n_white", "n_asian", "n_black", "n_mixed"),
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(labeltext = paste0("Overall population (n=", format(ethnicity.dataset.grouped%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                n_white = "White",
                n_asian = "South Asian",
                n_black = "Black",
                n_mixed = "Other/Mixed") %>%
  fp_set_style(box = c("#E64B35B2","#4DBBD5B2", "#00A087B2", "#3C5488B2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


pdf(width = 15, height = 6, "Plots/Paper/Sup_Mat/11.08.plot_sup_3.pdf")
plot_sup_3
dev.off()


#:--------------------------------------------------------------------------------
### Figure of risk of developing CKD (40% drop of Stage 5)

patient_prop_scores_qrisk <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/patient_prop_scores_qrisk.rds")

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_effects.rds"), by = c("patid", "pated"))

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.no_co.dataset,
  method = "nearest",
  distance = group.no_co.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

group.no_co.dataset.matched <- group.no_co.dataset %>%
  slice(which(matching_no_co$weights == 1))

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall.rds")

# predictions for the CKD - stage 3/4/5 outcomes in the population with no CVD/HF/CKD
predictions_no_co_egfr40_or_ckd5_stan_adjusted_full <- readRDS("Samples/SGLT2-GLP1/Aurum/additional_outcomes/predictions_no_co_egfr40_or_ckd5_stan_adjusted_full.rds")

# set axis limits
axis_min <- exp(min(c(predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                  predictions_no_co_egfr40_or_ckd5_stan_adjusted_full%>%unlist()%>%as.numeric())))

axis_max <- exp(max(c(predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall%>%select(-intervals)%>%unlist()%>%as.numeric(),
                      predictions_no_co_egfr40_or_ckd5_stan_adjusted_full%>%unlist()%>%as.numeric())))

# plot forest plot
## PSM
plot_psm <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "PSM (1:1)",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")


## PSM + adjusted
plot_psm_adusted <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "PSM (1:1) and adjustment",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset.matched),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")

## Adjusted
plot_adusted <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_overall %>% slice(4:6),
  predictions_no_co_egfr40_or_ckd5_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(),big.mark=",",scientific=FALSE), ", events = ", format(group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_egfr40_or_ckd5_censvar==1)%>%nrow(),big.mark=",",scientific=FALSE), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjustment",
             xlog = TRUE,
             clip = c(axis_min, axis_max),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazards ratio (95% CI, log scale)") %>%
  fp_add_header(paste0("Overall population (n=", format(nrow(group.no_co.dataset),big.mark=",",scientific=FALSE), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = "black")



pdf(width = 7, height = 12, "Plots/Paper/Sup_Mat/11.08.plot_sup_4.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3,
                                           ncol = 1, heights = unit(c(5, 5, 5), "null"))))

# first plot
pushViewport(viewport(layout.pos.row = 1,
                      layout.pos.col = 1))
plot_psm
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_psm_adusted
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_adusted
upViewport()

dev.off()


#:--------------------------------------------------------------------------------
### Compare predictions from BCF SGLT2 vs John et al. SGLT2

# read in predictions from BCF model
patient_predicted_outcomes <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_predicted_outcomes.rds")

# load the object for the Linear model. This model object is the object used
#   in the SGLT2vsDPP4 model by John Dennis et al.
load("m1_hba1cmodel_SGLT2_DPP4.Rdata")

# Collect the full cohort from CPRD Aurum that can be used to fit the SGLT2vsDPP4 linear model
full.cohort.updated <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  select(patid, pated, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth) %>%
  rename("prealtlog"="prealt",
         "prehba1cmmol"="prehba1c",
         "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog),
         ncurrtx = factor(ncurrtx, levels = c("1","2","3","4","5+"), labels = c("0","1","2","3","3")),
         drugline = factor(drugline, levels = c("2","3","4","5+"), labels = c("2","3","4","5")),
         hba1cmonth = ifelse(is.na(hba1cmonth), 12, hba1cmonth)) %>%
  drop_na() %>%
  # apply rules from the model
  filter(prehba1cmmol < 120) %>%
  filter(prehba1cmmol >= 53) %>%
  filter(egfr_ckdepi > 45)


# predict for SGLT2
dataset.sglt2 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("SGLT2", levels = c("SGLT2", "DPP4")))
predictions.sglt2 <- predict(m1, dataset.sglt2)



## Compare SGLT2 BCF vs Linear regression
interim.dataset <- full.cohort.updated %>%
  select(patid, pated) %>%
  cbind(pred.SGLT2.lm = predictions.sglt2) %>%
  left_join(patient_predicted_outcomes %>%
              select(patid, pated, pred.SGLT2) %>%
              rename("pred.SGLT2.bcf" = "pred.SGLT2"), by = c("patid", "pated")) %>%
  select(-patid, -pated)

# Plot the comparison of SGLT2 predictions from the linear SGLT2vsDPP4 model vs SGLT2vsGLP1 model
plot_sup_5.1 <- interim.dataset %>%
  ggplot(aes(y = pred.SGLT2.bcf, x = pred.SGLT2.lm)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red", linetype = "dashed") +
  stat_smooth() +
  ylab("SGLT2 predictions with BCF") +
  xlab("SGLT2 predictions with Linear Regression") +
  ylim(min(interim.dataset), max(interim.dataset)) +
  xlim(min(interim.dataset), max(interim.dataset)) +
  ggtitle(paste0("Comparison of SGLT2i predictions (mmol/mol) (n=", format(interim.dataset %>% nrow(),big.mark=",",scientific=FALSE), ")")) +
  theme_bw()

# Plot histogram of difference between both predictions
plot_sup_5.2 <- interim.dataset %>%
  mutate(diff = as.numeric(pred.SGLT2.bcf - pred.SGLT2.lm)) %>%
  ggplot(aes(x = diff)) +
  geom_density() +
  xlab("HbA1c difference (mmol/mol)") +
  ggtitle("Difference between HbA1c predictions (mmol/mol)") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


## PDF with the plot for the comparison  
pdf(width = 7, height = 7, "Plots/Paper/Sup_Mat/11.08.plot_sup_5.pdf")

plot_sup_5 <- patchwork::wrap_plots(list(plot_sup_5.1, plot_sup_5.2), ncol = 2, nrow = 1) +
  plot_annotation(tag_levels = "A")

dev.off()


#:--------------------------------------------------------------------------------
### Optimal therapy from both models combined

# predict for DPP4
dataset.dpp4 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("DPP4", levels = c("SGLT2", "DPP4")))
predictions.dpp4 <- predict(m1, dataset.dpp4)

interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "DPP4i", NA))),
         best_drug = factor(best_drug))

# Plot bar plot of the numbers of best drugs
plot_sup_6 <- interim.dataset %>%
  select(best_drug) %>%
  table() %>%
  as.data.frame() %>%
  rename("best_drug" = ".") %>%
  ggplot(aes(x = best_drug, y = Freq, fill = best_drug)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = format(Freq,big.mark=",",scientific=FALSE)), vjust = -0.5) +
  xlab("Optimal predicted therapy") +
  ylab("Number of patients") +
  ggtitle(paste0("Predicted optimal therapy (n=", format(interim.dataset%>%nrow(),big.mark=",",scientific=FALSE), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(label=comma)

## PDF with the plot for the best therapy
pdf(width = 7, height = 7, "Plots/Paper/Sup_Mat/11.08.plot_sup_6.pdf")
plot_sup_6
dev.off()


#:--------------------------------------------------------------------------------
### ROC AUC curve for propensity score model

require(ROCR)

# load PS model
bart_ps_model_final <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/bart_ps_model_final.rds")

# load training dataset
ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train")

# combine probabilities with real values
predictions.train <- bart_ps_model_final$p_hat_train

df.train <- as.data.frame(cbind(predictions.train, ps.model.train["drugclass"] %>%
                            mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
  set_names(c("predictions", "labels"))


performance.train <- pROC::roc(df.train$labels, df.train$predictions)

plot_sup_7.1 <- ggroc(performance.train, legacy.axes = "TRUE") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = paste0("Development cohort (n=", format(nrow(ps.model.train),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Sensitivity vs 1-Specificity"
  )

precision_recall.train <- performance.train %>%
  coords(ret = "all", transpose = FALSE) %>%
  select(precision, recall)

plot_sup_7.2 <- precision_recall.train %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = 0.5), lty = "dashed", colour = "grey") +
  geom_path() +
  theme_bw() +
  ylim(0, 1) +
  labs(
    x = "Recall",
    y = "Precision",
    title = paste0("Development cohort (n=", format(nrow(ps.model.train),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Precision vs Recall"
  )



# load PS model predictions for testing data
prop_score_testing_data <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/prop_score_testing_data.rds")

# load testing dataset
ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test")

# combined probabilities with real values
predictions.test <- prop_score_testing_data

df.test <- as.data.frame(cbind(predictions.test, ps.model.test["drugclass"] %>%
                            mutate(drugclass = ifelse(drugclass == "SGLT2", 0, 1)))) %>%
  set_names(c("predictions", "labels"))


performance.test <- pROC::roc(df.test$labels, df.test$predictions)

plot_sup_7.3 <- ggroc(performance.test, legacy.axes = "TRUE") +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dashed", colour = "red") +
  theme_bw() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = paste0("Validation cohort (n=", format(nrow(ps.model.test),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Sensitivity vs 1-Specificity"
  )

precision_recall.test <- performance.test %>%
  coords(ret = "all", transpose = FALSE) %>%
  select(precision, recall)

plot_sup_7.4 <- precision_recall.test %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = 0.5), lty = "dashed", colour = "grey") +
  geom_path() +
  theme_bw() +
  ylim(0, 1) +
  labs(
    x = "Recall",
    y = "Precision",
    title = paste0("Validation cohort (n=", format(nrow(ps.model.test),big.mark=",",scientific=FALSE), ")"),
    subtitle = "Precision vs Recall"
  )


plot_sup_7 <- patchwork::wrap_plots(list(plot_sup_7.1, plot_sup_7.2, plot_sup_7.3, plot_sup_7.4), ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))


## PDF with the plot for the best therapy
pdf(width = 8, height = 8, "Plots/Paper/Sup_Mat/11.08.plot_sup_7.pdf")
plot_sup_7
dev.off()



#:--------------------------------------------------------------------------------
### Variable selection for PS model

# load proportions of variable selection
vs_bart_ps_model <- readRDS("Samples/SGLT2-GLP1/Aurum/ps_model/vs_bart_ps_model.rds")

# dataset with values necessary for the plot
df <- t(apply(vs_bart_ps_model$permute_mat, 2, range)) %>%
  as.data.frame() %>%
  rename("min" = "V1",
         "max" = "V2") %>%
  cbind(labels = factor(names(vs_bart_ps_model$var_true_props_avg), 
                        levels = names(vs_bart_ps_model$var_true_props_av), 
                        labels = c("BMI", "Year of drug start", 
                                   "eGFR", "HbA1c", 
                                   "Number drugs ever prescribed - 5+", "Current age",
                                   "Duration of diabetes", "Number drugs ever prescribed - 2", 
                                   "Number drugs ever prescribed - 3", "Number other current drugs - 3",
                                   "Number other current drugs - 5+", "Heart failure - No",
                                   "Number other current drugs - 4", "Number drugs ever prescribed - 4", 
                                   "Heart failure - Yes", "Ethnicity - White",
                                   "Hospitalisation - No", "Chronic liver disease - Yes",
                                   "Sex - Female", "Number other current drugs - 2",
                                   "Hospitalisation - Yes", "Diabetic Nephropathy - Yes",
                                   "Chronic liver disease - No", "Number other current drugs - 1",
                                   "Retinopathy - No", "Sex - Male",
                                   "Retinopathy - Yes", "Neuropathy - No",
                                   "Ethnicity - South Asian", "Angina - Yes",
                                   "Atrial fibrillation - No", "Diabetic Nephropathy - No",
                                   "Neuropathy - Yes", "Peripheral arterial disease - No",
                                   "Angina - No", "Ethnicity - Mixed",
                                   "Transient ischaemic attack - No", "Ethnicity - Black",
                                   "Atrial fibrillation - Yes", "Ethnicity - Other",
                                   "Myocardial infarction - Yes", "Deprivation - 9",
                                   "Stroke - No", "Ischaemic heart disease - Yes",
                                   "Ischaemic heart disease - No", "Smoking status - Active",
                                   "Deprivation - 6", "Peripheral arterial disease - Yes",
                                   "Transient ischaemic attack - Yes", "Smoking status - Non-smoker",
                                   "Cardiac revascularisation - Yes", "Cardiac revascularisation - No",
                                   "Stroke - Yes", "Myocardial infarction - No",
                                   "Smoking status - Ex-smoker", "Hypertension - Yes",
                                   "Hypertension - No", "Deprivation - 8",
                                   "Deprivation - 1", "Deprivation - 5",
                                   "Deprivation - 3", "Deprivation - 4",
                                   "Deprivation - 7", "Deprivation - 10",
                                   "Deprivation - 2")),
        true.prop = vs_bart_ps_model$var_true_props_avg) %>%
  mutate(min = as.numeric(min),
         max = as.numeric(max),
         mean = (min+max)/2,
         max.values = ifelse(labels == "agetx", max, mean),
         colour.toggle = factor(ifelse(true.prop > max.values, 1, 0)),
         min = 0) %>%
  gather(label, plot.values, -labels, -max, -true.prop, -colour.toggle) %>%
  filter(labels != "Hospitalisation - No")  %>%
  filter(labels != "Hospitalisation - Yes") 

# plot the variable selection
plot_sup_8 <- df %>%
  ggplot(aes(x = labels, y = plot.values)) +
  geom_boxplot(width = 0, colour = "chartreuse4") +
  geom_point(aes(y = true.prop, shape = colour.toggle)) +
  scale_shape_manual(values = c(1, 16)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        axis.title.x = element_blank()) +
  labs(
    y = "Proportion included"
  )


## PDF with the plot for the best therapy
pdf(width = 11, height = 5, "Plots/Paper/Sup_Mat/11.08.plot_sup_8.pdf")
plot_sup_8
dev.off()



#:--------------------------------------------------------------------------------
### Comparing predictions of mu and tau for BCF prop / no prop

# comparing mu predictions
prop.score.comparison.mu <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/prop.score.comparison.mu.rds")

# comparing tau predictions
prop.score.comparison.tau <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/prop.score.comparison.rds")

# plot comparisons
plot_sup_9 <- patchwork::wrap_plots(
  
  # scatter plot of tau values
  prop.score.comparison.tau %>%
    as.data.frame() %>%
    ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red", lty = "dashed") +
    geom_point() +
    ggtitle("Predicted conditional average treatment effects (CATE)") +
    xlab("BCF model with propensity score") +
    ylab("BCF model without propensity score")
  
  ,
  
  # density plot of residuals for tau (BCF without prop - BCF with prop)
  prop.score.comparison.tau %>%
    as.data.frame() %>%
    mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
    ggplot() +
    geom_density(aes(x = effect.difference)) +
    labs(
      x = "Difference in predicted CATE (mmol/mol)",
      subtitle = "BCF without propensity score - BCF with propensity score"
    ) +
    theme(axis.title.y = element_blank())
  
  ,
  
  # scatter plot of mu values
  prop.score.comparison.mu %>%
    as.data.frame() %>%
    ggplot(aes(x = bcf_prop, y = bcf_no_prop)) +
    geom_abline(aes(intercept = 0, slope = 1), colour = "red", lty = "dashed") +
    geom_point() +
    ggtitle("Predicted HbA1c outcome (mmol/mol)") +
    xlab("BCF model with propensity score") +
    ylab("BCF model without propensity score")
  
  ,
  
  # density plot of residuals for mu (BCF without prop - BCF with prop)
  prop.score.comparison.mu %>%
    as.data.frame() %>%
    mutate(effect.difference = bcf_no_prop - bcf_prop) %>%
    ggplot() +
    geom_density(aes(x = effect.difference)) +
    labs(
      x = "Difference in predicted HbA1c outcome (mmol/mol)",
      subtitle = "BCF without propensity score - BCF with propensity score"
    ) +
    theme(axis.title.y = element_blank())
  
  , ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))


## PDF with the plot for the best therapy
pdf(width = 11, height = 11, "Plots/Paper/Sup_Mat/11.08.plot_sup_9.pdf")
plot_sup_9
dev.off()







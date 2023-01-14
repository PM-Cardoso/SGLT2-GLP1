####################
## Description:
##  - In this file we divide the population into ethnicity groups and validate each portion.
####################


## libraries
library(tidyverse)
library(forestplot)

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

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################


full.cohort <- set_up_data_sglt2_glp1(dataset.type="full.cohort")

## Read in propensity scores
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Read in treatment effects
treatment_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))



#----------------------------------------------------------

ethnicity.dataset <- full.cohort %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "prop.score", "effects"))) %>%
  drop_na() %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other", "Mixed"), labels = c("White", "South Asian", "Black", "Other/Mixed", "Other/Mixed")))

predicted_observed_white <- ethnicity.dataset %>%
  filter(ethnicity == "White") %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(hba1c_diff.q = ntile(hba1c_diff, 10))

predicted_observed_asian <- ethnicity.dataset %>%
  filter(ethnicity == "South Asian") %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(hba1c_diff.q = ntile(hba1c_diff, 5))

edicted_observed_black <- ethnicity.dataset %>%
  filter(ethnicity == "Black") %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(hba1c_diff.q = ntile(hba1c_diff, 5))

predicted_observed_mixed <- ethnicity.dataset %>%
  filter(ethnicity == "Other/Mixed") %>%
  rename("hba1c_diff" = "effects") %>%
  mutate(hba1c_diff.q = ntile(hba1c_diff, 5))

#### ---------------------- Propensity score matching
ATE_psm_1_1_white <- calc_ATE_validation_prop_matching(predicted_observed_white, "posthba1cfinal", predicted_observed_white%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", adjust = FALSE, order = "largest")

ATE_psm_1_1_asian <- calc_ATE_validation_prop_matching(predicted_observed_asian, "posthba1cfinal", predicted_observed_asian%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", adjust = FALSE, order = "largest")

ATE_psm_1_1_black <- calc_ATE_validation_prop_matching(edicted_observed_black, "posthba1cfinal", edicted_observed_black%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", adjust = FALSE, order = "largest")

ATE_psm_1_1_mixed <- calc_ATE_validation_prop_matching(predicted_observed_mixed, "posthba1cfinal", predicted_observed_mixed%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", adjust = FALSE, order = "largest")


#### ---------------------- Propensity score matching + adjusted
ATE_psm_1_1_adjust_white <- calc_ATE_validation_prop_matching(predicted_observed_white, "posthba1cfinal", predicted_observed_white%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_asian <- calc_ATE_validation_prop_matching(predicted_observed_asian, "posthba1cfinal", predicted_observed_asian%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_black <- calc_ATE_validation_prop_matching(edicted_observed_black, "posthba1cfinal", edicted_observed_black%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_mixed <- calc_ATE_validation_prop_matching(predicted_observed_mixed, "posthba1cfinal", predicted_observed_mixed%>%select(prop.score)%>%unlist(), quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")


#### ---------------------- Adjusted
ATE_adjust_white <- calc_ATE_validation_adjust(predicted_observed_white, "posthba1cfinal", quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_asian <- calc_ATE_validation_adjust(predicted_observed_asian, "posthba1cfinal", quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_black <- calc_ATE_validation_adjust(edicted_observed_black, "posthba1cfinal", quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_mixed <- calc_ATE_validation_adjust(predicted_observed_mixed, "posthba1cfinal", quantile_var = "hba1c_diff.q", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)




#:----------------------------------
### Plots
plot_ATE_psm_1_1_white <- ATE_plot(ATE_psm_1_1_white[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_white[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_asian <- ATE_plot(ATE_psm_1_1_asian[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_asian[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_black <- ATE_plot(ATE_psm_1_1_black[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_black[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_mixed <- ATE_plot(ATE_psm_1_1_mixed[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_mixed[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_adjust_white <- ATE_plot(ATE_psm_1_1_adjust_white[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_adjust_white[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_adjust_asian <- ATE_plot(ATE_psm_1_1_adjust_asian[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_adjust_asian[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_adjust_black <- ATE_plot(ATE_psm_1_1_adjust_black[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_adjust_black[["effects"]]$n_drug1)*2, ")"))
plot_ATE_psm_1_1_adjust_mixed <- ATE_plot(ATE_psm_1_1_adjust_mixed[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_psm_1_1_adjust_mixed[["effects"]]$n_drug1)*2, ")"))
plot_ATE_adjust_white <- ATE_plot(ATE_adjust_white[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_adjust_white[["effects"]]$N), ")"))
plot_ATE_adjust_asian <- ATE_plot(ATE_adjust_asian[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_adjust_asian[["effects"]]$N), ")"))
plot_ATE_adjust_black <- ATE_plot(ATE_adjust_black[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_adjust_black[["effects"]]$N), ")"))
plot_ATE_adjust_mixed <- ATE_plot(ATE_adjust_mixed[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci") +
  ggtitle(paste0("(n = ", sum(ATE_adjust_mixed[["effects"]]$N), ")"))



#:----------------------------------

### Combine plots
pdf(width = 10, height = 10, "Plots/11.10.plot_1.pdf")
# Propensity score matching
patchwork::wrap_plots(list(plot_ATE_psm_1_1_white, 
                           plot_ATE_psm_1_1_asian, 
                           plot_ATE_psm_1_1_black, 
                           plot_ATE_psm_1_1_mixed), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Treatment effect validation - propensity score matching", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5)))

# Propensity score matching + adjusted
patchwork::wrap_plots(list(plot_ATE_psm_1_1_adjust_white, 
                           plot_ATE_psm_1_1_adjust_asian, 
                           plot_ATE_psm_1_1_adjust_black, 
                           plot_ATE_psm_1_1_adjust_mixed), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Treatment effect validation - propensity score matching + adjusted", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5)))

# Adjusted
patchwork::wrap_plots(list(plot_ATE_adjust_white, 
                           plot_ATE_adjust_asian, 
                           plot_ATE_adjust_black, 
                           plot_ATE_adjust_mixed), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A",
                             title = "Treatment effect validation - adjusted", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5)))

dev.off()





#:---------------------------------------------------------------------------------------
# Split by subgroups

interval_breaks <- c(-5, -3, 0, 3, 5)

ethnicity.dataset <- full.cohort %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "prop.score", "effects"))) %>%
  drop_na() %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other", "Mixed"), labels = c("White", "South Asian", "Black", "Other/Mixed", "Other/Mixed")))


ethnicity.dataset.grouped <- group_values(data = ethnicity.dataset,
                                          variable = "effects",
                                          breaks = interval_breaks) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

predicted_observed_white <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "White")

predicted_observed_asian <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "South Asian")

predicted_observed_black <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Black")

predicted_observed_mixed <- ethnicity.dataset.grouped %>%
  filter(ethnicity == "Other/Mixed")


# #### ---------------------- Propensity score matching
# ATE_psm_1_1_white <- calc_ATE_validation_prop_matching(predicted_observed_white%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_white%>%select(prop.score)%>%unlist(), quantile_var = "intervals", adjust = FALSE, order = "largest")
# 
# ATE_psm_1_1_asian <- calc_ATE_validation_prop_matching(predicted_observed_asian%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_asian%>%select(prop.score)%>%unlist(), quantile_var = "intervals", adjust = FALSE, order = "largest")
# 
# ATE_psm_1_1_black <- calc_ATE_validation_prop_matching(edicted_observed_black%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", edicted_observed_black%>%select(prop.score)%>%unlist(), quantile_var = "intervals", adjust = FALSE, order = "largest")
# 
# ATE_psm_1_1_mixed <- calc_ATE_validation_prop_matching(predicted_observed_mixed%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_mixed%>%select(prop.score)%>%unlist(), quantile_var = "intervals", adjust = FALSE, order = "largest")
# 
# 
# #### ---------------------- Propensity score matching + adjusted
# ATE_psm_1_1_adjust_white <- calc_ATE_validation_prop_matching(predicted_observed_white%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_white%>%select(prop.score)%>%unlist(), quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")
# 
# ATE_psm_1_1_adjust_asian <- calc_ATE_validation_prop_matching(predicted_observed_asian%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_asian%>%select(prop.score)%>%unlist(), quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")
# 
# ATE_psm_1_1_adjust_black <- calc_ATE_validation_prop_matching(predicted_observed_black%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", edicted_observed_black%>%select(prop.score)%>%unlist(), quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")
# 
# ATE_psm_1_1_adjust_mixed <- calc_ATE_validation_prop_matching(predicted_observed_mixed%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", predicted_observed_mixed%>%select(prop.score)%>%unlist(), quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE, order = "largest")


#### ---------------------- Adjusted
ATE_adjust_white <- calc_ATE_validation_adjust(predicted_observed_white%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_asian <- calc_ATE_validation_adjust(predicted_observed_asian%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_black <- calc_ATE_validation_adjust(predicted_observed_black%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_mixed <- calc_ATE_validation_adjust(predicted_observed_mixed%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)


ATE_adjust_white_overall <- calc_ATE_validation_adjust(predicted_observed_white%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_asian_overall <- calc_ATE_validation_adjust(predicted_observed_asian%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_black_overall <- calc_ATE_validation_adjust(predicted_observed_black%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)

ATE_adjust_mixed_overall <- calc_ATE_validation_adjust(predicted_observed_mixed%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = unique(c(variables_mu, variables_tau)), adjust = TRUE)


hba1c_strata_axis_min <- plyr::round_any(floor(min(c(
  # ATE_psm_1_1_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_black[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
  # ATE_psm_1_1_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
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
  # ATE_psm_1_1_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_black[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  # ATE_psm_1_1_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_white[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_white_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_asian_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_black_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
  ATE_adjust_mixed_overall[["effects"]] %>% select(c("obs","lci","uci")) %>% max()
))), 2, f = ceiling)




plot_psm_1_1 <- rbind(
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
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ethnicity.dataset.grouped%>%filter(intervals==levels(ethnicity.dataset.grouped$intervals)[6])%>%nrow(),")"), intervals)))))),
         group = factor(group)) %>%
  rename("lower" = "lci", "upper" = "uci", "mean" = "obs", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             # fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(hba1c_strata_axis_min, hba1c_strata_axis_max),
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", ethnicity.dataset.grouped%>%nrow(), ")")) %>%
  fp_set_style(box = c("#E64B35B2","#4DBBD5B2", "#00A087B2", "#3C5488B2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



pdf(width = 9, height = 7, "Plots/11.10.plot_2.pdf")
plot_psm_1_1
dev.off()





#:---------------------------------------------------------------------------------------
# Treatment effect histogram split by Ethnicity 

ethnicity.dataset <- full.cohort %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau)), "prop.score", "effects"))) %>%
  drop_na() %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "South Asian", "Black", "Other", "Mixed"), labels = c("White", "South Asian", "Black", "Other/Mixed", "Other/Mixed"))) %>%
  rename("mean" = "effects")

plot_hist_white <- hist_plot(ethnicity.dataset %>% filter(ethnicity == "White"), title = "Ethnicity: white", xmin = -15, xmax = 20)

plot_hist_asian <- hist_plot(ethnicity.dataset %>% filter(ethnicity == "South Asian"), title = "Ethnicity: south asian", xmin = -15, xmax = 20)

plot_hist_black <- hist_plot(ethnicity.dataset %>% filter(ethnicity == "Black"), title = "Ethnicity: black", xmin = -15, xmax = 20)

plot_hist_mixed <- hist_plot(ethnicity.dataset %>% filter(ethnicity == "Other/Mixed"), title = "Ethnicity: other/mixed", xmin = -15, xmax = 20)


plot_3 <- patchwork::wrap_plots(list(plot_hist_white, plot_hist_asian, plot_hist_black, plot_hist_mixed), ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                             title = "Treatment effect heterogeneity", # title of full plot
                             theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot


pdf(width = 5, height = 12, "Plots/11.10.plot_3.pdf")
plot_3
dev.off()







# Strata population
qpdf::pdf_combine(input = c("Plots/11.10.plot_1.pdf",
                            "Plots/11.10.plot_2.pdf",
                            "Plots/11.10.plot_3.pdf"),
                  output = "Plots/11.10.ethnicity_validation.pdf")


file.remove(c("Plots/11.10.plot_1.pdf",
              "Plots/11.10.plot_2.pdf",
              "Plots/11.10.plot_3.pdf"))





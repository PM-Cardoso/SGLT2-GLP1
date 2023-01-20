####################
## Description:
##  - In this file we check:
##    - Weight change
##    - eGFR change
##    - Discontinuation
####################

###
## Need to try discontinuationw without rcs
###

# library(qpdf)
library(tidyverse)
library(margins)
library(forestplot)
library(grid)
library(rms)
library(cowplot)

# Analysis of outcomes
# library(rstanarm)
# library(brms)
# library(survminer)

## increase memery usage to 50gb of RAM
# options(java.parameters = "-Xmx100g")
# 
# library(bartMachine)

marginal_distribution <- function(x,var) {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 64, alpha = 0.4, position = "identity") +
    guides(fill = "none") +
    #theme_void() +
    theme(legend.title = element_blank(), panel.background = element_rect( fill = "white",color = "grey50")) +  
    #scale_x_continuous(breaks = seq(0,50,5)) +
    xlab(expression("SGLT2i HbA1c benefit")) +
    #theme(plot.margin = margin()) +
    theme(text = element_text(size = 14),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="grey50"))
  
}



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
dir.create(paste0(output_path, "/additional_outcomes"))

## make directory for outputs
dir.create(paste0(output_path, "/additional_outcomes/trees_weight"))

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

interval_breaks <- c(-5, -3, 0, 3, 5)

## Read in propensity scores
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Read in treatment effects
treatment_effects <- readRDS(paste0(output_path, "/response_model_bcf/patient_effects.rds"))

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

model_variables <- unique(c(variables_mu, variables_tau))[which(unique(c(variables_mu, variables_tau)) != "sex")]





#:--------------------------------------------------------
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

group.hba1c.dataset.male <- group.hba1c.dataset %>% filter(sex == "Male")
group.hba1c.dataset.female <- group.hba1c.dataset %>% filter(sex == "Female")

# Propensity score matching
ATE_psm_1_1_hba1c <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c_male <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c_female <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

ATE_psm_1_1_hba1c_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", adjust = FALSE, order = "largest")

# Propensity score matching + adjust
ATE_psm_1_1_adjust_hba1c <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c_male <- calc_ATE_validation_prop_matching(group.hba1c.dataset.male%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.male$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c_female <- calc_ATE_validation_prop_matching(group.hba1c.dataset.female%>%mutate(intervals = as.numeric(intervals)), "posthba1cfinal", group.hba1c.dataset.female$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

ATE_psm_1_1_adjust_hba1c_full <- calc_ATE_validation_prop_matching(group.hba1c.dataset%>%mutate(intervals = as.numeric(1)), "posthba1cfinal", group.hba1c.dataset$prop.score, quantile_var = "intervals", breakdown = model_variables, adjust = TRUE, order = "largest")

# Adjust
ATE_adjust_hba1c <- calc_ATE_validation_adjust(group.hba1c.dataset, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_male <- calc_ATE_validation_adjust(group.hba1c.dataset.male, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_female <- calc_ATE_validation_adjust(group.hba1c.dataset.female, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_adjust_hba1c_full <- calc_ATE_validation_adjust(group.hba1c.dataset %>% mutate(intervals = as.numeric(1)), "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)


# Axis
hba1c_overall_axis_min <- plyr::round_any(floor(min(c(ATE_psm_1_1_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_psm_1_1_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_psm_1_1_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                      ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_overall_axis_max <- plyr::round_any(ceiling(max(c(ATE_psm_1_1_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                        ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                        ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                        ATE_psm_1_1_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_psm_1_1_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                        ATE_adjust_hba1c[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)


hba1c_strata_axis_min <- plyr::round_any(floor(min(c(ATE_psm_1_1_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_psm_1_1_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_psm_1_1_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_psm_1_1_adjust_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_psm_1_1_adjust_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_adjust_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                     ATE_adjust_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_strata_axis_max <- plyr::round_any(ceiling(max(c(ATE_psm_1_1_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                       ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                       ATE_adjust_hba1c_full[["effects"]] %>% select(c("obs","lci","uci")) %>% min(),
                                                       ATE_psm_1_1_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_psm_1_1_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_psm_1_1_adjust_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_psm_1_1_adjust_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_adjust_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% max(),
                                                       ATE_adjust_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)



# Propensity score matching
plot_psm_1_1_overall_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "Propensity score matching",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", sum(ATE_psm_1_1_hba1c[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_psm_1_1_male_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Male (n=", sum(ATE_psm_1_1_hba1c_male[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_psm_1_1_female_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "Propensity score matching",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Female (n=", sum(ATE_psm_1_1_hba1c_female[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Propensity score matching + adjust
plot_psm_1_1_adjust_overall_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "Propensity score matching + adjust",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", sum(ATE_psm_1_1_adjust_hba1c[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_psm_1_1_adjust_male_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_male[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Male (n=", sum(ATE_psm_1_1_adjust_hba1c_male[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_psm_1_1_adjust_female_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_psm_1_1_adjust_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_psm_1_1_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == 1, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(1)%>%unlist()*2, ")"), 
                            ifelse(intervals == 2, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(2)%>%unlist()*2, ")"), 
                                   ifelse(intervals == 3, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(3)%>%unlist()*2, ")"), 
                                          ifelse(intervals == 4, paste0("0-3 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(4)%>%unlist()*2, ")"), 
                                                 ifelse(intervals == 5, paste0("3-5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(5)%>%unlist()*2, ")"),
                                                        ifelse(intervals == 6, paste0(">5 mmol/mol (n=", ATE_psm_1_1_adjust_hba1c_female[["effects"]]%>%select(n_drug1)%>%slice(6)%>%unlist()*2, ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "Propensity score matching + adjust",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Female (n=", sum(ATE_psm_1_1_adjust_hba1c_female[["effects"]]$n_drug1)*2, ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Adjust
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
             title = "Adjust",
             xticks = seq(hba1c_overall_axis_min, hba1c_overall_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Overall population (n=", group.hba1c.dataset%>%nrow(), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_adjust_male_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c_male[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[1])%>%nrow(), ")"), 
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[3])%>%nrow(), ")"), 
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[4])%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[5])%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hba1c.dataset.male%>%filter(intervals == levels(group.hba1c.dataset.male$intervals)[6])%>%nrow(), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Male (n=", group.hba1c.dataset.male%>%nrow(), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_adjust_female_hba1c <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", obs = NA, lci = NA, uci = NA),
  ATE_adjust_hba1c_female[["effects"]] %>% select(intervals, obs, lci, uci) %>% slice(4:6),
  ATE_adjust_hba1c_full[["effects"]] %>% select(intervals, obs, lci, uci) %>% mutate(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hba1c.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[1])%>%nrow(), ")"), 
                            ifelse(intervals == levels(group.hba1c.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.hba1c.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[3])%>%nrow(), ")"), 
                                          ifelse(intervals == levels(group.hba1c.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[4])%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.hba1c.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[5])%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hba1c.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hba1c.dataset.female%>%filter(intervals == levels(group.hba1c.dataset.female$intervals)[6])%>%nrow(), ")"), intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             title = "Adjust",
             xticks = seq(hba1c_strata_axis_min, hba1c_strata_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header(paste0("Female (n=", group.hba1c.dataset.female%>%nrow(), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



pdf(width = 7, height = 12, "Plots/hba1c_grouping_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HbA1c treatment effect", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_psm_1_1_overall_hba1c
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_psm_1_1_adjust_overall_hba1c
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_adjust_overall_hba1c
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/hba1c_grouping_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HbA1c treatment effect", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_psm_1_1_female_hba1c
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_psm_1_1_male_hba1c
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_psm_1_1_adjust_female_hba1c
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_psm_1_1_adjust_male_hba1c
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_adjust_female_hba1c
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_adjust_male_hba1c
upViewport()

dev.off()






#:--------------------------------------------------------
# Weight change

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

matching_weight <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ preweight + agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.weight.dataset,
  method = "nearest",
  distance = group.weight.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_weight, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_weight$match.matrix))

group.weight.dataset.matched <- group.weight.dataset %>%
  slice(which(matching_weight$weights == 1))


## Propensity score matching
if (class(try(
  
  # predictions for the weight adjusted model sex strata
  predictions_weight_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_psm_1_1 <- vector()
  
  formula <- "w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preweight, 3) + sex*factor(drugclass)"
  
  models_weight_psm_1_1 <- stan_glm(formula, 
                                    data = group.weight.dataset.matched,
                                    family = gaussian(link = "identity"),
                                    prior = normal(0, 2),
                                    prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_psm_1_1, paste0(output_path, "/additional_outcomes/models_weight_psm_1_1.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                                    filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                                    filter(sex == "Male") %>%
                                    mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                                  type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                    uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                    sex = "Male",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                                   filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                                   filter(sex == "Male") %>%
                                   mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                                 type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                    lci = quantile(male_glp1, probs = c(0.025)),
                                                                                    uci = quantile(male_glp1, probs = c(0.975)),
                                                                                    sex = "Male",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                                      filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                                      filter(sex == "Female") %>%
                                      mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                                    type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                    uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                    sex = "Female",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                                     filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                                     filter(sex == "Female") %>%
                                     mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                                   type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                    lci = quantile(female_glp1, probs = c(0.025)),
                                                                                    uci = quantile(female_glp1, probs = c(0.975)),
                                                                                    sex = "Female",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1 <- predictions_weight_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1.rds"))
  
  
}

if (class(try(
  
  # predictions for the weight adjusted model
  predictions_weight_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_psm_1_1.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_psm_1_1_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                               filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                               mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                             type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_overall <- rbind(predictions_weight_stan_psm_1_1_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                    lci = quantile(sglt2, probs = c(0.025)),
                                                                                                    uci = quantile(sglt2, probs = c(0.975)),
                                                                                                    drugclass = "SGLT2",
                                                                                                    intervals = levels(group.weight.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                              filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                              mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_overall <- rbind(predictions_weight_stan_psm_1_1_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                    lci = quantile(glp1, probs = c(0.025)),
                                                                                                    uci = quantile(glp1, probs = c(0.975)),
                                                                                                    drugclass = "GLP1",
                                                                                                    intervals = levels(group.weight.dataset.matched$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1_overall <- predictions_weight_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_overall.rds"))
  
  
}
if (class(try(
  
  # predictions for the weight adjusted model full
  predictions_weight_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_weight_stan_psm_1_1_full <- vector()
  
  models_weight_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_psm_1_1.rds"))
  
  sglt2 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                             mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                           type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_psm_1_1_full <- rbind(predictions_weight_stan_psm_1_1_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                            lci = quantile(sglt2, probs = c(0.025)),
                                                                                            uci = quantile(sglt2, probs = c(0.975)),
                                                                                            drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                            mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                          type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_psm_1_1_full <- rbind(predictions_weight_stan_psm_1_1_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                            lci = quantile(glp1, probs = c(0.025)),
                                                                                            uci = quantile(glp1, probs = c(0.975)),
                                                                                            drugclass = "GLP1"))
  
  predictions_weight_stan_psm_1_1_full <- predictions_weight_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_full.rds"))
  
  
}


## Propensity score matching + adjustment
if (class(try(
  
  # predictions for the weight adjusted model sex strata
  predictions_weight_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_psm_1_1_adjusted <- vector()
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)
  
  formula <- paste0("w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(preweight, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_weight_psm_1_1_adjusted <- stan_glm(formula, 
                                             data = group.weight.dataset.matched,
                                             family = gaussian(link = "identity"),
                                             prior = normal(0, 2),
                                             prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_weight_psm_1_1_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                            filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                          type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                                      uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                           filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                         type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(male_glp1, probs = c(0.025)),
                                                                                                      uci = quantile(male_glp1, probs = c(0.975)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                              filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                                      uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                             filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                           type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(female_glp1, probs = c(0.025)),
                                                                                                      uci = quantile(female_glp1, probs = c(0.975)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1_adjusted <- predictions_weight_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted.rds"))
  
  
  
}


if (class(try(
  
  # predictions for the weight adjusted model
  predictions_weight_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_psm_1_1_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_psm_1_1_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                       filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                       mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                     type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted_overall <- rbind(predictions_weight_stan_psm_1_1_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                      lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                      uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                      drugclass = "SGLT2",
                                                                                                                      intervals = levels(group.weight.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                      filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                      mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                    type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_psm_1_1_adjusted_overall <- rbind(predictions_weight_stan_psm_1_1_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                      lci = quantile(glp1, probs = c(0.025)),
                                                                                                                      uci = quantile(glp1, probs = c(0.975)),
                                                                                                                      drugclass = "GLP1",
                                                                                                                      intervals = levels(group.weight.dataset.matched$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1_adjusted_overall <- predictions_weight_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the weight adjusted model full
  predictions_weight_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_weight_stan_psm_1_1_adjusted_full <- vector()
  
  models_weight_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_psm_1_1_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                     mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                   type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_psm_1_1_adjusted_full <- rbind(predictions_weight_stan_psm_1_1_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                              lci = quantile(sglt2, probs = c(0.025)),
                                                                                                              uci = quantile(sglt2, probs = c(0.975)),
                                                                                                              drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                    mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                  type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_psm_1_1_adjusted_full <- rbind(predictions_weight_stan_psm_1_1_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                              lci = quantile(glp1, probs = c(0.025)),
                                                                                                              uci = quantile(glp1, probs = c(0.975)),
                                                                                                              drugclass = "GLP1"))
  
  predictions_weight_stan_psm_1_1_adjusted_full <- predictions_weight_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted_full.rds"))
  
  
}



## Adjustment
if (class(try(
  
  # predictions for the weight adjusted model sex strata
  predictions_weight_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_adjusted <- vector()
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)
  
  formula <- paste0("w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(preweight, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_weight_adjusted <- stan_glm(formula, 
                                     data = group.weight.dataset,
                                     family = gaussian(link = "identity"),
                                     prior = normal(0, 2),
                                     prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_adjusted, paste0(output_path, "/additional_outcomes/models_weight_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                            filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                          type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                      uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                      sex = "Male",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                           filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                         type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                      lci = quantile(male_glp1, probs = c(0.025)),
                                                                                      uci = quantile(male_glp1, probs = c(0.975)),
                                                                                      sex = "Male",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                              filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                      uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                      sex = "Female",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                             filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                           type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                      lci = quantile(female_glp1, probs = c(0.025)),
                                                                                      uci = quantile(female_glp1, probs = c(0.975)),
                                                                                      sex = "Female",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_adjusted <- predictions_weight_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted.rds"))
  
  
  
}

if (class(try(
  
  # predictions for the weight adjusted model
  predictions_weight_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.weight.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_weight_stan_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                       filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                       mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset.matched$drugclass))),
                     type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted_overall <- rbind(predictions_weight_stan_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(sglt2, probs = c(0.025)),
                                                                                                      uci = quantile(sglt2, probs = c(0.975)),
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    glp1 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                      filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                      mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset.matched$drugclass))),
                    type = "response") %>%
      rowMeans()
    
    predictions_weight_stan_adjusted_overall <- rbind(predictions_weight_stan_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                      lci = quantile(glp1, probs = c(0.025)),
                                                                                                      uci = quantile(glp1, probs = c(0.975)),
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_adjusted_overall <- predictions_weight_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the weight adjusted model full
  predictions_weight_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_weight_stan_adjusted_full <- vector()
  
  models_weight_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_weight_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                     mutate(drugclass = factor("SGLT2", levels = levels(group.weight.dataset$drugclass))),
                   type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_adjusted_full <- rbind(predictions_weight_stan_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                              lci = quantile(sglt2, probs = c(0.025)),
                                                                                              uci = quantile(sglt2, probs = c(0.975)),
                                                                                              drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_weight_adjusted, newdata = group.weight.dataset %>%
                    mutate(drugclass = factor("GLP1", levels = levels(group.weight.dataset$drugclass))),
                  type = "response") %>%
    rowMeans()
  
  predictions_weight_stan_adjusted_full <- rbind(predictions_weight_stan_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                              lci = quantile(glp1, probs = c(0.025)),
                                                                                              uci = quantile(glp1, probs = c(0.975)),
                                                                                              drugclass = "GLP1"))
  
  predictions_weight_stan_adjusted_full <- predictions_weight_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted_full.rds"))
  
  
}



#:--------------------
#:---- PLOTS
#:--------------------

## limits

weight_overall_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_weight_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_weight_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                       predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

weight_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                         predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                         predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                         predictions_weight_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                         predictions_weight_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                         predictions_weight_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


weight_strata_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_weight_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_weight_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                      predictions_weight_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

weight_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                        predictions_weight_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                        predictions_weight_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                        predictions_weight_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                        predictions_weight_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                        predictions_weight_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


#:---- PSM 1:1
plot_weight_psm_1_1_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(weight_overall_axis_min, weight_overall_axis_max),
             xticks = seq(weight_overall_axis_min, weight_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", group.weight.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_weight_psm_1_1_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Male (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_weight_psm_1_1_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Female (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



#:---- PSM 1:1 + adjusted
plot_weight_psm_1_1_adjusted_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(weight_overall_axis_min, weight_overall_axis_max),
             xticks = seq(weight_overall_axis_min, weight_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", group.weight.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_weight_psm_1_1_adjusted_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Male (n=", group.weight.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_weight_psm_1_1_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_weight_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Female (n=", group.weight.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


#:---- Adjusted
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
             title = "Adjusted",
             clip = c(weight_overall_axis_min, weight_overall_axis_max),
             xticks = seq(weight_overall_axis_min, weight_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Overall population (n=", group.weight.dataset%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_weight_adjusted_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_weight_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Male (n=", group.weight.dataset%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_weight_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_weight_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.weight.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.weight.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.weight.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.weight.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.weight.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.weight.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.weight.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.weight.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(weight_strata_axis_min, weight_strata_axis_max),
             xticks = seq(weight_strata_axis_min, weight_strata_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted weight change (kg)") %>%
  fp_add_header(paste0("Female (n=", group.weight.dataset%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/weight_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("Weight change", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_weight_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_weight_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_weight_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/weight_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("Weight change", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_weight_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_weight_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_weight_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_weight_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_weight_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_weight_adjusted_male
upViewport()

dev.off()





#:--------------------------------------------------------
# eGFR change

## Read in data for egfr
egfr.dataset <- set_up_data_sglt2_glp1(dataset.type = "egfr.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(egfr.change = postegfr - preegfr)


group.egfr.dataset <- group_values(data = egfr.dataset,
                                   variable = "effects",
                                   breaks = interval_breaks) %>%
  drop_na(intervals)


breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.egfr.dataset[,breakdown_adjust], is.factor)

matching_egfr <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.egfr.dataset,
  method = "nearest",
  distance = group.egfr.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_egfr, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_egfr$match.matrix))

group.egfr.dataset.matched <- group.egfr.dataset %>%
  slice(which(matching_egfr$weights == 1))


## Propensity score matching
if (class(try(
  
  # predictions for the egfr adjusted model sex strata
  predictions_egfr_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_psm_1_1 <- vector()
  
  formula <- "egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preegfr, 3) + sex*factor(drugclass)"
  
  models_egfr_psm_1_1 <- stan_glm(formula, 
                                  data = group.egfr.dataset.matched,
                                  family = gaussian(link = "identity"),
                                  prior = normal(0, 2),
                                  prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_psm_1_1, paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                            filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                          type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                sex = "Male",
                                                                                drugclass = "SGLT2",
                                                                                intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                           filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                         type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                lci = quantile(male_glp1, probs = c(0.025)),
                                                                                uci = quantile(male_glp1, probs = c(0.975)),
                                                                                sex = "Male",
                                                                                drugclass = "GLP1",
                                                                                intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                              filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                sex = "Female",
                                                                                drugclass = "SGLT2",
                                                                                intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                             filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                           type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                lci = quantile(female_glp1, probs = c(0.025)),
                                                                                uci = quantile(female_glp1, probs = c(0.975)),
                                                                                sex = "Female",
                                                                                drugclass = "GLP1",
                                                                                intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1 <- predictions_egfr_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1.rds"))
  
  
}

if (class(try(
  
  # predictions for the egfr adjusted model
  predictions_egfr_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_psm_1_1_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                       filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                       mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                     type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_overall <- rbind(predictions_egfr_stan_psm_1_1_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                lci = quantile(sglt2, probs = c(0.025)),
                                                                                                uci = quantile(sglt2, probs = c(0.975)),
                                                                                                drugclass = "SGLT2",
                                                                                                intervals = levels(group.egfr.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                      filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                      mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                    type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_overall <- rbind(predictions_egfr_stan_psm_1_1_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                lci = quantile(glp1, probs = c(0.025)),
                                                                                                uci = quantile(glp1, probs = c(0.975)),
                                                                                                drugclass = "GLP1",
                                                                                                intervals = levels(group.egfr.dataset.matched$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1_overall <- predictions_egfr_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the egfr adjusted model sex strata
  predictions_egfr_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_egfr_stan_psm_1_1_full <- vector()
  
  models_egfr_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1.rds"))
  
  sglt2 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                     mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                   type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_psm_1_1_full <- rbind(predictions_egfr_stan_psm_1_1_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                        lci = quantile(sglt2, probs = c(0.025)),
                                                                                        uci = quantile(sglt2, probs = c(0.975)),
                                                                                        drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                    mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                  type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_psm_1_1_full <- rbind(predictions_egfr_stan_psm_1_1_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                        lci = quantile(glp1, probs = c(0.025)),
                                                                                        uci = quantile(glp1, probs = c(0.975)),
                                                                                        drugclass = "GLP1"))
  
  predictions_egfr_stan_psm_1_1_full <- predictions_egfr_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_full.rds"))
  
  
}

## Propensity score matching + adjustment
if (class(try(
  
  # predictions for the egfr adjusted model sex strata
  predictions_egfr_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_psm_1_1_adjusted <- vector()
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)
  
  formula <- paste0("egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_egfr_psm_1_1_adjusted <- stan_glm(formula, 
                                           data = group.egfr.dataset.matched,
                                           family = gaussian(link = "identity"),
                                           prior = normal(0, 2),
                                           prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                            filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                          type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                  lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                                  uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                                  sex = "Male",
                                                                                                  drugclass = "SGLT2",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                           filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                         type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                  lci = quantile(male_glp1, probs = c(0.025)),
                                                                                                  uci = quantile(male_glp1, probs = c(0.975)),
                                                                                                  sex = "Male",
                                                                                                  drugclass = "GLP1",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                              filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                  lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                                  uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                                  sex = "Female",
                                                                                                  drugclass = "SGLT2",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                             filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                           type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                  lci = quantile(female_glp1, probs = c(0.025)),
                                                                                                  uci = quantile(female_glp1, probs = c(0.975)),
                                                                                                  sex = "Female",
                                                                                                  drugclass = "GLP1",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1_adjusted <- predictions_egfr_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted.rds"))
  
  
}

if (class(try(
  
  # predictions for the egfr adjusted model
  predictions_egfr_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_psm_1_1_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                       filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                       mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                     type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted_overall <- rbind(predictions_egfr_stan_psm_1_1_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                  lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                  uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                  drugclass = "SGLT2",
                                                                                                                  intervals = levels(group.egfr.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                      filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                      mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                    type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_psm_1_1_adjusted_overall <- rbind(predictions_egfr_stan_psm_1_1_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                  lci = quantile(glp1, probs = c(0.025)),
                                                                                                                  uci = quantile(glp1, probs = c(0.975)),
                                                                                                                  drugclass = "GLP1",
                                                                                                                  intervals = levels(group.egfr.dataset.matched$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1_adjusted_overall <- predictions_egfr_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_overall.rds"))
  
}
if (class(try(
  
  # predictions for the egfr adjusted model sex strata
  predictions_egfr_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_egfr_stan_psm_1_1_adjusted_full <- vector()
  
  models_egfr_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                     mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset.matched$drugclass))),
                   type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_psm_1_1_adjusted_full <- rbind(predictions_egfr_stan_psm_1_1_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                          lci = quantile(sglt2, probs = c(0.025)),
                                                                                                          uci = quantile(sglt2, probs = c(0.975)),
                                                                                                          drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                    mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset.matched$drugclass))),
                  type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_psm_1_1_adjusted_full <- rbind(predictions_egfr_stan_psm_1_1_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                          lci = quantile(glp1, probs = c(0.025)),
                                                                                                          uci = quantile(glp1, probs = c(0.975)),
                                                                                                          drugclass = "GLP1"))
  
  predictions_egfr_stan_psm_1_1_adjusted_full <- predictions_egfr_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted_full.rds"))
  
  
}



## Adjustment
if (class(try(
  
  # predictions for the egfr adjustment model sex strata
  predictions_egfr_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_adjusted <- vector()
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.egfr.dataset[,breakdown_adjust], is.factor)
  
  formula <- paste0("egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_egfr_adjusted <- stan_glm(formula, 
                                   data = group.egfr.dataset,
                                   family = gaussian(link = "identity"),
                                   prior = normal(0, 2),
                                   prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_adjusted, paste0(output_path, "/additional_outcomes/models_egfr_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                            filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset$drugclass))),
                          type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                  lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                  uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                  sex = "Male",
                                                                                  drugclass = "SGLT2",
                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                           filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset$drugclass))),
                         type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                  lci = quantile(male_glp1, probs = c(0.025)),
                                                                                  uci = quantile(male_glp1, probs = c(0.975)),
                                                                                  sex = "Male",
                                                                                  drugclass = "GLP1",
                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                              filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                  lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                  uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                  sex = "Female",
                                                                                  drugclass = "SGLT2",
                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                             filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset$drugclass))),
                           type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                  lci = quantile(female_glp1, probs = c(0.025)),
                                                                                  uci = quantile(female_glp1, probs = c(0.975)),
                                                                                  sex = "Female",
                                                                                  drugclass = "GLP1",
                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_adjusted <- predictions_egfr_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the egfr adjusted model
  predictions_egfr_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_egfr_stan_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                       filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                       mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset$drugclass))),
                     type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted_overall <- rbind(predictions_egfr_stan_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                  lci = quantile(sglt2, probs = c(0.025)),
                                                                                                  uci = quantile(sglt2, probs = c(0.975)),
                                                                                                  drugclass = "SGLT2",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
    glp1 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                      filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                      mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset$drugclass))),
                    type = "response") %>%
      rowMeans()
    
    predictions_egfr_stan_adjusted_overall <- rbind(predictions_egfr_stan_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                  lci = quantile(glp1, probs = c(0.025)),
                                                                                                  uci = quantile(glp1, probs = c(0.975)),
                                                                                                  drugclass = "GLP1",
                                                                                                  intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_adjusted_overall <- predictions_egfr_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the egfr adjusted model sex strata
  predictions_egfr_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_egfr_stan_adjusted_full <- vector()
  
  models_egfr_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_egfr_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                     mutate(drugclass = factor("SGLT2", levels = levels(group.egfr.dataset$drugclass))),
                   type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_adjusted_full <- rbind(predictions_egfr_stan_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                          lci = quantile(sglt2, probs = c(0.025)),
                                                                                          uci = quantile(sglt2, probs = c(0.975)),
                                                                                          drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                    mutate(drugclass = factor("GLP1", levels = levels(group.egfr.dataset$drugclass))),
                  type = "response") %>%
    rowMeans()
  
  predictions_egfr_stan_adjusted_full <- rbind(predictions_egfr_stan_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                          lci = quantile(glp1, probs = c(0.025)),
                                                                                          uci = quantile(glp1, probs = c(0.975)),
                                                                                          drugclass = "GLP1"))
  
  predictions_egfr_stan_adjusted_full <- predictions_egfr_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted_full.rds"))
  
  
}



#:--------------------
#:---- PLOTS
#:--------------------

## limits

egfr_overall_axis_min <- plyr::round_any(floor(min(c(predictions_egfr_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                     predictions_egfr_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                     predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                     predictions_egfr_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                     predictions_egfr_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                     predictions_egfr_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

egfr_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_egfr_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_egfr_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                       predictions_egfr_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                       predictions_egfr_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                       predictions_egfr_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


egfr_strata_axis_min <- plyr::round_any(floor(min(c(predictions_egfr_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                    predictions_egfr_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                    predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                    predictions_egfr_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                    predictions_egfr_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                                    predictions_egfr_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

egfr_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_egfr_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_egfr_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_egfr_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                                      predictions_egfr_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                      predictions_egfr_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                      predictions_egfr_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)



#:---- PSM 1:1
plot_egfr_psm_1_1_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1_overall %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Overall population (n=", group.egfr.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_egfr_psm_1_1_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Male (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_egfr_psm_1_1_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Female (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


#:---- PSM 1:1 + adjusted
plot_egfr_psm_1_1_adjusted_overall <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1_adjusted_overall %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1_adjusted_overall %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Overall population (n=", group.egfr.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_egfr_psm_1_1_adjusted_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Male (n=", group.egfr.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_egfr_psm_1_1_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset.matched$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Female (n=", group.egfr.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


#:---- Adjusted
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
             title = "Adjusted",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Overall population (n=", group.egfr.dataset%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_egfr_adjusted_male <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.egfr.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Male (n=", group.egfr.dataset%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_egfr_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)),
  predictions_egfr_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.egfr.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[1])%>%nrow(),")"),
                            ifelse(intervals == levels(group.egfr.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[2])%>%nrow(), ")"), 
                                   ifelse(intervals == levels(group.egfr.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[3])%>%nrow(),")"),
                                          ifelse(intervals == levels(group.egfr.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[4])%>%nrow(),")"), 
                                                 ifelse(intervals == levels(group.egfr.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[5])%>%nrow(),")"),
                                                        ifelse(intervals == levels(group.egfr.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.egfr.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.egfr.dataset$intervals)[6])%>%nrow(),")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(egfr_overall_axis_min, egfr_overall_axis_max),
             xticks = seq(egfr_overall_axis_min, egfr_overall_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Predicted eGFR change") %>%
  fp_add_header(paste0("Female (n=", group.egfr.dataset%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/egfr_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("eGFR change", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_egfr_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_egfr_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_egfr_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/egfr_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("eGFR change", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_egfr_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_egfr_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_egfr_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_egfr_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_egfr_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_egfr_adjusted_male
upViewport()

dev.off()






#:--------------------------------------------------------
# Discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(stopdrug_6m_3mFU = factor(stopdrug_6m_3mFU))

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                              variable = "effects",
                                              breaks = interval_breaks) %>%
  drop_na(intervals)

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

matching_discontinuation <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration +  + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
  data = group.discontinuation.dataset,
  method = "nearest",
  distance = group.discontinuation.dataset[,"prop.score"], 
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_discontinuation, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_discontinuation$match.matrix))

group.discontinuation.dataset.matched <- group.discontinuation.dataset %>%
  slice(which(matching_discontinuation$weights == 1))



## Propensity score matching
if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_psm_1_1 <- vector()
  
  formula <- "stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass)"
  
  models_discontinuation_psm_1_1 <- stan_glm(formula, 
                                             data = group.discontinuation.dataset.matched,
                                             family = binomial(link = "logit"),
                                             prior = normal(0, 2),
                                             prior_intercept = normal(0, 2))
  
  saveRDS(models_discontinuation_psm_1_1, paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                                    filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                    filter(sex == "Male") %>%
                                    mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))) %>%
                                    slice(1),
                                  type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                                      uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                                   filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                   filter(sex == "Male") %>%
                                   mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))) %>%
                                   slice(1),
                                 type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(male_glp1, probs = c(0.025)),
                                                                                                      uci = quantile(male_glp1, probs = c(0.975)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                                      filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                      filter(sex == "Female") %>%
                                      mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))) %>%
                                      slice(1),
                                    type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                                      uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                                     filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                     filter(sex == "Female") %>%
                                     mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))) %>%
                                     slice(1),
                                   type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(female_glp1, probs = c(0.025)),
                                                                                                      uci = quantile(female_glp1, probs = c(0.975)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1 <- predictions_discontinuation_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation adjusted model
  predictions_discontinuation_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_discontinuation_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_psm_1_1_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                               filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                               mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                             type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_overall <- rbind(predictions_discontinuation_stan_psm_1_1_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                      lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                      uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                      drugclass = "SGLT2",
                                                                                                                      intervals = levels(group.discontinuation.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                              filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                              mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_overall <- rbind(predictions_discontinuation_stan_psm_1_1_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                      lci = quantile(glp1, probs = c(0.025)),
                                                                                                                      uci = quantile(glp1, probs = c(0.975)),
                                                                                                                      drugclass = "GLP1",
                                                                                                                      intervals = levels(group.discontinuation.dataset.matched$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1_overall <- predictions_discontinuation_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_discontinuation_stan_psm_1_1_full <- vector()
  
  models_discontinuation_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1.rds"))
  
  sglt2 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                             mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                           type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_psm_1_1_full <- rbind(predictions_discontinuation_stan_psm_1_1_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                              lci = quantile(sglt2, probs = c(0.025)),
                                                                                                              uci = quantile(sglt2, probs = c(0.975)),
                                                                                                              drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                            mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                          type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_psm_1_1_full <- rbind(predictions_discontinuation_stan_psm_1_1_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                              lci = quantile(glp1, probs = c(0.025)),
                                                                                                              uci = quantile(glp1, probs = c(0.975)),
                                                                                                              drugclass = "GLP1"))
  
  
  predictions_discontinuation_stan_psm_1_1_full <- predictions_discontinuation_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_full.rds"))
  
  
}

## Propensity score matching + adjustment
if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_psm_1_1_adjusted <- vector()
  
  formula <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(breakdown_adjust[factors], collapse = "+"))
  
  models_discontinuation_psm_1_1_adjusted <- stan_glm(formula, 
                                                      data = group.discontinuation.dataset.matched,
                                                      family = binomial(link = "logit"),
                                                      prior = normal(0, 2),
                                                      prior_intercept = normal(0, 2))
  
  saveRDS(models_discontinuation_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1_adjusted.rds"))
  
  # #not splines version
  # formula  <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + agetx + t2dmduration + prehba1c + preegfr + prealt + prebmi +", paste(breakdown_adjust[factors], collapse = "+"))
  # 
  # models_discontinuation_psm_1_1_adjusted_no_rcs <- stan_glm(formula,
  #                                                     data = group.discontinuation.dataset.matched,
  #                                                     family = binomial(link = "logit"),
  #                                                     prior = normal(0, 2),
  #                                                     prior_intercept = normal(0, 2))
  # 
  # saveRDS(models_discontinuation_psm_1_1_adjusted_no_rcs, paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1_adjusted_no_rcs.rds"))
  
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                                    filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                    filter(sex == "Male") %>%
                                    mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                                  type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                                        lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                                                        uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                                                        sex = "Male",
                                                                                                                        drugclass = "SGLT2",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                                   filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                   filter(sex == "Male") %>%
                                   mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                                 type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                                        lci = quantile(male_glp1, probs = c(0.025)),
                                                                                                                        uci = quantile(male_glp1, probs = c(0.975)),
                                                                                                                        sex = "Male",
                                                                                                                        drugclass = "GLP1",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                                      filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                      filter(sex == "Female") %>%
                                      mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                                    type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                                        lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                                                        uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                                                        sex = "Female",
                                                                                                                        drugclass = "SGLT2",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                                     filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                                     filter(sex == "Female") %>%
                                     mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                                   type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                                        lci = quantile(female_glp1, probs = c(0.025)),
                                                                                                                        uci = quantile(female_glp1, probs = c(0.975)),
                                                                                                                        sex = "Female",
                                                                                                                        drugclass = "GLP1",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1_adjusted <- predictions_discontinuation_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation adjusted model
  predictions_discontinuation_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_discontinuation_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_psm_1_1_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                               filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                               mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                             type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted_overall <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                                        lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                                        uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                                        drugclass = "SGLT2",
                                                                                                                                        intervals = levels(group.discontinuation.dataset.matched$intervals)[i]))
    
    glp1 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                              filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                              mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_psm_1_1_adjusted_overall <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                                        lci = quantile(glp1, probs = c(0.025)),
                                                                                                                                        uci = quantile(glp1, probs = c(0.975)),
                                                                                                                                        drugclass = "GLP1",
                                                                                                                                        intervals = levels(group.discontinuation.dataset.matched$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1_adjusted_overall <- predictions_discontinuation_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_discontinuation_stan_psm_1_1_adjusted_full <- vector()
  
  models_discontinuation_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                             mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                           type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_psm_1_1_adjusted_full <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                                lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                                uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                                drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                            mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset.matched$drugclass))),
                          type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_psm_1_1_adjusted_full <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                                lci = quantile(glp1, probs = c(0.025)),
                                                                                                                                uci = quantile(glp1, probs = c(0.975)),
                                                                                                                                drugclass = "GLP1"))
  
  
  
  
  
  predictions_discontinuation_stan_psm_1_1_adjusted_full <- predictions_discontinuation_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted_full.rds"))
  
  
}

## Adjustment
if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_adjusted <- vector()
  
  formula <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(breakdown_adjust[factors], collapse = "+"))
  
  models_discontinuation_adjusted <- stan_glm(formula, 
                                              data = group.discontinuation.dataset,
                                              family = binomial(link = "logit"),
                                              prior = normal(0, 2),
                                              prior_intercept = normal(0, 2))
  
  saveRDS(models_discontinuation_adjusted, paste0(output_path, "/additional_outcomes/models_discontinuation_adjusted.rds"))
  
  # #not splines version
  # formula  <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + agetx + t2dmduration + prehba1c + preegfr + prealt + prebmi +", paste(breakdown_adjust[factors], collapse = "+"))
  # 
  # models_discontinuation_adjusted_no_rcs <- stan_glm(formula,
  #                                                     data = group.discontinuation.dataset,
  #                                                     family = binomial(link = "logit"),
  #                                                     prior = normal(0, 2),
  #                                                     prior_intercept = normal(0, 2))
  # 
  # saveRDS(models_discontinuation_adjusted_no_rcs, paste0(output_path, "/additional_outcomes/models_discontinuation_adjusted_no_rcs.rds"))
  
  
  for (i in mnumber) {
    
    male_sglt2 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                                    filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                                    filter(sex == "Male") %>%
                                    mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset$drugclass))),
                                  type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                        lci = quantile(male_sglt2, probs = c(0.025)),
                                                                                                        uci = quantile(male_sglt2, probs = c(0.975)),
                                                                                                        sex = "Male",
                                                                                                        drugclass = "SGLT2",
                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                                   filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                                   filter(sex == "Male") %>%
                                   mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset$drugclass))),
                                 type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                        lci = quantile(male_glp1, probs = c(0.025)),
                                                                                                        uci = quantile(male_glp1, probs = c(0.975)),
                                                                                                        sex = "Male",
                                                                                                        drugclass = "GLP1",
                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                                      filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                                      filter(sex == "Female") %>%
                                      mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset$drugclass))),
                                    type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                        lci = quantile(female_sglt2, probs = c(0.025)),
                                                                                                        uci = quantile(female_sglt2, probs = c(0.975)),
                                                                                                        sex = "Female",
                                                                                                        drugclass = "SGLT2",
                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                                     filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                                     filter(sex == "Female") %>%
                                     mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset$drugclass))),
                                   type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                        lci = quantile(female_glp1, probs = c(0.025)),
                                                                                                        uci = quantile(female_glp1, probs = c(0.975)),
                                                                                                        sex = "Female",
                                                                                                        drugclass = "GLP1",
                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_adjusted <- predictions_discontinuation_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation adjusted model
  predictions_discontinuation_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_discontinuation_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_adjusted.rds"))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_discontinuation_stan_adjusted_overall <- vector()
  
  for (i in mnumber) {
    
    sglt2 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                               filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                               mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset$drugclass))),
                             type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted_overall <- rbind(predictions_discontinuation_stan_adjusted_overall, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                        lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                        uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                        drugclass = "SGLT2",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    glp1 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                              filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                              mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset$drugclass))),
                            type = "response") %>%
      rowMeans()
    
    predictions_discontinuation_stan_adjusted_overall <- rbind(predictions_discontinuation_stan_adjusted_overall, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                        lci = quantile(glp1, probs = c(0.025)),
                                                                                                                        uci = quantile(glp1, probs = c(0.975)),
                                                                                                                        drugclass = "GLP1",
                                                                                                                        intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_adjusted_overall <- predictions_discontinuation_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the discontinuation model sex strata
  predictions_discontinuation_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  predictions_discontinuation_stan_adjusted_full <- vector()
  
  models_discontinuation_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/models_discontinuation_adjusted.rds"))
  
  sglt2 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                             mutate(drugclass = factor("SGLT2", levels = levels(group.discontinuation.dataset$drugclass))),
                           type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_adjusted_full <- rbind(predictions_discontinuation_stan_adjusted_full, cbind(mean = mean(sglt2, na.rm = TRUE),
                                                                                                                lci = quantile(sglt2, probs = c(0.025)),
                                                                                                                uci = quantile(sglt2, probs = c(0.975)),
                                                                                                                drugclass = "SGLT2"))
  
  glp1 <- posterior_epred(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                            mutate(drugclass = factor("GLP1", levels = levels(group.discontinuation.dataset$drugclass))),
                          type = "response") %>%
    rowMeans()
  
  predictions_discontinuation_stan_adjusted_full <- rbind(predictions_discontinuation_stan_adjusted_full, cbind(mean = mean(glp1, na.rm = TRUE),
                                                                                                                lci = quantile(glp1, probs = c(0.025)),
                                                                                                                uci = quantile(glp1, probs = c(0.975)),
                                                                                                                drugclass = "GLP1"))
  
  predictions_discontinuation_stan_adjusted_full <- predictions_discontinuation_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted_full.rds"))
  
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits
discontinuation_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                  predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                  predictions_discontinuation_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                  predictions_discontinuation_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                  predictions_discontinuation_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)

discontinuation_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                 predictions_discontinuation_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                 predictions_discontinuation_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                 predictions_discontinuation_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                                                 predictions_discontinuation_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                                                 predictions_discontinuation_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 5, f = ceiling)



#:---- PSM 1:1
plot_discontinuation_psm_1_1_overall <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", group.discontinuation.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_discontinuation_psm_1_1_male <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Male") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Male (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_discontinuation_psm_1_1_female <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect", sex = "Female") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Female (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


#:---- PSM 1:1 + adjusted
plot_discontinuation_psm_1_1_adjusted_overall <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted_overall %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted_overall %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), 
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(0, discontinuation_overall_axis_max),
             xticks = seq(0, discontinuation_overall_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Overall population (n=", group.discontinuation.dataset.matched%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


plot_discontinuation_psm_1_1_adjusted_male <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Male (n=", group.discontinuation.dataset.matched%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_discontinuation_psm_1_1_adjusted_female <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset.matched$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Female (n=", group.discontinuation.dataset.matched%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



#:---- Adjusted
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
             title = "Adjusted",
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


plot_discontinuation_adjusted_male <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Male") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Male")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Male (n=", group.discontinuation.dataset%>%filter(sex=="Male")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_discontinuation_adjusted_female <- rbind(
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 150, lci = 150, uci = 150, sex = "Female", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6)) %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100),
  predictions_discontinuation_stan_adjusted_full %>% cbind(intervals = "Average treatment effect", sex = "Female") %>%
    mutate(mean = as.numeric(mean)*100,
           lci = as.numeric(lci)*100,
           uci = as.numeric(uci)*100)
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == levels(group.discontinuation.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[1])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.discontinuation.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[2])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.discontinuation.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[3])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.discontinuation.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[4])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.discontinuation.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[5])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.discontinuation.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%nrow(), ", event = ", group.discontinuation.dataset%>%filter(sex=="Female")%>%filter(intervals==levels(group.discontinuation.dataset$intervals)[6])%>%filter(stopdrug_6m_3mFU==1)%>%nrow(), ")"), intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(0, discontinuation_strata_axis_max),
             xticks = seq(0, discontinuation_strata_axis_max, 5),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation (%)") %>%
  fp_add_header(paste0("Female (n=", group.discontinuation.dataset%>%filter(sex=="Female")%>%nrow(), ")")) %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))



pdf(width = 7, height = 12, "Plots/discontinuation_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("Discontinuation", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_discontinuation_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_discontinuation_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_discontinuation_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/discontinuation_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("Discontinuation", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_discontinuation_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_discontinuation_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_discontinuation_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_discontinuation_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_discontinuation_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_discontinuation_adjusted_male
upViewport()

dev.off()




#:--------------------------------------------------------
#:-------- COMORBIDITIES
# 
# # I have to refit the propensity score model whilst including qrisk2 as a variable.
# # read variables used in initial propensity score
# vs_bart_ps_model <- readRDS(paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))
# 
# variables_chosen <- unique(gsub("_.*", "", vs_bart_ps_model$important_vars_local_names))
# 
# # fit new propensity score model
# ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train") %>%
#   left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
#               select(patid, pated, qrisk2_10yr_score), by = c("patid", "pated"))
# ## Fit initial model using all the available variables to estimate HbA1c outcome
# if (class(try(
#   
#   bart_qrisk_ps_model <- readRDS(paste0(output_path, "/additional_outcomes/bart_qrisk_ps_model.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bart_qrisk_ps_model <- bartMachine::bartMachine(X = ps.model.train %>%
#                                                     select(
#                                                       all_of(c(variables_chosen, "qrisk2_10yr_score"))
#                                                     ),
#                                                   y = ps.model.train[,"drugclass"] %>%
#                                                     unlist(),
#                                                   num_trees = 50,
#                                                   use_missing_data = TRUE,
#                                                   num_burn_in = 5000,
#                                                   num_iterations_after_burn_in = 5000,
#                                                   serialize = TRUE)
#   
#   saveRDS(bart_qrisk_ps_model, paste0(output_path, "/additional_outcomes/bart_qrisk_ps_model.rds"))
#   
# }
# 
# 
# # Prop scores for test dataset
# 
# ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test") %>%
#   left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
#               select(patid, pated, qrisk2_10yr_score), by = c("patid", "pated"))
# 
# # calculate prop score
# if (class(try(
#   
#   prop_score_qrisk_testing_data <- readRDS(paste0(output_path, "/additional_outcomes/prop_score_qrisk_testing_data.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   prop_score_qrisk_testing_data <- predict(bart_qrisk_ps_model, ps.model.test %>%
#                                        select(
#                                          colnames(bart_qrisk_ps_model$X)
#                                        ))
#   
#   saveRDS(prop_score_qrisk_testing_data, paste0(output_path, "/additional_outcomes/prop_score_qrisk_testing_data.rds"))
#   
# }
# 
# 
# if (class(try(
#   
#   patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   # Prop scores for train dataset
#   patient_prop_scores_qrisk <- ps.model.train %>%
#     select(patid, pated) %>%
#     cbind(prop.score = bart_qrisk_ps_model$p_hat_train)
#   
#   
#   patient_prop_scores_qrisk <- patient_prop_scores_qrisk %>%
#     rbind(
#       ps.model.test %>%
#         select(patid, pated) %>%
#         cbind(prop.score = prop_score_qrisk_testing_data)
#     ) %>%
#     as.data.frame()
#   
#   saveRDS(patient_prop_scores_qrisk, paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))
#   
# }

#:--------------------------
# No comorbidities

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

no_co.dataset <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

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

# require(cobalt)
# cobalt::love.plot(matching_no_co, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_no_co$match.matrix))

group.no_co.dataset.matched <- group.no_co.dataset %>%
  slice(which(matching_no_co$weights == 1))


# CVD survival analysis


## Propensity score matching

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.no_co.dataset.matched)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#            data = group.no_co.dataset.matched,
#            # fun = "event",
#            palette = c("dodgerblue2", "#f1a340"),
#            surv.median.line = "none",
#            ylim = c(0.6, 1),
#            censor = FALSE,
#            conf.int = TRUE,
#            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#            legend.title = "Therapy",
#            legend = "bottom",
#            title = paste0("CVD outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#            subtitle = "Propensity score matching")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.no_co.dataset.matched$intervals)
# 
# 
# kaplan_meier_no_co_cvd_psm_1_1 <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                          labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                             data = group.no_co.dataset.matched,
#                             fun = "event",
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0, 0.4),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("CVD outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Propensity score matching")
# 
# cumulative_incidence_no_co_cvd_psm_1_1 <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                           labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_psm_1_1_effect <- coxph(formula(formula_freq_1),
                                                        data = group.no_co.dataset.matched)

formula_freq_2 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_psm_1_1_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                                 data = group.no_co.dataset.matched)

formula_freq_3 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_psm_1_1_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                                 data = group.no_co.dataset.matched)

# anova(models_no_co_cvd_spline_effects_psm_1_1_effect, models_no_co_cvd_spline_effects_psm_1_1_spline_effect_3, models_no_co_cvd_spline_effects_psm_1_1_spline_effect_5)


dataset.used <- group.no_co.dataset.matched %>%
  select(postdrug_mace_censtime_yrs, postdrug_mace_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')


m1 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_cvd_psm_1_1_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                  rel_heights = c(1,0.4), rel_widths = c(1,1))


m2 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_cvd_psm_1_1_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                  rel_heights = c(1,0.4), rel_widths = c(1,1))




#--- Cox  
if (class(try(
  
  # predictions for the CVD outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_cvd_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1 <- vector()
  
  
  ###:-------------------------
  ## test - comparison with male, interval 1
  
  # interim.dataset <- group.no_co.dataset.matched %>%
  #   select(drugclass, intervals, sex, qrisk2_10yr_score) %>%
  #   cbind(time = group.no_co.dataset.matched$postdrug_mace_censtime_yrs,
  #         censored = group.no_co.dataset.matched$postdrug_mace_censvar,
  #         qrisk2_10yr_score_1 = rcs(group.no_co.dataset.matched$qrisk2_10yr_score, 3)[,1],
  #         qrisk2_10yr_score_2 = rcs(group.no_co.dataset.matched$qrisk2_10yr_score, 3)[,2]) %>%
  #   as.data.frame() %>%
  #   drop_na() %>%
  #   mutate(sex = relevel(sex, ref = "Male"),
  #          intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[1]))
  # 
  # model.interim.dataset <- model.matrix(~time + censored + drugclass, data = interim.dataset) %>%
  #   as.data.frame()
  # 
  # N <- nrow(model.interim.dataset)
  # X <- as.matrix(pull(model.interim.dataset, drugclassSGLT2))
  # is_censored <- pull(model.interim.dataset,censored)==0
  # times <- pull(model.interim.dataset,time)
  # msk_censored <- is_censored == 1
  # N_censored <- sum(msk_censored)
  # 
  # stan_data <- list(N_uncensored=N-N_censored, 
  #                   N_censored=N_censored, 
  #                   X_censored=as.matrix(X[msk_censored,]),
  #                   X_uncensored=as.matrix(X[!msk_censored,]),
  #                   times_censored=times[msk_censored],
  #                   times_uncensored = times[!msk_censored],
  #                   NC=ncol(X)
  # )
  # 
  # fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
  # 
  # 
  # formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass)"
  # 
  # fit_lin <- coxph(formula(formula_freq),
  #       data = group.no_co.dataset.matched %>%
  #         mutate(sex = relevel(sex, ref = "Male"),
  #                intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[1])))
  

  # formula <- "time | cens(censored) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + qrisk2_10yr_score_1 + qrisk2_10yr_score_2"
  # 
  # models_test <- brms::brm(formula = formula(formula),
  #                          data = interim.dataset,
  #                          family = "exponential",
  #                          chains = 1)

  ###:-------------------------
  
  # formula <- 'postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)'
  # 
  # models_no_co_cvd_psm_1_1 <- brms::brm(formula = formula(formula),
  #                                       data = group.no_co.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)"
  
  models_no_co_cvd_psm_1_1_male <- vector("list", quantiles)
  
  models_no_co_cvd_psm_1_1_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_psm_1_1_male[[i]] <- coxph(formula(formula_freq),
                                                data = group.no_co.dataset.matched %>%
                                                  mutate(sex = relevel(sex, ref = "Male"),
                                                         intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1 <- rbind(predictions_no_co_cvd_stan_psm_1_1, cbind(mean = models_no_co_cvd_psm_1_1_male[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_cvd_psm_1_1_male[[i]])[1,1],
                                                                                          uci = confint(models_no_co_cvd_psm_1_1_male[[i]])[1,2],
                                                                                          sex = "Male",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    models_no_co_cvd_psm_1_1_female[[i]] <- coxph(formula(formula_freq),
                                                  data = group.no_co.dataset.matched %>%
                                                    mutate(sex = relevel(sex, ref = "Female"),
                                                           intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1 <- rbind(predictions_no_co_cvd_stan_psm_1_1, cbind(mean = models_no_co_cvd_psm_1_1_female[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_cvd_psm_1_1_female[[i]])[1,1],
                                                                                          uci = confint(models_no_co_cvd_psm_1_1_female[[i]])[1,2],
                                                                                          sex = "Female",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_psm_1_1 <- predictions_no_co_cvd_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1.rds"))
  
}

# No sex strata
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1_overall <- vector()
  
  models_no_co_cvd_psm_1_1_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_psm_1_1_overall[[i]] <- coxph(formula(formula_freq),
                                                   data = group.no_co.dataset.matched %>%
                                                     mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1_overall <- rbind(predictions_no_co_cvd_stan_psm_1_1_overall, cbind(mean = models_no_co_cvd_psm_1_1_overall[[i]]$coefficients[1],
                                                                                                          lci = confint(models_no_co_cvd_psm_1_1_overall[[i]])[1,1],
                                                                                                          uci = confint(models_no_co_cvd_psm_1_1_overall[[i]])[1,2],
                                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_psm_1_1_overall <- predictions_no_co_cvd_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_overall.rds"))
  
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1_full <- vector()
  
  models_no_co_cvd_psm_1_1_full <- coxph(formula(formula_freq),
                                         data = group.no_co.dataset.matched)
  
  predictions_no_co_cvd_stan_psm_1_1_full <- rbind(predictions_no_co_cvd_stan_psm_1_1_full, cbind(mean = models_no_co_cvd_psm_1_1_full$coefficients[1],
                                                                                                  lci = confint(models_no_co_cvd_psm_1_1_full)[1,1],
                                                                                                  uci = confint(models_no_co_cvd_psm_1_1_full)[1,2]))
  
  predictions_no_co_cvd_stan_psm_1_1_full <- predictions_no_co_cvd_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_full.rds"))
  
}

## Propensity score matching + adjusted
if (class(try(
  
  # predictions for the CVD outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_cvd_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_no_co_cvd_psm_1_1_adjusted <- brms::brm(formula = formula(formula),
  #                                       data = group.no_co.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_no_co_cvd_psm_1_1_adjusted_male <- vector("list", quantiles)
  
  models_no_co_cvd_psm_1_1_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_psm_1_1_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                         data = group.no_co.dataset.matched %>%
                                                           mutate(sex = relevel(sex, ref = "Male"),
                                                                  intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1_adjusted <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_cvd_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_cvd_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                            sex = "Male",
                                                                                                            intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    models_no_co_cvd_psm_1_1_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                           data = group.no_co.dataset.matched %>%
                                                             mutate(sex = relevel(sex, ref = "Female"),
                                                                    intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1_adjusted <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_cvd_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_cvd_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                            sex = "Female",
                                                                                                            intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted <- predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- vector()
  
  models_no_co_cvd_psm_1_1_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_psm_1_1_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                            data = group.no_co.dataset.matched %>%
                                                              mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted_overall, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_overall[[i]]$coefficients[1],
                                                                                                                            lci = confint(models_no_co_cvd_psm_1_1_adjusted_overall[[i]])[1,1],
                                                                                                                            uci = confint(models_no_co_cvd_psm_1_1_adjusted_overall[[i]])[1,2],
                                                                                                                            intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall <- predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_overall.rds"))
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- vector()
  
  models_no_co_cvd_psm_1_1_adjusted_full <- coxph(formula(formula_freq),
                                                  data = group.no_co.dataset.matched)
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted_full, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_full$coefficients[1],
                                                                                                                    lci = confint(models_no_co_cvd_psm_1_1_adjusted_full)[1,1],
                                                                                                                    uci = confint(models_no_co_cvd_psm_1_1_adjusted_full)[1,2]))
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full <- predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted_full.rds"))
  
}


## Adjusted

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.no_co.dataset)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.no_co.dataset,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.60, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("CVD outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Overall population")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.no_co.dataset$intervals)
# 
# 
# kaplan_meier_no_co_cvd_adjusted <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                   labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.no_co.dataset,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.4),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("CVD outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Overall population")
# 
# cumulative_incidence_no_co_cvd_adjusted <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                    labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_adjusted_effect <- coxph(formula(formula_freq_1),
                                                         data = group.no_co.dataset)

formula_freq_2 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_adjusted_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                                  data = group.no_co.dataset)

formula_freq_3 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_no_co_cvd_spline_effects_adjusted_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                                  data = group.no_co.dataset)

# anova(models_no_co_cvd_spline_effects_adjusted_effect, models_no_co_cvd_spline_effects_adjusted_spline_effect_3, models_no_co_cvd_spline_effects_adjusted_spline_effect_5)


dataset.used <- group.no_co.dataset %>%
  select(postdrug_mace_censtime_yrs, postdrug_mace_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Overall population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_cvd_adjusted_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                   rel_heights = c(1,0.4), rel_widths = c(1,1))


m2 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Overall population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_cvd_adjusted_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                   rel_heights = c(1,0.4), rel_widths = c(1,1))



#--- Cox  
if (class(try(
  
  # predictions for the CVD outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_cvd_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_no_co_cvd_adjusted <- brms::brm(formula = formula(formula),
  #                                        data = group.no_co.dataset,
  #                                        family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_no_co_cvd_adjusted_male <- vector("list", quantiles)
  
  models_no_co_cvd_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                 data = group.no_co.dataset %>%
                                                   mutate(sex = relevel(sex, ref = "Male"),
                                                          intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_cvd_stan_adjusted <- rbind(predictions_no_co_cvd_stan_adjusted, cbind(mean = models_no_co_cvd_adjusted_male[[i]]$coefficients[1],
                                                                                            lci = confint(models_no_co_cvd_adjusted_male[[i]])[1,1],
                                                                                            uci = confint(models_no_co_cvd_adjusted_male[[i]])[1,2],
                                                                                            sex = "Male",
                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))
    
    models_no_co_cvd_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                   data = group.no_co.dataset %>%
                                                     mutate(sex = relevel(sex, ref = "Female"),
                                                            intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_cvd_stan_adjusted <- rbind(predictions_no_co_cvd_stan_adjusted, cbind(mean = models_no_co_cvd_adjusted_female[[i]]$coefficients[1],
                                                                                            lci = confint(models_no_co_cvd_adjusted_female[[i]])[1,1],
                                                                                            uci = confint(models_no_co_cvd_adjusted_female[[i]])[1,2],
                                                                                            sex = "Female",
                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_adjusted <- predictions_no_co_cvd_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_adjusted_overall <- vector()
  
  models_no_co_cvd_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_cvd_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                    data = group.no_co.dataset %>%
                                                      mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_cvd_stan_adjusted_overall <- rbind(predictions_no_co_cvd_stan_adjusted_overall, cbind(mean = models_no_co_cvd_adjusted_overall[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_cvd_adjusted_overall[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_cvd_adjusted_overall[[i]])[1,2],
                                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))
    
  }
  
  predictions_no_co_cvd_stan_adjusted_overall <- predictions_no_co_cvd_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_overall.rds"))
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_cvd_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_cvd_stan_adjusted_full <- vector()
  
  models_no_co_cvd_adjusted_full <- coxph(formula(formula_freq),
                                          data = group.no_co.dataset)
  
  predictions_no_co_cvd_stan_adjusted_full <- rbind(predictions_no_co_cvd_stan_adjusted_full, cbind(mean = models_no_co_cvd_adjusted_full$coefficients[1],
                                                                                                    lci = confint(models_no_co_cvd_adjusted_full)[1,1],
                                                                                                    uci = confint(models_no_co_cvd_adjusted_full)[1,2]))
  
  predictions_no_co_cvd_stan_adjusted_full <- predictions_no_co_cvd_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_cvd_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted_full.rds"))
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits
# no_co_cvd_overall_axis_min <- plyr::round_any(floor(min(c(predictions_no_co_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_cvd_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                        predictions_no_co_cvd_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# no_co_cvd_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_no_co_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                             predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                             predictions_no_co_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                             predictions_no_co_cvd_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                          predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                          predictions_no_co_cvd_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)
# 
# 
# no_co_cvd_strata_axis_min <- plyr::round_any(floor(min(c(predictions_no_co_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                       predictions_no_co_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                       predictions_no_co_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# no_co_cvd_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_no_co_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                         predictions_no_co_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                         predictions_no_co_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)


# Propensity score matching
plot_no_co_cvd_psm_1_1_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_psm_1_1_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_psm_1_1_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Propensity score matching + adjusted
plot_no_co_cvd_psm_1_1_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_psm_1_1_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_psm_1_1_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Adjusted
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
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_cvd_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_cvd_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/no_co_cvd_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("CVD outcomes for no CVD/HF/CKD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_no_co_cvd_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_no_co_cvd_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_no_co_cvd_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/no_co_cvd_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("CVD outcomes for no CVD/HF/CKD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_no_co_cvd_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_no_co_cvd_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_no_co_cvd_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_no_co_cvd_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_no_co_cvd_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_no_co_cvd_adjusted_male
upViewport()

dev.off()


# Heart failure survival analysis

## Propensity score matching

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.no_co.dataset.matched)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.no_co.dataset.matched,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.50, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("HF outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Propensity score matching")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.no_co.dataset.matched$intervals)
# 
# 
# kaplan_meier_no_co_hf_psm_1_1 <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                  labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.no_co.dataset.matched,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.5),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("HF outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Propensity score matching")
# 
# cumulative_incidence_no_co_hf_psm_1_1 <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                          labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_psm_1_1_effect <- coxph(formula(formula_freq_1),
                                                       data = group.no_co.dataset.matched)

formula_freq_2 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_psm_1_1_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                                data = group.no_co.dataset.matched)

formula_freq_3 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_psm_1_1_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                                data = group.no_co.dataset.matched)

# anova(models_no_co_hf_spline_effects_psm_1_1_effect, models_no_co_hf_spline_effects_psm_1_1_spline_effect_3, models_no_co_hf_spline_effects_psm_1_1_spline_effect_5)


dataset.used <- group.no_co.dataset.matched %>%
  select(postdrug_hf_censtime_yrs, postdrug_hf_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_hf_psm_1_1_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                 rel_heights = c(1,0.4), rel_widths = c(1,1))


m2 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.no_co.dataset.matched%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_hf_psm_1_1_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                 rel_heights = c(1,0.4), rel_widths = c(1,1))


#--- Cox
if (class(try(
  
  # predictions for the HF outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_hf_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1 <- vector()
  
  # formula <- 'postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)'
  # 
  # models_no_co_hf_psm_1_1 <- brms::brm(formula = formula(formula),
  #                                       data = group.no_co.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)"
  
  models_no_co_hf_psm_1_1_male <- vector("list", quantiles)
  
  models_no_co_hf_psm_1_1_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_psm_1_1_male[[i]] <- coxph(formula(formula_freq),
                                               data = group.no_co.dataset.matched %>%
                                                 mutate(sex = relevel(sex, ref = "Male"),
                                                        intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1 <- rbind(predictions_no_co_hf_stan_psm_1_1, cbind(mean = models_no_co_hf_psm_1_1_male[[i]]$coefficients[1],
                                                                                        lci = confint(models_no_co_hf_psm_1_1_male[[i]])[1,1],
                                                                                        uci = confint(models_no_co_hf_psm_1_1_male[[i]])[1,2],
                                                                                        sex = "Male",
                                                                                        intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    models_no_co_hf_psm_1_1_female[[i]] <- coxph(formula(formula_freq),
                                                 data = group.no_co.dataset.matched %>%
                                                   mutate(sex = relevel(sex, ref = "Female"),
                                                          intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1 <- rbind(predictions_no_co_hf_stan_psm_1_1, cbind(mean = models_no_co_hf_psm_1_1_female[[i]]$coefficients[1],
                                                                                        lci = confint(models_no_co_hf_psm_1_1_female[[i]])[1,1],
                                                                                        uci = confint(models_no_co_hf_psm_1_1_female[[i]])[1,2],
                                                                                        sex = "Female",
                                                                                        intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_psm_1_1 <- predictions_no_co_hf_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1_overall <- vector()
  
  models_no_co_hf_psm_1_1_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_psm_1_1_overall[[i]] <- coxph(formula(formula_freq),
                                                  data = group.no_co.dataset.matched %>%
                                                    mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1_overall <- rbind(predictions_no_co_hf_stan_psm_1_1_overall, cbind(mean = models_no_co_hf_psm_1_1_overall[[i]]$coefficients[1],
                                                                                                        lci = confint(models_no_co_hf_psm_1_1_overall[[i]])[1,1],
                                                                                                        uci = confint(models_no_co_hf_psm_1_1_overall[[i]])[1,2],
                                                                                                        intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_psm_1_1_overall <- predictions_no_co_hf_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_overall.rds"))
  
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1_full <- vector()
  
  models_no_co_hf_psm_1_1_full <- coxph(formula(formula_freq),
                                        data = group.no_co.dataset.matched)
  
  predictions_no_co_hf_stan_psm_1_1_full <- rbind(predictions_no_co_hf_stan_psm_1_1_full, cbind(mean = models_no_co_hf_psm_1_1_full$coefficients[1],
                                                                                                lci = confint(models_no_co_hf_psm_1_1_full)[1,1],
                                                                                                uci = confint(models_no_co_hf_psm_1_1_full)[1,2]))
  
  predictions_no_co_hf_stan_psm_1_1_full <- predictions_no_co_hf_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_full.rds"))
  
}

## Propensity score matching + adjusted
if (class(try(
  
  # predictions for the HF outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_hf_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_no_co_hf_psm_1_1_adjusted <- brms::brm(formula = formula(formula),
  #                                       data = group.no_co.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_no_co_hf_psm_1_1_adjusted_male <- vector("list", quantiles)
  
  models_no_co_hf_psm_1_1_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_psm_1_1_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                        data = group.no_co.dataset.matched %>%
                                                          mutate(sex = relevel(sex, ref = "Male"),
                                                                 intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1_adjusted <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted, cbind(mean = models_no_co_hf_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                          lci = confint(models_no_co_hf_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                          uci = confint(models_no_co_hf_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                          sex = "Male",
                                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    models_no_co_hf_psm_1_1_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                          data = group.no_co.dataset.matched %>%
                                                            mutate(sex = relevel(sex, ref = "Female"),
                                                                   intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1_adjusted <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted, cbind(mean = models_no_co_hf_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                          lci = confint(models_no_co_hf_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                          uci = confint(models_no_co_hf_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                          sex = "Female",
                                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_psm_1_1_adjusted <- predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- vector()
  
  models_no_co_hf_psm_1_1_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_psm_1_1_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                           data = group.no_co.dataset.matched %>%
                                                             mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset.matched$intervals)[i])))
    
    predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted_overall, cbind(mean = models_no_co_hf_psm_1_1_adjusted_overall[[i]]$coefficients[1],
                                                                                                                          lci = confint(models_no_co_hf_psm_1_1_adjusted_overall[[i]])[1,1],
                                                                                                                          uci = confint(models_no_co_hf_psm_1_1_adjusted_overall[[i]])[1,2],
                                                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall <- predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_overall.rds"))
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_psm_1_1_adjusted_full <- vector()
  
  models_no_co_hf_psm_1_1_adjusted_full <- coxph(formula(formula_freq),
                                                 data = group.no_co.dataset.matched)
  
  predictions_no_co_hf_stan_psm_1_1_adjusted_full <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted_full, cbind(mean = models_no_co_hf_psm_1_1_adjusted_full$coefficients[1],
                                                                                                                  lci = confint(models_no_co_hf_psm_1_1_adjusted_full)[1,1],
                                                                                                                  uci = confint(models_no_co_hf_psm_1_1_adjusted_full)[1,2]))
  
  predictions_no_co_hf_stan_psm_1_1_adjusted_full <- predictions_no_co_hf_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted_full.rds"))
  
}



## Adjusted

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.no_co.dataset)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.no_co.dataset,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.50, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("HF outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Overall population")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.no_co.dataset$intervals)
# 
# 
# kaplan_meier_no_co_hf_adjusted <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                  labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.no_co.dataset,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.5),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("HF outcomes for no CVD/HF/CKD population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Overall population")
# 
# cumulative_incidence_no_co_hf_adjusted <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                           labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_adjusted_effect <- coxph(formula(formula_freq_1),
                                                        data = group.no_co.dataset)

formula_freq_2 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_adjusted_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                                 data = group.no_co.dataset)

formula_freq_3 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_no_co_hf_spline_effects_adjusted_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                                 data = group.no_co.dataset)

# anova(models_no_co_hf_spline_effects_adjusted_effect, models_no_co_hf_spline_effects_adjusted_spline_effect_3, models_no_co_hf_spline_effects_adjusted_spline_effect_5)


dataset.used <- group.no_co.dataset %>%
  select(postdrug_hf_censtime_yrs, postdrug_hf_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Overall population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_hf_adjusted_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                  rel_heights = c(1,0.4), rel_widths = c(1,1))

m2 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.no_co.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.no_co.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no CVD/HF/CKD population") +
  labs(subtitle = paste0("Overall population (n=", group.no_co.dataset%>%nrow(), ", event = ", group.no_co.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.no_co.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_no_co_hf_adjusted_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                  rel_heights = c(1,0.4), rel_widths = c(1,1))

#--- Cox
if (class(try(
  
  # predictions for the HF outcomes in population with no CVD/HF/CKD strata sex and intervals
  predictions_no_co_hf_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_no_co_hf_adjusted <- brms::brm(formula = formula(formula),
  #                                        data = group.no_co.dataset,
  #                                        family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_no_co_hf_adjusted_male <- vector("list", quantiles)
  
  models_no_co_hf_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                data = group.no_co.dataset %>%
                                                  mutate(sex = relevel(sex, ref = "Male"),
                                                         intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_hf_stan_adjusted <- rbind(predictions_no_co_hf_stan_adjusted, cbind(mean = models_no_co_hf_adjusted_male[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_hf_adjusted_male[[i]])[1,1],
                                                                                          uci = confint(models_no_co_hf_adjusted_male[[i]])[1,2],
                                                                                          sex = "Male",
                                                                                          intervals = levels(group.no_co.dataset$intervals)[i]))
    
    models_no_co_hf_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                  data = group.no_co.dataset %>%
                                                    mutate(sex = relevel(sex, ref = "Female"),
                                                           intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_hf_stan_adjusted <- rbind(predictions_no_co_hf_stan_adjusted, cbind(mean = models_no_co_hf_adjusted_female[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_hf_adjusted_female[[i]])[1,1],
                                                                                          uci = confint(models_no_co_hf_adjusted_female[[i]])[1,2],
                                                                                          sex = "Female",
                                                                                          intervals = levels(group.no_co.dataset$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_adjusted <- predictions_no_co_hf_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_adjusted_overall <- vector()
  
  models_no_co_hf_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_no_co_hf_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                   data = group.no_co.dataset %>%
                                                     mutate(intervals = relevel(intervals, ref = levels(group.no_co.dataset$intervals)[i])))
    
    predictions_no_co_hf_stan_adjusted_overall <- rbind(predictions_no_co_hf_stan_adjusted_overall, cbind(mean = models_no_co_hf_adjusted_overall[[i]]$coefficients[1],
                                                                                                          lci = confint(models_no_co_hf_adjusted_overall[[i]])[1,1],
                                                                                                          uci = confint(models_no_co_hf_adjusted_overall[[i]])[1,2],
                                                                                                          intervals = levels(group.no_co.dataset$intervals)[i]))
    
  }
  
  predictions_no_co_hf_stan_adjusted_overall <- predictions_no_co_hf_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_overall.rds"))
  
}
# Full population (interval as a continuous, without subgrouping individuals)
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD/HF/CKD
  predictions_no_co_hf_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_no_co_hf_stan_adjusted_full <- vector()
  
  models_no_co_hf_adjusted_full <- coxph(formula(formula_freq),
                                         data = group.no_co.dataset)
  
  predictions_no_co_hf_stan_adjusted_full <- rbind(predictions_no_co_hf_stan_adjusted_full, cbind(mean = models_no_co_hf_adjusted_full$coefficients[1],
                                                                                                  lci = confint(models_no_co_hf_adjusted_full)[1,1],
                                                                                                  uci = confint(models_no_co_hf_adjusted_full)[1,2]))
  
  predictions_no_co_hf_stan_adjusted_full <- predictions_no_co_hf_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_no_co_hf_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted_full.rds"))
  
}



#:--------------------
#:---- PLOTS
#:--------------------

## limits

# no_co_hf_overall_axis_min <- plyr::round_any(floor(min(c(predictions_no_co_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_hf_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                           predictions_no_co_hf_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# no_co_hf_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_no_co_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                            predictions_no_co_hf_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                             predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                             predictions_no_co_hf_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)
# 
# 
# no_co_hf_strata_axis_min <- plyr::round_any(floor(min(c(predictions_no_co_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_no_co_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_no_co_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_no_co_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                          predictions_no_co_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# no_co_hf_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_no_co_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_no_co_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                            predictions_no_co_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                            predictions_no_co_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)



# Propensity score matching
plot_no_co_hf_psm_1_1_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_overall %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_psm_1_1_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_psm_1_1_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Propensity score matching + adjusted
plot_no_co_hf_psm_1_1_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_psm_1_1_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_psm_1_1_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset.matched%>%filter(intervals==levels(group.no_co.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Adjusted
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
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.no_co.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.no_co.dataset%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_no_co_hf_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_no_co_hf_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_no_co_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.no_co.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.no_co.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.no_co.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.no_co.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.no_co.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.no_co.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.no_co.dataset%>%filter(intervals==levels(group.no_co.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.no_co.dataset%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/no_co_hf_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HF outcomes for no CVD/HF/CKD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_no_co_hf_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_no_co_hf_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_no_co_hf_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/no_co_hf_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HF outcomes for no CVD/HF/CKD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_no_co_hf_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_no_co_hf_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_no_co_hf_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_no_co_hf_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_no_co_hf_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_no_co_hf_adjusted_male
upViewport()

dev.off()



#:--------------------------
# Cardiovascular outcome

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

cvd.dataset <- set_up_data_sglt2_glp1(dataset.type="cvd.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.cvd.dataset <- group_values(data = cvd.dataset,
                                  variable = "effects",
                                  breaks = interval_breaks) %>%
  drop_na(intervals)

matching_cvd <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.cvd.dataset,
  method = "nearest",
  distance = group.cvd.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_cvd, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_cvd$match.matrix))

group.cvd.dataset.matched <- group.cvd.dataset %>%
  slice(which(matching_cvd$weights == 1))

# CVD survival analysis

## Propensity score matching

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.cvd.dataset.matched)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.cvd.dataset.matched,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.80, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("CVD outcomes for CVD population (n=", group.cvd.dataset.matched%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Propensity score matching")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.cvd.dataset.matched$intervals)
# 
# 
# kaplan_meier_cvd_cvd_psm_1_1 <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                   labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.cvd.dataset.matched,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.2),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("CVD outcomes for CVD population (n=", group.cvd.dataset.matched%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Propensity score matching")
# 
# cumulative_incidence_cvd_cvd_psm_1_1 <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                         labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_psm_1_1_effect <- coxph(formula(formula_freq_1),
                                                      data = group.cvd.dataset.matched)

formula_freq_2 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_psm_1_1_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                               data = group.cvd.dataset.matched)

formula_freq_3 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_psm_1_1_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                               data = group.cvd.dataset.matched)

# anova(models_cvd_cvd_spline_effects_psm_1_1_effect, models_cvd_cvd_spline_effects_psm_1_1_spline_effect_3, models_cvd_cvd_spline_effects_psm_1_1_spline_effect_5)


dataset.used <- group.cvd.dataset.matched %>%
  select(postdrug_mace_censtime_yrs, postdrug_mace_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.cvd.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.cvd.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.cvd.dataset.matched%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.cvd.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_cvd_cvd_psm_1_1_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                rel_heights = c(1,0.4), rel_widths = c(1,1))

m2 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.cvd.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.cvd.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.cvd.dataset.matched%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.cvd.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_cvd_cvd_psm_1_1_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                rel_heights = c(1,0.4), rel_widths = c(1,1))

#--- Cox  
if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD strata sex and intervals
  predictions_cvd_cvd_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1 <- vector()
  
  # formula <- 'postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)'
  # 
  # models_cvd_cvd_psm_1_1 <- brms::brm(formula = formula(formula),
  #                                       data = group.cvd.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)"
  
  models_cvd_cvd_psm_1_1_male <- vector("list", quantiles)
  
  models_cvd_cvd_psm_1_1_female <- vector("list", quantiles)
  
  
  for (i in mnumber) {
    
    models_cvd_cvd_psm_1_1_male[[i]] <- coxph(formula(formula_freq),
                                              data = group.cvd.dataset.matched %>%
                                                mutate(sex = relevel(sex, ref = "Male"),
                                                       intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1 <- rbind(predictions_cvd_cvd_stan_psm_1_1, cbind(mean = models_cvd_cvd_psm_1_1_male[[i]]$coefficients[1],
                                                                                      lci = confint(models_cvd_cvd_psm_1_1_male[[i]])[1,1],
                                                                                      uci = confint(models_cvd_cvd_psm_1_1_male[[i]])[1,2],
                                                                                      sex = "Male",
                                                                                      intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    models_cvd_cvd_psm_1_1_female[[i]] <- coxph(formula(formula_freq),
                                                data = group.cvd.dataset.matched %>%
                                                  mutate(sex = relevel(sex, ref = "Female"),
                                                         intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1 <- rbind(predictions_cvd_cvd_stan_psm_1_1, cbind(mean = models_cvd_cvd_psm_1_1_female[[i]]$coefficients[1],
                                                                                      lci = confint(models_cvd_cvd_psm_1_1_female[[i]])[1,1],
                                                                                      uci = confint(models_cvd_cvd_psm_1_1_female[[i]])[1,2],
                                                                                      sex = "Female",
                                                                                      intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_psm_1_1 <- predictions_cvd_cvd_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1_overall <- vector()
  
  models_cvd_cvd_psm_1_1_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_cvd_cvd_psm_1_1_overall[[i]] <- coxph(formula(formula_freq),
                                                 data = group.cvd.dataset.matched %>%
                                                   mutate(intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1_overall <- rbind(predictions_cvd_cvd_stan_psm_1_1_overall, cbind(mean = models_cvd_cvd_psm_1_1_overall[[i]]$coefficients[1],
                                                                                                      lci = confint(models_cvd_cvd_psm_1_1_overall[[i]])[1,1],
                                                                                                      uci = confint(models_cvd_cvd_psm_1_1_overall[[i]])[1,2],
                                                                                                      intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_psm_1_1_overall <- predictions_cvd_cvd_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1_full <- vector()
  
  models_cvd_cvd_psm_1_1_full <- coxph(formula(formula_freq),
                                       data = group.cvd.dataset.matched)
  
  predictions_cvd_cvd_stan_psm_1_1_full <- rbind(predictions_cvd_cvd_stan_psm_1_1_full, cbind(mean = models_cvd_cvd_psm_1_1_full$coefficients[1],
                                                                                              lci = confint(models_cvd_cvd_psm_1_1_full)[1,1],
                                                                                              uci = confint(models_cvd_cvd_psm_1_1_full)[1,2]))
  
  predictions_cvd_cvd_stan_psm_1_1_full <- predictions_cvd_cvd_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_full.rds"))
  
  
}


## Propensity score matching + adjusted
if (class(try(
  
  # predictions for the CVD outcomes in population with no CVD strata sex and intervals
  predictions_cvd_cvd_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset.matched[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_cvd_cvd_psm_1_1_adjusted <- brms::brm(formula = formula(formula),
  #                                       data = group.cvd.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_cvd_cvd_psm_1_1_adjusted_male <- vector("list", quantiles)
  
  models_cvd_cvd_psm_1_1_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_cvd_cvd_psm_1_1_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                       data = group.cvd.dataset.matched %>%
                                                         mutate(sex = relevel(sex, ref = "Male"),
                                                                intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1_adjusted <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                        lci = confint(models_cvd_cvd_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                        uci = confint(models_cvd_cvd_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                        sex = "Male",
                                                                                                        intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    models_cvd_cvd_psm_1_1_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                         data = group.cvd.dataset.matched %>%
                                                           mutate(sex = relevel(sex, ref = "Female"),
                                                                  intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1_adjusted <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                        lci = confint(models_cvd_cvd_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                        uci = confint(models_cvd_cvd_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                        sex = "Female",
                                                                                                        intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted <- predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1_adjusted_overall <- vector()
  
  models_cvd_cvd_psm_1_1_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_cvd_cvd_psm_1_1_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                          data = group.cvd.dataset.matched %>%
                                                            mutate(intervals = relevel(intervals, ref = levels(group.cvd.dataset.matched$intervals)[i])))
    
    predictions_cvd_cvd_stan_psm_1_1_adjusted_overall <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted_overall, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_overall[[i]]$coefficients[1],
                                                                                                                        lci = confint(models_cvd_cvd_psm_1_1_adjusted_overall[[i]])[1,1],
                                                                                                                        uci = confint(models_cvd_cvd_psm_1_1_adjusted_overall[[i]])[1,2],
                                                                                                                        intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted_overall <- predictions_cvd_cvd_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted_overall.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full <- vector()
  
  models_cvd_cvd_psm_1_1_adjusted_full <- coxph(formula(formula_freq),
                                                data = group.cvd.dataset.matched)
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted_full, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_full$coefficients[1],
                                                                                                                lci = confint(models_cvd_cvd_psm_1_1_adjusted_full)[1,1],
                                                                                                                uci = confint(models_cvd_cvd_psm_1_1_adjusted_full)[1,2]))
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full <- predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted_full.rds"))
  
  
}


## Adjusted
# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.cvd.dataset)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.cvd.dataset,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.80, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("CVD outcomes for CVD population (n=", group.cvd.dataset%>%nrow(), ", event = ", group.cvd.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Overall population")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.cvd.dataset$intervals)
# 
# 
# kaplan_meier_cvd_cvd_adjusted <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                 labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.cvd.dataset,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.2),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("CVD outcomes for CVD population (n=", group.cvd.dataset%>%nrow(), ", event = ", group.cvd.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Overall population")
# 
# cumulative_incidence_cvd_cvd_adjusted <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                          labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_adjusted_effect <- coxph(formula(formula_freq_1),
                                                       data = group.cvd.dataset)

formula_freq_2 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_adjusted_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                                data = group.cvd.dataset)

formula_freq_3 <- "Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_cvd_cvd_spline_effects_adjusted_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                                data = group.cvd.dataset)

# anova(models_cvd_cvd_spline_effects_adjusted_effect, models_cvd_cvd_spline_effects_adjusted_spline_effect_3, models_cvd_cvd_spline_effects_adjusted_spline_effect_5)

dataset.used <- group.cvd.dataset %>%
  select(postdrug_mace_censtime_yrs, postdrug_mace_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.cvd.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.cvd.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD population") +
  labs(subtitle = paste0("Overall population (n=", group.cvd.dataset%>%nrow(), ", event = ", group.cvd.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.cvd.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_cvd_cvd_adjusted_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                 rel_heights = c(1,0.4), rel_widths = c(1,1))

m2 <- cph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.cvd.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.cvd.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("CVD outcomes for no CVD population") +
  labs(subtitle = paste0("Overall population (n=", group.cvd.dataset%>%nrow(), ", event = ", group.cvd.dataset%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"))


hist.dta <- group.cvd.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_cvd_cvd_adjusted_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                                 rel_heights = c(1,0.4), rel_widths = c(1,1))

#--- Cox  
if (class(try(
  
  # predictions for the CVD outcomes in population with no CVD strata sex and intervals
  predictions_cvd_cvd_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_mace_censtime_yrs | cens(postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_cvd_cvd_adjusted <- brms::brm(formula = formula(formula),
  #                                        data = group.cvd.dataset,
  #                                        family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_cvd_cvd_adjusted_male <- vector("list", quantiles)
  
  models_cvd_cvd_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_cvd_cvd_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                               data = group.cvd.dataset %>%
                                                 mutate(sex = relevel(sex, ref = "Male"),
                                                        intervals = relevel(intervals, ref = levels(group.cvd.dataset$intervals)[i])))
    
    predictions_cvd_cvd_stan_adjusted <- rbind(predictions_cvd_cvd_stan_adjusted, cbind(mean = models_cvd_cvd_adjusted_male[[i]]$coefficients[1],
                                                                                        lci = confint(models_cvd_cvd_adjusted_male[[i]])[1,1],
                                                                                        uci = confint(models_cvd_cvd_adjusted_male[[i]])[1,2],
                                                                                        sex = "Male",
                                                                                        intervals = levels(group.cvd.dataset$intervals)[i]))
    
    models_cvd_cvd_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                 data = group.cvd.dataset %>%
                                                   mutate(sex = relevel(sex, ref = "Female"),
                                                          intervals = relevel(intervals, ref = levels(group.cvd.dataset$intervals)[i])))
    
    predictions_cvd_cvd_stan_adjusted <- rbind(predictions_cvd_cvd_stan_adjusted, cbind(mean = models_cvd_cvd_adjusted_female[[i]]$coefficients[1],
                                                                                        lci = confint(models_cvd_cvd_adjusted_female[[i]])[1,1],
                                                                                        uci = confint(models_cvd_cvd_adjusted_female[[i]])[1,2],
                                                                                        sex = "Female",
                                                                                        intervals = levels(group.cvd.dataset$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_adjusted <- predictions_cvd_cvd_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_adjusted_overall <- vector()
  
  models_cvd_cvd_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_cvd_cvd_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                  data = group.cvd.dataset %>%
                                                    mutate(intervals = relevel(intervals, ref = levels(group.cvd.dataset$intervals)[i])))
    
    predictions_cvd_cvd_stan_adjusted_overall <- rbind(predictions_cvd_cvd_stan_adjusted_overall, cbind(mean = models_cvd_cvd_adjusted_overall[[i]]$coefficients[1],
                                                                                                        lci = confint(models_cvd_cvd_adjusted_overall[[i]])[1,1],
                                                                                                        uci = confint(models_cvd_cvd_adjusted_overall[[i]])[1,2],
                                                                                                        intervals = levels(group.cvd.dataset$intervals)[i]))
    
  }
  
  predictions_cvd_cvd_stan_adjusted_overall <- predictions_cvd_cvd_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted_overall.rds"))
  
}

if (class(try(
  
  # predictions for the CVD outcomes in the population with no CVD
  predictions_cvd_cvd_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.cvd.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_cvd_cvd_stan_adjusted_full <- vector()
  
  models_cvd_cvd_adjusted_full <- coxph(formula(formula_freq),
                                        data = group.cvd.dataset)
  
  predictions_cvd_cvd_stan_adjusted_full <- rbind(predictions_cvd_cvd_stan_adjusted_full, cbind(mean = models_cvd_cvd_adjusted_full$coefficients[1],
                                                                                                lci = confint(models_cvd_cvd_adjusted_full)[1,1],
                                                                                                uci = confint(models_cvd_cvd_adjusted_full)[1,2]))
  
  predictions_cvd_cvd_stan_adjusted_full <- predictions_cvd_cvd_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_cvd_cvd_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted_full.rds"))
  
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

# cvd_cvd_overall_axis_min <- plyr::round_any(floor(min(c(predictions_cvd_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_cvd_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_cvd_cvd_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_cvd_cvd_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                           predictions_cvd_cvd_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# cvd_cvd_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_cvd_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_cvd_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                           predictions_cvd_cvd_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                             predictions_cvd_cvd_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                             predictions_cvd_cvd_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)
# 
# 
# cvd_cvd_strata_axis_min <- plyr::round_any(floor(min(c(predictions_cvd_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_cvd_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_cvd_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_cvd_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                          predictions_cvd_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# cvd_cvd_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_cvd_cvd_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_cvd_cvd_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_cvd_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                            predictions_cvd_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                            predictions_cvd_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)



# Propensity score matching
plot_cvd_cvd_psm_1_1_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_overall %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.cvd.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_psm_1_1_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.cvd.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_psm_1_1_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.cvd.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Propensity score matching + adjusted
plot_cvd_cvd_psm_1_1_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.cvd.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_psm_1_1_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.cvd.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_psm_1_1_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset.matched%>%filter(intervals==levels(group.cvd.dataset.matched$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.cvd.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Adjusted
plot_cvd_cvd_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted_overall %>% slice(4:6),
  predictions_cvd_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.cvd.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.cvd.dataset%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_cvd_cvd_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_cvd_cvd_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_cvd_cvd_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.cvd.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[1])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.cvd.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[2])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.cvd.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[3])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.cvd.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[4])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.cvd.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[5])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.cvd.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.cvd.dataset%>%filter(intervals==levels(group.cvd.dataset$intervals)[6])%>%filter(postdrug_mace_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.cvd.dataset%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/cvd_cvd_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("CVD outcomes for no CVD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_cvd_cvd_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_cvd_cvd_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_cvd_cvd_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/cvd_cvd_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("CVD outcomes for no CVD population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_cvd_cvd_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_cvd_cvd_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_cvd_cvd_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_cvd_cvd_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_cvd_cvd_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_cvd_cvd_adjusted_male
upViewport()

dev.off()




#:--------------------------
# Heart failure outcome

patient_prop_scores_qrisk <- readRDS(paste0(output_path, "/additional_outcomes/patient_prop_scores_qrisk.rds"))

hf.dataset <- set_up_data_sglt2_glp1(dataset.type="hf.dataset") %>%
  left_join(patient_prop_scores_qrisk, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

group.hf.dataset <- group_values(data = hf.dataset,
                                 variable = "effects",
                                 breaks = interval_breaks) %>%
  drop_na(intervals)

matching_hf <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
  data = group.hf.dataset,
  method = "nearest",
  distance = group.hf.dataset[,"prop.score"],
  replace = FALSE,
  m.order = "largest",
  caliper = 0.05,
  mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
)

# require(cobalt)
# cobalt::love.plot(matching_hf, binary = "std", thresholds = c(m = .1), sample.names = c("Unmatched", "Matched"), title = "Full dataset")

n_drug <- 2*sum(!is.na(matching_hf$match.matrix))

group.hf.dataset.matched <- group.hf.dataset %>%
  slice(which(matching_hf$weights == 1))


# Heart failure survival analysis

## Propensity score matching

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.hf.dataset.matched)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.hf.dataset.matched,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.70, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("HF outcomes for HF population (n=", group.hf.dataset.matched%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Propensity score matching")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.hf.dataset.matched$intervals)
# 
# 
# kaplan_meier_hf_hf_psm_1_1 <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                  labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.hf.dataset.matched,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.3),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("HF outcomes for HF population (n=", group.hf.dataset.matched%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Propensity score matching")
# 
# cumulative_incidence_hf_hf_psm_1_1 <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                       labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_psm_1_1_effect <- coxph(formula(formula_freq_1),
                                                    data = group.hf.dataset.matched)

formula_freq_2 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_psm_1_1_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                             data = group.hf.dataset.matched)

formula_freq_3 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_psm_1_1_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                             data = group.hf.dataset.matched)

# anova(models_hf_hf_spline_effects_psm_1_1_effect, models_hf_hf_spline_effects_psm_1_1_spline_effect_3, models_hf_hf_spline_effects_psm_1_1_spline_effect_5)


dataset.used <- group.hf.dataset.matched %>%
  select(postdrug_hf_censtime_yrs, postdrug_hf_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.hf.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.hf.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no HF population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.hf.dataset.matched%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.hf.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_hf_hf_psm_1_1_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                              rel_heights = c(1,0.4), rel_widths = c(1,1))

m2 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.hf.dataset.matched$effects, .01, na.rm=TRUE)
c99 <- quantile(group.hf.dataset.matched$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no HF population") +
  labs(subtitle = paste0("Propensity score matching (n=", group.hf.dataset.matched%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.hf.dataset.matched %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_hf_hf_psm_1_1_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                              rel_heights = c(1,0.4), rel_widths = c(1,1))

#--- Cox
if (class(try(
  
  # predictions for the HF outcomes in population with no HF strata sex and intervals
  predictions_hf_hf_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1 <- vector()
  
  # formula <- 'postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)'
  # 
  # models_hf_hf_psm_1_1 <- brms::brm(formula = formula(formula),
  #                                       data = group.hf.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3)"
  
  models_hf_hf_psm_1_1_male <- vector("list", quantiles)
  
  models_hf_hf_psm_1_1_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_psm_1_1_male[[i]] <- coxph(formula(formula_freq),
                                            data = group.hf.dataset.matched %>%
                                              mutate(sex = relevel(sex, ref = "Male"),
                                                     intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1 <- rbind(predictions_hf_hf_stan_psm_1_1, cbind(mean = models_hf_hf_psm_1_1_male[[i]]$coefficients[1],
                                                                                  lci = confint(models_hf_hf_psm_1_1_male[[i]])[1,1],
                                                                                  uci = confint(models_hf_hf_psm_1_1_male[[i]])[1,2],
                                                                                  sex = "Male",
                                                                                  intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
    models_hf_hf_psm_1_1_female[[i]] <- coxph(formula(formula_freq),
                                              data = group.hf.dataset.matched %>%
                                                mutate(sex = relevel(sex, ref = "Female"),
                                                       intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1 <- rbind(predictions_hf_hf_stan_psm_1_1, cbind(mean = models_hf_hf_psm_1_1_female[[i]]$coefficients[1],
                                                                                  lci = confint(models_hf_hf_psm_1_1_female[[i]])[1,1],
                                                                                  uci = confint(models_hf_hf_psm_1_1_female[[i]])[1,2],
                                                                                  sex = "Female",
                                                                                  intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_psm_1_1 <- predictions_hf_hf_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_psm_1_1_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1_overall <- vector()
  
  models_hf_hf_psm_1_1_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_psm_1_1_overall[[i]] <- coxph(formula(formula_freq),
                                               data = group.hf.dataset.matched %>%
                                                 mutate(intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1_overall <- rbind(predictions_hf_hf_stan_psm_1_1_overall, cbind(mean = models_hf_hf_psm_1_1_overall[[i]]$coefficients[1],
                                                                                                  lci = confint(models_hf_hf_psm_1_1_overall[[i]])[1,1],
                                                                                                  uci = confint(models_hf_hf_psm_1_1_overall[[i]])[1,2],
                                                                                                  intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_psm_1_1_overall <- predictions_hf_hf_stan_psm_1_1_overall %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1_overall, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_overall.rds"))
  
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_psm_1_1_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  formula_freq <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + qrisk2_10yr_score"
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1_full <- vector()
  
  models_hf_hf_psm_1_1_full <- coxph(formula(formula_freq),
                                     data = group.hf.dataset.matched)
  
  predictions_hf_hf_stan_psm_1_1_full <- rbind(predictions_hf_hf_stan_psm_1_1_full, cbind(mean = models_hf_hf_psm_1_1_full$coefficients[1],
                                                                                          lci = confint(models_hf_hf_psm_1_1_full)[1,1],
                                                                                          uci = confint(models_hf_hf_psm_1_1_full)[1,2]))
  
  predictions_hf_hf_stan_psm_1_1_full <- predictions_hf_hf_stan_psm_1_1_full %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1_full, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_full.rds"))
  
  
}

## Propensity score matching + adjusted
if (class(try(
  
  # predictions for the HF outcomes in population with no HF strata sex and intervals
  predictions_hf_hf_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset.matched[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_hf_hf_psm_1_1_adjusted <- brms::brm(formula = formula(formula),
  #                                       data = group.hf.dataset.matched,
  #                                       family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_hf_hf_psm_1_1_adjusted_male <- vector("list", quantiles)
  
  models_hf_hf_psm_1_1_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_psm_1_1_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                                     data = group.hf.dataset.matched %>%
                                                       mutate(sex = relevel(sex, ref = "Male"),
                                                              intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1_adjusted <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted, cbind(mean = models_hf_hf_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                    lci = confint(models_hf_hf_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                    uci = confint(models_hf_hf_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                    sex = "Male",
                                                                                                    intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
    models_hf_hf_psm_1_1_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                                       data = group.hf.dataset.matched %>%
                                                         mutate(sex = relevel(sex, ref = "Female"),
                                                                intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1_adjusted <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted, cbind(mean = models_hf_hf_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                    lci = confint(models_hf_hf_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                    uci = confint(models_hf_hf_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                    sex = "Female",
                                                                                                    intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_psm_1_1_adjusted <- predictions_hf_hf_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_psm_1_1_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1_adjusted_overall <- vector()
  
  models_hf_hf_psm_1_1_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_psm_1_1_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                        data = group.hf.dataset.matched %>%
                                                          mutate(intervals = relevel(intervals, ref = levels(group.hf.dataset.matched$intervals)[i])))
    
    predictions_hf_hf_stan_psm_1_1_adjusted_overall <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted_overall, cbind(mean = models_hf_hf_psm_1_1_adjusted_overall[[i]]$coefficients[1],
                                                                                                                    lci = confint(models_hf_hf_psm_1_1_adjusted_overall[[i]])[1,1],
                                                                                                                    uci = confint(models_hf_hf_psm_1_1_adjusted_overall[[i]])[1,2],
                                                                                                                    intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_psm_1_1_adjusted_overall <- predictions_hf_hf_stan_psm_1_1_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted_overall.rds"))
  
}


if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_psm_1_1_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset.matched[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_psm_1_1_adjusted_full <- vector()
  
  models_hf_hf_psm_1_1_adjusted_full <- coxph(formula(formula_freq),
                                              data = group.hf.dataset.matched)
  
  predictions_hf_hf_stan_psm_1_1_adjusted_full <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted_full, cbind(mean = models_hf_hf_psm_1_1_adjusted_full$coefficients[1],
                                                                                                            lci = confint(models_hf_hf_psm_1_1_adjusted_full)[1,1],
                                                                                                            uci = confint(models_hf_hf_psm_1_1_adjusted_full)[1,2]))
  
  predictions_hf_hf_stan_psm_1_1_adjusted_full <- predictions_hf_hf_stan_psm_1_1_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_psm_1_1_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted_full.rds"))
  
  
}


## Adjusted

# #--- Kaplan-meier curve
# fit_kaplan <- survfit(formula("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + sex + intervals"),
#                       data = group.hf.dataset)
# 
# ggsurvfitplot <- ggsurvplot(fit_kaplan,
#                             data = group.hf.dataset.matched,
#                             palette = c("dodgerblue2", "#f1a340"),
#                             surv.median.line = "none",
#                             ylim = c(0.70, 1),
#                             censor = FALSE,
#                             conf.int = TRUE,
#                             legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                             legend.title = "Therapy",
#                             legend = "bottom",
#                             title = paste0("HF outcomes for HF population (n=", group.hf.dataset%>%nrow(), ", event = ", group.hf.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                             subtitle = "Overall population")
# 
# interval.labs <- c(paste0("Predicted HbA1c benefit on SGLT2i >5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on SGLT2i 0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA 3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                    paste0("Predicted HbA1c benefit on GLP1-RA >5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))
# 
# names(interval.labs) <- levels(group.hf.dataset$intervals)
# 
# 
# kaplan_meier_hf_hf_adjusted <- ggsurvfitplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                               labeller = labeller(intervals = interval.labs))
# 
# #---- Cumulative incidence
# ggcumincplot <- ggsurvplot(fit_kaplan,
#                            data = group.hf.dataset.matched,
#                            fun = "event",
#                            palette = c("dodgerblue2", "#f1a340"),
#                            surv.median.line = "none",
#                            ylim = c(0, 0.3),
#                            censor = FALSE,
#                            conf.int = TRUE,
#                            legend.labs = c(rep(c("GLP1-RA", "SGLT2i"), each = 12)),
#                            legend.title = "Therapy",
#                            legend = "bottom",
#                            title = paste0("HF outcomes for HF population (n=", group.hf.dataset%>%nrow(), ", event = ", group.hf.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
#                            subtitle = "Overall population")
# 
# cumulative_incidence_hf_hf_adjusted <- ggcumincplot$plot + facet_wrap(intervals~sex, ncol = 2,
#                                                                        labeller = labeller(intervals = interval.labs))

#--- Continuous drugclass*spline(hba1c.diff)

formula_freq_1 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*effects + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_adjusted_effect <- coxph(formula(formula_freq_1),
                                                     data = group.hf.dataset)

formula_freq_2 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,3) + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_adjusted_spline_effect_3 <- coxph(formula(formula_freq_2),
                                                              data = group.hf.dataset)

formula_freq_3 <- "Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + factor(drugclass)*rcs(effects,5) + rcs(qrisk2_10yr_score, 3)"

models_hf_hf_spline_effects_adjusted_spline_effect_5 <- coxph(formula(formula_freq_3),
                                                              data = group.hf.dataset)

# anova(models_hf_hf_spline_effects_adjusted_effect, models_hf_hf_spline_effects_adjusted_spline_effect_3, models_hf_hf_spline_effects_adjusted_spline_effect_5)


dataset.used <- group.hf.dataset %>%
  select(postdrug_hf_censtime_yrs, postdrug_hf_censvar, drugclass, effects, qrisk2_10yr_score)
ddist <- datadist(dataset.used); options(datadist='ddist')

m1 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 5) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.hf.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.hf.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m1,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no HF population") +
  labs(subtitle = paste0("Overall population (n=", group.hf.dataset%>%nrow(), ", event = ", group.hf.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.hf.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_hf_hf_adjusted_spline_5 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                               rel_heights = c(1,0.4), rel_widths = c(1,1))


m2 <- cph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass*rcs(effects, 3) + rcs(qrisk2_10yr_score, 3),
          data = dataset.used,x=T,y=T)

c1 <- quantile(group.hf.dataset$effects, .01, na.rm=TRUE)
c99 <- quantile(group.hf.dataset$effects, .99, na.rm=TRUE)

# run the contrast calculation BMI by casevalid
contrast_spline.1 <- rms::contrast(m2,list(drugclass = "SGLT2", effects = seq(c1,c99,by=0.05)),list(drugclass = "GLP1", effects = seq(c1,c99,by=0.05)))
# save the contrast calculations in a dataframe
contrast_spline_df <- as.data.frame(contrast_spline.1[c('effects','Contrast','Lower','Upper')])

contrast_spline_plot_1 <- ggplot(data=contrast_spline_df,aes(x=effects, y=exp(Contrast))) +
  geom_line(data=contrast_spline_df,aes(x=effects, y=exp(Contrast)), size=1) +
  xlab(expression("SGLT2i HbA1c benefit")) +
  ylab("HR") +
  # scale_x_continuous(breaks = seq(0,50,5)) +
  #scale_y_continuous(breaks = seq(0.8,1.6,0.1), limits = c(0.8,1.6)) +
  geom_ribbon(data=contrast_spline_df,aes(x=effects,ymin=exp(Lower),ymax=exp(Upper)),alpha=0.5) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  theme(legend.position=c(0.8, 0.1)) +
  theme(legend.title = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.line = element_line(colour =  "grey50" ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.text = element_text(colour="black", size=rel(1))) +
  ggtitle("HF outcomes for no HF population") +
  labs(subtitle = paste0("Overall population (n=", group.hf.dataset%>%nrow(), ", event = ", group.hf.dataset%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"))


hist.dta <- group.hf.dataset %>% filter(effects>=c1 &  effects <= c99)
hist.dta$dummy <- 1
x_hist <- marginal_distribution(hist.dta, "effects")


# Arranging the plot using cowplot
continuous_effects_hazard_hf_hf_adjusted_spline_3 <- plot_grid(contrast_spline_plot_1, x_hist, ncol = 1,align = 'hv',
                                                               rel_heights = c(1,0.4), rel_widths = c(1,1))

#--- Cox
if (class(try(
  
  # predictions for the HF outcomes in population with no HF strata sex and intervals
  predictions_hf_hf_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_adjusted <- vector()
  
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset[,breakdown_adjust], is.factor)
  
  # formula <- paste0("postdrug_hf_censtime_yrs | cens(postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  # 
  # models_hf_hf_adjusted <- brms::brm(formula = formula(formula),
  #                                        data = group.hf.dataset,
  #                                        family = "weibull")
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex*factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  models_hf_hf_adjusted_male <- vector("list", quantiles)
  
  models_hf_hf_adjusted_female <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_adjusted_male[[i]] <- coxph(formula(formula_freq),
                                             data = group.hf.dataset %>%
                                               mutate(sex = relevel(sex, ref = "Male"),
                                                      intervals = relevel(intervals, ref = levels(group.hf.dataset$intervals)[i])))
    
    predictions_hf_hf_stan_adjusted <- rbind(predictions_hf_hf_stan_adjusted, cbind(mean = models_hf_hf_adjusted_male[[i]]$coefficients[1],
                                                                                    lci = confint(models_hf_hf_adjusted_male[[i]])[1,1],
                                                                                    uci = confint(models_hf_hf_adjusted_male[[i]])[1,2],
                                                                                    sex = "Male",
                                                                                    intervals = levels(group.hf.dataset$intervals)[i]))
    
    models_hf_hf_adjusted_female[[i]] <- coxph(formula(formula_freq),
                                               data = group.hf.dataset %>%
                                                 mutate(sex = relevel(sex, ref = "Female"),
                                                        intervals = relevel(intervals, ref = levels(group.hf.dataset$intervals)[i])))
    
    predictions_hf_hf_stan_adjusted <- rbind(predictions_hf_hf_stan_adjusted, cbind(mean = models_hf_hf_adjusted_female[[i]]$coefficients[1],
                                                                                    lci = confint(models_hf_hf_adjusted_female[[i]])[1,1],
                                                                                    uci = confint(models_hf_hf_adjusted_female[[i]])[1,2],
                                                                                    sex = "Female",
                                                                                    intervals = levels(group.hf.dataset$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_adjusted <- predictions_hf_hf_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted.rds"))
  
}

if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_adjusted_overall <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted_overall.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_adjusted_overall <- vector()
  
  models_hf_hf_adjusted_overall <- vector("list", quantiles)
  
  for (i in mnumber) {
    
    models_hf_hf_adjusted_overall[[i]] <- coxph(formula(formula_freq),
                                                data = group.hf.dataset %>%
                                                  mutate(intervals = relevel(intervals, ref = levels(group.hf.dataset$intervals)[i])))
    
    predictions_hf_hf_stan_adjusted_overall <- rbind(predictions_hf_hf_stan_adjusted_overall, cbind(mean = models_hf_hf_adjusted_overall[[i]]$coefficients[1],
                                                                                                    lci = confint(models_hf_hf_adjusted_overall[[i]])[1,1],
                                                                                                    uci = confint(models_hf_hf_adjusted_overall[[i]])[1,2],
                                                                                                    intervals = levels(group.hf.dataset$intervals)[i]))
    
  }
  
  predictions_hf_hf_stan_adjusted_overall <- predictions_hf_hf_stan_adjusted_overall %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_adjusted_overall, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted_overall.rds"))
  
}


if (class(try(
  
  # predictions for the HF outcomes in the population with no HF
  predictions_hf_hf_stan_adjusted_full <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted_full.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  breakdown_adjust <- unique(c(variables_mu, variables_tau))
  # categorical variables in breakdown
  factors <- sapply(group.hf.dataset[,breakdown_adjust], is.factor)
  
  formula_freq <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ factor(drugclass) + rcs(qrisk2_10yr_score, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))
  
  # maximum number of deciles being tested
  quantiles <- length(levels(group.hf.dataset[,"intervals"]))
  # create lists with results
  mnumber = c(1:quantiles)
  predictions_hf_hf_stan_adjusted_full <- vector()
  
  models_hf_hf_adjusted_full <- coxph(formula(formula_freq),
                                      data = group.hf.dataset)
  
  predictions_hf_hf_stan_adjusted_full <- rbind(predictions_hf_hf_stan_adjusted_full, cbind(mean = models_hf_hf_adjusted_full$coefficients[1],
                                                                                            lci = confint(models_hf_hf_adjusted_full)[1,1],
                                                                                            uci = confint(models_hf_hf_adjusted_full)[1,2]))
  
  predictions_hf_hf_stan_adjusted_full <- predictions_hf_hf_stan_adjusted_full %>%
    as.data.frame()
  
  saveRDS(predictions_hf_hf_stan_adjusted_full, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted_full.rds"))
  
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

# hf_hf_overall_axis_min <- plyr::round_any(floor(min(c(predictions_hf_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                       predictions_hf_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                       predictions_hf_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                       predictions_hf_hf_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                          predictions_hf_hf_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                          predictions_hf_hf_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# hf_hf_overall_axis_max <- plyr::round_any(ceiling(max(c(predictions_hf_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_hf_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_hf_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_hf_hf_stan_psm_1_1_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                            predictions_hf_hf_stan_psm_1_1_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                            predictions_hf_hf_stan_adjusted_overall %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)
# 
# 
# hf_hf_strata_axis_min <- plyr::round_any(floor(min(c(predictions_hf_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                      predictions_hf_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                      predictions_hf_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                      predictions_hf_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                         predictions_hf_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(),
#                                                         predictions_hf_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp()))), 2, f = floor)
# 
# hf_hf_strata_axis_max <- plyr::round_any(ceiling(max(c(predictions_hf_hf_stan_psm_1_1_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_hf_hf_stan_psm_1_1_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_hf_hf_stan_adjusted_full %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min() %>% exp(), 
#                                                        predictions_hf_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(), 
#                                                           predictions_hf_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp(),
#                                                           predictions_hf_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() %>% exp()))), 2, f = ceiling)


# Propensity score matching
plot_hf_hf_psm_1_1_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_overall %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.hf.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_psm_1_1_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1 %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.hf.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_psm_1_1_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1 %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.hf.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Propensity score matching + adjusted
plot_hf_hf_psm_1_1_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted_overall %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.hf.dataset.matched), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_psm_1_1_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.hf.dataset.matched%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_psm_1_1_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_psm_1_1_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_psm_1_1_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset.matched$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset.matched$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset.matched$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset.matched$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset.matched$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset.matched$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset.matched%>%filter(intervals==levels(group.hf.dataset.matched$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Propensity score matching + adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.hf.dataset.matched%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

# Adjusted
plot_hf_hf_adjusted_overall <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted_overall %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted_overall %>% slice(4:6),
  predictions_hf_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Overall population (n=", nrow(group.hf.dataset), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_adjusted_male<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted %>% filter(sex == "Male") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(sex=="Male")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Male")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Male (n=", nrow(group.hf.dataset%>%filter(sex=="Male")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))

plot_hf_hf_adjusted_female<- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(1:3),
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", mean = NA, lci = NA, uci = NA),
  predictions_hf_hf_stan_adjusted %>% filter(sex == "Female") %>% select(-sex) %>% slice(4:6),
  predictions_hf_hf_stan_adjusted_full %>% cbind(intervals = "Average treatment effect")
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == levels(group.hf.dataset$intervals)[1], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[1])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                            ifelse(intervals == levels(group.hf.dataset$intervals)[2], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[2])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                   ifelse(intervals == levels(group.hf.dataset$intervals)[3], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[3])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                          ifelse(intervals == levels(group.hf.dataset$intervals)[4], paste0("0-3 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[4])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                 ifelse(intervals == levels(group.hf.dataset$intervals)[5], paste0("3-5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[5])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"),
                                                        ifelse(intervals == levels(group.hf.dataset$intervals)[6], paste0(">5 mmol/mol (n=", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(sex=="Female")%>%nrow(), ", event = ", group.hf.dataset%>%filter(intervals==levels(group.hf.dataset$intervals)[6])%>%filter(postdrug_hf_censvar==1)%>%filter(sex=="Female")%>%nrow(), ")"), intervals))))))) %>%
  rename("lower" = "lci", "upper" = "uci") %>%
  mutate(mean = exp(mean),
         lower = exp(lower),
         upper = exp(upper)) %>%
  forestplot(labeltext = intervals,
             fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawNormalCI, fpDrawDiamondCI),
             ci.vertices = TRUE,
             zero = 1,
             title = "Adjusted",
             xlog = TRUE,
             ci.vertices.height = 0.05,
             boxsize = .1,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Hazard Ratio") %>%
  fp_add_header(paste0("Female (n=", nrow(group.hf.dataset%>%filter(sex=="Female")), ")")) %>%
  fp_add_lines(h_2 = "black",
               h_6 = gpar(col = "black", lty = 2),
               h_10 = gpar(col = "black", lty = 2))


pdf(width = 7, height = 12, "Plots/hf_hf_overall.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 1, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HF outcomes for no HF population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_hf_hf_psm_1_1_overall
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_hf_hf_psm_1_1_adjusted_overall
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_hf_hf_adjusted_overall
upViewport()

dev.off()

pdf(width = 14, height = 12, "Plots/hf_hf_strata.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HF outcomes for no HF population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_hf_hf_psm_1_1_female
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_hf_hf_psm_1_1_male
upViewport()
# third plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 1))
plot_hf_hf_psm_1_1_adjusted_female
upViewport()
# forth plot
pushViewport(viewport(layout.pos.row = 3,
                      layout.pos.col = 2))
plot_hf_hf_psm_1_1_adjusted_male
upViewport()
# fifth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 1))
plot_hf_hf_adjusted_female
upViewport()
# sixth plot
pushViewport(viewport(layout.pos.row = 4,
                      layout.pos.col = 2))
plot_hf_hf_adjusted_male
upViewport()

dev.off()


# #### Kaplan-meier curves
# pdf(width = 12, height = 14, "Plots/11.06.kaplan_meier_plots.pdf")
# 
# # No comorbidities, CVD outcomes, Propensity score matching
# kaplan_meier_no_co_cvd_psm_1_1
# # No comorbidities, CVD outcomes, Overall population
# kaplan_meier_no_co_cvd_adjusted
# # No comorbidities, HF outcomes, Propensity score matching
# kaplan_meier_no_co_hf_psm_1_1
# # No comorbidities, HF outcomes, Overall population
# kaplan_meier_no_co_hf_adjusted
# # No comorbidities, CVD outcomes, Propensity score matching
# kaplan_meier_cvd_cvd_psm_1_1
# # No comorbidities, CVD outcomes, Overall population
# kaplan_meier_cvd_cvd_adjusted
# # No comorbidities, HF outcomes, Propensity score matching
# kaplan_meier_hf_hf_psm_1_1
# # No comorbidities, HF outcomes, Overall population
# kaplan_meier_hf_hf_adjusted
# 
# dev.off()


# #### Kaplan-meier cumulative incidence curves
# pdf(width = 12, height = 14, "Plots/11.06.cumulative_incidence_plots.pdf")
# 
# # No comorbidities, CVD outcomes, Propensity score matching
# cumulative_incidence_no_co_cvd_psm_1_1
# # No comorbidities, CVD outcomes, Overall population
# cumulative_incidence_no_co_cvd_adjusted
# # No comorbidities, HF outcomes, Propensity score matching
# cumulative_incidence_no_co_hf_psm_1_1
# # No comorbidities, HF outcomes, Overall population
# cumulative_incidence_no_co_hf_adjusted
# # No comorbidities, CVD outcomes, Propensity score matching
# cumulative_incidence_cvd_cvd_psm_1_1
# # No comorbidities, CVD outcomes, Overall population
# cumulative_incidence_cvd_cvd_adjusted
# # No comorbidities, HF outcomes, Propensity score matching
# cumulative_incidence_hf_hf_psm_1_1
# # No comorbidities, HF outcomes, Overall population
# cumulative_incidence_hf_hf_adjusted
# 
# dev.off()

### Continuous effects hazard plot
pdf(width = 9, height = 9, "Plots/11.06.continuous_hazard_plots_spline_5.pdf")
# No comorbidities, CVD outcomes, PSM
continuous_effects_hazard_no_co_cvd_psm_1_1_spline_5
# No comorbidities, CVD outcomes, adjusted
continuous_effects_hazard_no_co_cvd_adjusted_spline_5
# No comorbidities, HF outcomes, PSM
continuous_effects_hazard_no_co_hf_psm_1_1_spline_5
# No comorbidities, HF outcomes, adjusted
continuous_effects_hazard_no_co_hf_adjusted_spline_5
# No CVD, CVD outcomes, PSM
continuous_effects_hazard_cvd_cvd_psm_1_1_spline_5
# No CVD, CVD outcomes, adjusted
continuous_effects_hazard_cvd_cvd_adjusted_spline_5
# No HF, HF outcomes, PSM
continuous_effects_hazard_hf_hf_psm_1_1_spline_5
# No HF, HF outcomes, adjusted
continuous_effects_hazard_hf_hf_adjusted_spline_5

dev.off()

pdf(width = 9, height = 9, "Plots/11.06.continuous_hazard_plots_spline_3.pdf")
# No comorbidities, CVD outcomes, PSM
continuous_effects_hazard_no_co_cvd_psm_1_1_spline_3
# No comorbidities, CVD outcomes, adjusted
continuous_effects_hazard_no_co_cvd_adjusted_spline_3
# No comorbidities, HF outcomes, PSM
continuous_effects_hazard_no_co_hf_psm_1_1_spline_3
# No comorbidities, HF outcomes, adjusted
continuous_effects_hazard_no_co_hf_adjusted_spline_3
# No CVD, CVD outcomes, PSM
continuous_effects_hazard_cvd_cvd_psm_1_1_spline_3
# No CVD, CVD outcomes, adjusted
continuous_effects_hazard_cvd_cvd_adjusted_spline_3
# No HF, HF outcomes, PSM
continuous_effects_hazard_hf_hf_psm_1_1_spline_3
# No HF, HF outcomes, adjusted
continuous_effects_hazard_hf_hf_adjusted_spline_3

dev.off()



#### Merge PDF's together
# Overall population
qpdf::pdf_combine(input = c("Plots/hba1c_grouping_overall.pdf",
                            "Plots/weight_overall.pdf",
                            "Plots/egfr_overall.pdf",
                            "Plots/discontinuation_overall.pdf",
                            "Plots/no_co_cvd_overall.pdf",
                            "Plots/no_co_hf_overall.pdf",
                            "Plots/cvd_cvd_overall.pdf",
                            "Plots/hf_hf_overall.pdf"),
                  output = "Plots/11.06.overall_additional_outcomes.pdf")

# Strata population
qpdf::pdf_combine(input = c("Plots/hba1c_grouping_strata.pdf",
                            "Plots/weight_strata.pdf",
                            "Plots/egfr_strata.pdf",
                            "Plots/discontinuation_strata.pdf",
                            "Plots/no_co_cvd_strata.pdf",
                            "Plots/no_co_hf_strata.pdf",
                            "Plots/cvd_cvd_strata.pdf",
                            "Plots/hf_hf_strata.pdf"),
                  output = "Plots/11.06.strata_additional_outcomes.pdf")


file.remove(c("Plots/hba1c_grouping_overall.pdf",
              "Plots/weight_overall.pdf",
              "Plots/egfr_overall.pdf",
              "Plots/discontinuation_overall.pdf",
              "Plots/no_co_cvd_overall.pdf",
              "Plots/no_co_hf_overall.pdf",
              "Plots/cvd_cvd_overall.pdf",
              "Plots/hf_hf_overall.pdf"))

file.remove(c("Plots/hba1c_grouping_strata.pdf",
              "Plots/weight_strata.pdf",
              "Plots/egfr_strata.pdf",
              "Plots/discontinuation_strata.pdf",
              "Plots/no_co_cvd_strata.pdf",
              "Plots/no_co_hf_strata.pdf",
              "Plots/cvd_cvd_strata.pdf",
              "Plots/hf_hf_strata.pdf"))


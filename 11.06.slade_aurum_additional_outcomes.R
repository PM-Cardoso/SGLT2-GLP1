####################
## Description:
##  - In this file we check:
##    - Weight change
##    - eGFR change
##    - Discontinuation
####################

library(tidyverse)
library(margins)
library(forestplot)
library(grid)
library(rms)

## Analysis of outcomes
# library(rstanarm)
# library(bcf)


## increase memery usage to 50gb of RAM
# options(java.parameters = "-Xmx100g")
# 
# library(bartMachine)



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
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))

breaks_hba1c <- c(-5, -3, 0, 3, 5)

group.hba1c.dataset <- group_values(data = hba1c.train,
                                    variable = "effects",
                                    breaks = breaks_hba1c) %>%
  drop_na(intervals) %>%
  rename("hba1c_diff" = "effects")

group.hba1c.dataset.male <- group.hba1c.dataset %>% filter(sex == "Male")
group.hba1c.dataset.female <- group.hba1c.dataset %>% filter(sex == "Female")

ATE_hba1c_male <- calc_ATE_validation_adjust(group.hba1c.dataset.male, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

ATE_hba1c_female <- calc_ATE_validation_adjust(group.hba1c.dataset.female, "posthba1cfinal", quantile_var = "intervals", breakdown = model_variables, adjust = TRUE)

hba1c_axis_min <- plyr::round_any(floor(min(c(ATE_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% min(), ATE_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% min()))), 2, f = floor)

hba1c_axis_max <- plyr::round_any(ceiling(max(c(ATE_hba1c_male[["effects"]] %>% select(c("obs","lci","uci")) %>% max(), ATE_hba1c_female[["effects"]] %>% select(c("obs","lci","uci")) %>% max()))), 2, f = ceiling)



plot_hba1c_male <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", N = NA, hba1c_diff.pred = NA, obs = NA, lci = NA, uci = NA),
  ATE_hba1c_male[["effects"]][1:3,],
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", N = NA, hba1c_diff.pred = NA, obs = NA, lci = NA, uci = NA),
  ATE_hba1c_male[["effects"]][4:6,]
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1602)", ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=2315)", ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=6474)", ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5285)", ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1084)", ifelse(intervals == "(5,31]", ">5 mmol/mol (n=569)", intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             xticks = seq(hba1c_axis_min, hba1c_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header("Male (n=17329")
  
plot_hba1c_female <- rbind(
  cbind(intervals = "Predicted HbA1c benefit on SGLT2i", N = NA, hba1c_diff.pred = NA, obs = NA, lci = NA, uci = NA),
  ATE_hba1c_female[["effects"]][1:3,],
  cbind(intervals = "Predicted HbA1c benefit on GLP1-RA", N = NA, hba1c_diff.pred = NA, obs = NA, lci = NA, uci = NA),
  ATE_hba1c_female[["effects"]][4:6,]
) %>%
  as.data.frame() %>%
  mutate(obs = as.numeric(obs),
         lci = as.numeric(lci),
         uci = as.numeric(uci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=591)", ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=538)", ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2267)", ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=3779)", ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=2694)", ifelse(intervals == "(5,31]", ">5 mmol/mol (n=1487)", intervals))))))) %>%
  rename("mean" = "obs", "lower" = "lci", "upper" = "uci") %>%
  forestplot(labeltext = intervals,
             ci.vertices = TRUE,
             xticks = seq(hba1c_axis_min, hba1c_axis_max, 2),
             ci.vertices.height = 0.05,
             boxsize = .1,
             xlab = "Average treatment effect (mmol/mol)") %>%
  fp_add_header("Female (n=11356)")


pdf(width = 14, height = 3.5, "Plots/hba1c_grouping.pdf")

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           ncol = 2, heights = unit(c(0.5, 5), "null"))))
# title
grid.text("HbA1c treatment effect", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

# first plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 1))
plot_hba1c_male
upViewport()
# second plot
pushViewport(viewport(layout.pos.row = 2,
                      layout.pos.col = 2))
plot_hba1c_female

dev.off()

#:--------------------------------------------------------
# Weight change

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


breaks_weight <- c(-5, -3, 0, 3, 5)

group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = breaks_weight) %>%
  drop_na(intervals)


## Fit a bcf for the weight change
# if (class(try(
#   
#   weight_bcf <- readRDS(paste0(output_path, "/additional_outcomes/weight_bcf.rds"))
#   
#   , silent = TRUE)) == "try-error") {
#   
#   weight_bcf <- bcf::bcf(y = group.weight.dataset$w.change,
#                          z = group.weight.dataset %>%
#                            select(drugclass) %>%
#                            mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
#                            unlist(),
#                          x_control = group.weight.dataset %>%
#                            select(
#                              all_of(c(variables_mu, variables_tau)), "intervals", "preweight"
#                            ) %>%
#                            mutate_all(funs(as.numeric(.))) %>%
#                            as.matrix(),
#                          x_moderate = group.weight.dataset %>%
#                            select(
#                              all_of(c(variables_mu, variables_tau)), "intervals", "preweight"
#                            ) %>%
#                            mutate_all(funs(as.numeric(.))) %>%
#                            as.matrix(),
#                          pihat = rep(0.5, nrow(group.weight.dataset)),
#                          nburn = 40000,
#                          nsim = 2000,
#                          nthin = 5,
#                          n_chains = 2,
#                          # n_threads was max((RcppParallel::defaultNumThreads() - 2)/n_cores, 1) (this uses all of the server)
#                          n_threads = 4,
#                          update_interval = 500,
#                          ntree_control = 200,
#                          sd_control = 2 * sd(group.weight.dataset$w.change),
#                          base_control = 0.95,
#                          power_control = 2,
#                          ntree_moderate = 200,
#                          sd_moderate = 2 * sd(group.weight.dataset$w.change),
#                          base_moderate = 0.95,
#                          power_moderate = 2,
#                          use_muscale = FALSE,
#                          use_tauscale = FALSE,
#                          include_pi = "none",
#                          save_tree_directory = paste0(output_path, "/additional_outcomes/trees_weight"))
#   
#   saveRDS(weight_bcf, paste0(output_path, "/additional_outcomes/weight_bcf.rds"))
#   
#   
# }
# 
# plot_sigma_weight_bcf <- ggplot() +
#   geom_path(aes(x = 1:length(weight_bcf$sigma), y = weight_bcf$sigma)) +
#   ggtitle("Sigma") +
#   ylab("Sigma") +
#   theme(axis.title.x = element_blank())
# 
# 
# w.changes_bcf <- colMeans(weight_bcf$mu)
# 
# # maximum number of deciles being tested
# quantiles <- length(levels(group.weight.dataset[,"intervals"]))
# predictions_weight_bcf <- vector()
# mnumber = c(1:quantiles)
# 
# for (i in mnumber) {
# 
#   group_w.changes_male_sglt2 <- w.changes_bcf[which(group.weight.dataset$intervals == levels(group.weight.dataset$intervals)[i] &
#                                                 group.weight.dataset$sex == "Male" &
#                                                 group.weight.dataset$drugclass == "SGLT2")]
# 
#   predictions_weight_bcf <- rbind(predictions_weight_bcf, cbind(mean = mean(group_w.changes_male_sglt2, na.rm = TRUE),
#                                           lci = quantile(group_w.changes_male_sglt2, probs = c(0.05)),
#                                           uci = quantile(group_w.changes_male_sglt2, probs = c(0.95)),
#                                           sex = "Male",
#                                           drugclass = "SGLT2i",
#                                           intervals = levels(group.weight.dataset$intervals)[i]))
# 
#   group_w.changes_male_glp1 <- w.changes_bcf[which(group.weight.dataset$intervals == levels(group.weight.dataset$intervals)[i] &
#                                                       group.weight.dataset$sex == "Male" &
#                                                       group.weight.dataset$drugclass == "GLP1")]
# 
#   predictions_weight_bcf <- rbind(predictions_weight_bcf, cbind(mean = mean(group_w.changes_male_glp1, na.rm = TRUE),
#                                           lci = quantile(group_w.changes_male_glp1, probs = c(0.05)),
#                                           uci = quantile(group_w.changes_male_glp1, probs = c(0.95)),
#                                           sex = "Male",
#                                           drugclass = "GLP1-RA",
#                                           intervals = levels(group.weight.dataset$intervals)[i]))
# 
#   group_w.changes_female_sglt2 <- w.changes_bcf[which(group.weight.dataset$intervals == levels(group.weight.dataset$intervals)[i] &
#                                                       group.weight.dataset$sex == "Female" &
#                                                       group.weight.dataset$drugclass == "SGLT2")]
# 
#   predictions_weight_bcf <- rbind(predictions_weight_bcf, cbind(mean = mean(group_w.changes_female_sglt2, na.rm = TRUE),
#                                           lci = quantile(group_w.changes_female_sglt2, probs = c(0.05)),
#                                           uci = quantile(group_w.changes_female_sglt2, probs = c(0.95)),
#                                           sex = "Female",
#                                           drugclass = "SGLT2i",
#                                           intervals = levels(group.weight.dataset$intervals)[i]))
# 
# 
#   group_w.changes_female_glp1 <- w.changes_bcf[which(group.weight.dataset$intervals == levels(group.weight.dataset$intervals)[i] &
#                                                      group.weight.dataset$sex == "Female" &
#                                                      group.weight.dataset$drugclass == "GLP1")]
# 
#   predictions_weight_bcf <- rbind(predictions_weight_bcf, cbind(mean = mean(group_w.changes_female_glp1, na.rm = TRUE),
#                                           lci = quantile(group_w.changes_female_glp1, probs = c(0.05)),
#                                           uci = quantile(group_w.changes_female_glp1, probs = c(0.95)),
#                                           sex = "Female",
#                                           drugclass = "GLP1-RA",
#                                           intervals = levels(group.weight.dataset$intervals)[i]))
# 
# }
# 
# 
# predictions_weight_bcf <- predictions_weight_bcf %>%
#   as.data.frame()
# 
# plot_weight_male <- rbind(
#   # cbind(mean = NA, lci = NA, uci = NA, sex = "Male", drugclass = "SGLT2i", intervals = "Predicted HbA1c benefit on SGLT2i"),
#   predictions_weight_bcf %>%
#     filter(sex == "Male") %>%
#     slice(c(1:6)),
#   # cbind(mean = NA, lci = NA, uci = NA, sex = "Male", drugclass = NA, intervals = "Predicted HbA1c benefit on GLP1-RA"),
#   predictions_weight_bcf %>%
#     filter(sex == "Male") %>%
#     slice(-c(1:6))
# ) %>%
#   as.data.frame() %>%
#   mutate(mean = as.numeric(mean),
#          lci = as.numeric(lci),
#          intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=2155)", ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=3252)", ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=9124)", ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=7276)", ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1484)", ifelse(intervals == "(5,31]", ">5 mmol/mol (n=731)", intervals)))))),
#          uci = as.numeric(uci)) %>%
#   rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
#   group_by(group) %>%
#   forestplot(ci.vertices = TRUE,
#              ci.vertices.height = 0.1,
#              boxsize = .2,
#              xlab = "Weight change (kg)") %>%
#   fp_add_header("Male (n=24022)") %>%
#   fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
#                default = gpar(vertices = TRUE))
# 
# 
# plot_weight_female <- rbind(
#   # cbind(mean = NA, lci = NA, uci = NA, sex = "Male", drugclass = "SGLT2i", intervals = "Predicted HbA1c benefit on SGLT2i"),
#   predictions_weight_bcf %>%
#     filter(sex == "Female") %>%
#     slice(c(1:6)),
#   # cbind(mean = NA, lci = NA, uci = NA, sex = "Male", drugclass = NA, intervals = "Predicted HbA1c benefit on GLP1-RA"),
#   predictions_weight_bcf %>%
#     filter(sex == "Female") %>%
#     slice(-c(1:6))
# ) %>%
#   as.data.frame() %>%
#   mutate(mean = as.numeric(mean),
#          lci = as.numeric(lci),
#          intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=815)", ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=704)", ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3155)", ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5285)", ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3554)", ifelse(intervals == "(5,31]", ">5 mmol/mol (n=2019)", intervals)))))),
#          uci = as.numeric(uci)) %>%
#   rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
#   group_by(group) %>%
#   forestplot(ci.vertices = TRUE,
#              ci.vertices.height = 0.1,
#              boxsize = .2,
#              xlab = "Weight change (kg)") %>%
#   fp_add_header("Female (n=15532)") %>%
#   fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
#                default = gpar(vertices = TRUE))
# 
# 
# pdf(width = 14, height = 4, "Plots/bcf_weight.pdf")
# # grid.newpage()
# # pushViewport(viewport(layout = grid.layout(nrow = 1,
# #                                            ncol = 2)))
# # # first plot
# # pushViewport(viewport(layout.pos.row = 1,
# #                       layout.pos.col = 1))
# # plot_weight_male
# # upViewport()
# # # second plot
# # pushViewport(viewport(layout.pos.row = 1,
# #                       layout.pos.col = 2))
# # plot_weight_female
# dev.off()

## Adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.weight.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_weight_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

formula <- paste0("w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preweight, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = " + "))

## Fit a bcf for the weight change
if (class(try(
  
  predictions_weight_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_adjusted <- stan_glm(formula, 
                     data = group.weight.dataset,
                     family = gaussian(link = "identity"),
                     prior = normal(0, 2),
                     prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_adjusted, paste0(output_path, "/additional_outcomes/models_weight_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_weight_adjusted, newdata = group.weight.dataset %>%
                            filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                    lci = quantile(male_sglt2, probs = c(0.05)),
                                                                    uci = quantile(male_sglt2, probs = c(0.95)),
                                                                    sex = "Male",
                                                                    drugclass = "SGLT2",
                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_weight_adjusted, newdata = group.weight.dataset %>%
                           filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                    lci = quantile(male_glp1, probs = c(0.05)),
                                                                    uci = quantile(male_glp1, probs = c(0.95)),
                                                                    sex = "Male",
                                                                    drugclass = "GLP1",
                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_weight_adjusted, newdata = group.weight.dataset %>%
                              filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                    lci = quantile(female_sglt2, probs = c(0.05)),
                                                                    uci = quantile(female_sglt2, probs = c(0.95)),
                                                                    sex = "Female",
                                                                    drugclass = "SGLT2",
                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_weight_adjusted, newdata = group.weight.dataset %>%
                             filter(intervals == levels(group.weight.dataset$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_weight_stan_adjusted <- rbind(predictions_weight_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                    lci = quantile(female_glp1, probs = c(0.05)),
                                                                    uci = quantile(female_glp1, probs = c(0.95)),
                                                                    sex = "Female",
                                                                    drugclass = "GLP1",
                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_adjusted <- predictions_weight_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_weight_stan_adjusted.rds"))
  
}


## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.weight.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_weight_stan_psm_1_1 <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

formula <- "w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preweight, 3) + sex"


matching_weight <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ preweight + agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_weight_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_psm_1_1 <- stan_glm(formula, 
                                     data = group.weight.dataset.matched,
                                     family = gaussian(link = "identity"),
                                     prior = normal(0, 2),
                                     prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_psm_1_1, paste0(output_path, "/additional_outcomes/models_weight_psm_1_1.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                            filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                      uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                      sex = "Male",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                           filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                      lci = quantile(male_glp1, probs = c(0.05)),
                                                                                      uci = quantile(male_glp1, probs = c(0.95)),
                                                                                      sex = "Male",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                              filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                      uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                      sex = "Female",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_weight_psm_1_1, newdata = group.weight.dataset.matched %>%
                             filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_weight_stan_psm_1_1 <- rbind(predictions_weight_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                      lci = quantile(female_glp1, probs = c(0.05)),
                                                                                      uci = quantile(female_glp1, probs = c(0.95)),
                                                                                      sex = "Female",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1 <- predictions_weight_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1.rds"))
  
}


## Propensity score matching + adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.weight.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_weight_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.weight.dataset[,breakdown_adjust], is.factor)

formula <- paste0("w.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preweight, 3) + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) +", paste(breakdown_adjust[factors], collapse = "+"))

matching_weight <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ preweight + agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_weight_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_weight_psm_1_1_adjusted <- stan_glm(formula, 
                                    data = group.weight.dataset.matched,
                                    family = gaussian(link = "identity"),
                                    prior = normal(0, 2),
                                    prior_intercept = normal(0, 2))
  
  saveRDS(models_weight_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_weight_psm_1_1_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                            filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                    uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                    sex = "Male",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                           filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                    lci = quantile(male_glp1, probs = c(0.05)),
                                                                                    uci = quantile(male_glp1, probs = c(0.95)),
                                                                                    sex = "Male",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                              filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                    uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                    sex = "Female",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_weight_psm_1_1_adjusted, newdata = group.weight.dataset.matched %>%
                             filter(intervals == levels(group.weight.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_weight_stan_psm_1_1_adjusted <- rbind(predictions_weight_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                    lci = quantile(female_glp1, probs = c(0.05)),
                                                                                    uci = quantile(female_glp1, probs = c(0.95)),
                                                                                    sex = "Female",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.weight.dataset$intervals)[i]))
    
  }
  
  predictions_weight_stan_psm_1_1_adjusted <- predictions_weight_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_weight_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_weight_stan_psm_1_1_adjusted.rds"))
  
}

#:--------------------
#:---- PLOTS
#:--------------------

## limits

weight_axis_min <- plyr::round_any(floor(min(c(predictions_weight_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                               predictions_weight_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                               predictions_weight_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

weight_axis_max <- plyr::round_any(ceiling(max(c(predictions_weight_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                 predictions_weight_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                 predictions_weight_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


#:---- Adjusted

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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=2155)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=3252)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=9124)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=7276)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1484)", 
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=731)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Male (n=24022)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_weight_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=815)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=704)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3155)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5285)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3554)", 
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=2019)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Female (n=15532)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1

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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=987)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1406)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3298)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2219)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=479)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=260)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Male (n=8649)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_weight_psm_1_1_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=439)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=356)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1510)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2173)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1348)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=823)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Female (n=6649)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1 adjusted

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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=987)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1406)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3298)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2219)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=479)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=260)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Male (n=8649)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_weight_psm_1_1_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_weight_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=439)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=356)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1510)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2173)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1348)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=823)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(weight_axis_min, weight_axis_max),
             xticks = seq(weight_axis_min, weight_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Weight change (kg)") %>%
  fp_add_header("Female (n=6649)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


pdf(width = 14, height = 12, "Plots/weight.pdf")
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
# fourth plot
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

## Read in data for weight
egfr.dataset <- set_up_data_sglt2_glp1(dataset.type = "egfr.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(egfr.change = postegfr - preegfr)


breaks_egfr <- c(-5, -3, 0, 3, 5)

group.egfr.dataset <- group_values(data = egfr.dataset,
                                     variable = "effects",
                                     breaks = breaks_egfr) %>%
  drop_na(intervals)


## Adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_egfr_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.egfr.dataset[,breakdown_adjust], is.factor)

formula <- paste0("egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) + ", paste(breakdown_adjust[factors], collapse = "+"))

## Fit a bcf for the egfr change
if (class(try(
  
  predictions_egfr_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_adjusted <- stan_glm(formula, 
                                     data = group.egfr.dataset,
                                     family = gaussian(link = "identity"),
                                     prior = normal(0, 2),
                                     prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_adjusted, paste0(output_path, "/additional_outcomes/models_egfr_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                            filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                      uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                      sex = "Male",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                           filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                      lci = quantile(male_glp1, probs = c(0.05)),
                                                                                      uci = quantile(male_glp1, probs = c(0.95)),
                                                                                      sex = "Male",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                              filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                      lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                      uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                      sex = "Female",
                                                                                      drugclass = "SGLT2",
                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_egfr_adjusted, newdata = group.egfr.dataset %>%
                             filter(intervals == levels(group.egfr.dataset$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_adjusted <- rbind(predictions_egfr_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                      lci = quantile(female_glp1, probs = c(0.05)),
                                                                                      uci = quantile(female_glp1, probs = c(0.95)),
                                                                                      sex = "Female",
                                                                                      drugclass = "GLP1",
                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_adjusted <- predictions_egfr_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_adjusted.rds"))
  
}


## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_egfr_stan_psm_1_1 <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.egfr.dataset[,breakdown_adjust], is.factor)

formula <- "egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(preegfr, 3) + sex"

matching_egfr <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_egfr_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_psm_1_1 <- stan_glm(formula, 
                                    data = group.egfr.dataset.matched,
                                    family = gaussian(link = "identity"),
                                    prior = normal(0, 2),
                                    prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_psm_1_1, paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                            filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                    uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                    sex = "Male",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                           filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                    lci = quantile(male_glp1, probs = c(0.05)),
                                                                                    uci = quantile(male_glp1, probs = c(0.95)),
                                                                                    sex = "Male",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                              filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                    lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                    uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                    sex = "Female",
                                                                                    drugclass = "SGLT2",
                                                                                    intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_egfr_psm_1_1, newdata = group.egfr.dataset.matched %>%
                             filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_psm_1_1 <- rbind(predictions_egfr_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                    lci = quantile(female_glp1, probs = c(0.05)),
                                                                                    uci = quantile(female_glp1, probs = c(0.95)),
                                                                                    sex = "Female",
                                                                                    drugclass = "GLP1",
                                                                                    intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1 <- predictions_egfr_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1.rds"))
  
}


## Propensity score matching + adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.egfr.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_egfr_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.egfr.dataset[,breakdown_adjust], is.factor)

formula <- paste0("egfr.change ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(breakdown_adjust[factors], collapse = "+"))

matching_egfr <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_egfr_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_egfr_psm_1_1_adjusted <- stan_glm(formula, 
                                             data = group.egfr.dataset.matched,
                                             family = gaussian(link = "identity"),
                                             prior = normal(0, 2),
                                             prior_intercept = normal(0, 2))
  
  saveRDS(models_egfr_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_egfr_psm_1_1_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                            filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                                      uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                           filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(male_glp1, probs = c(0.05)),
                                                                                                      uci = quantile(male_glp1, probs = c(0.95)),
                                                                                                      sex = "Male",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                              filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"))
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                      lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                                      uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "SGLT2",
                                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_egfr_psm_1_1_adjusted, newdata = group.egfr.dataset.matched %>%
                             filter(intervals == levels(group.egfr.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"))
    
    predictions_egfr_stan_psm_1_1_adjusted <- rbind(predictions_egfr_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                      lci = quantile(female_glp1, probs = c(0.05)),
                                                                                                      uci = quantile(female_glp1, probs = c(0.95)),
                                                                                                      sex = "Female",
                                                                                                      drugclass = "GLP1",
                                                                                                      intervals = levels(group.egfr.dataset$intervals)[i]))
    
  }
  
  predictions_egfr_stan_psm_1_1_adjusted <- predictions_egfr_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_egfr_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_egfr_stan_psm_1_1_adjusted.rds"))
  
}

#:--------------------
#:---- PLOTS
#:--------------------

## limits

egfr_axis_min <- plyr::round_any(floor(min(c(predictions_egfr_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(), 
                                               predictions_egfr_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
                                               predictions_egfr_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()))), 2, f = floor)

egfr_axis_max <- plyr::round_any(ceiling(max(c(predictions_egfr_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(), 
                                                 predictions_egfr_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
                                                 predictions_egfr_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()))), 2, f = ceiling)


#:---- Adjusted


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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=2399)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=3541)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=9997)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=7947)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1678)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=870)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Male (n=26432)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_egfr_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=872)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=794)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3389)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5781)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3964)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=2294)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Female (n=17094)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1


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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1080)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1506)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3513)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2320)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=511)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=300)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Male (n=9230)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_egfr_psm_1_1_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=433)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=380)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1568)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2298)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1408)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=947)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Female (n=7034)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1 adjusted


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
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1080)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1506)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3513)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2320)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=511)", 
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=300)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Male (n=9230)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_egfr_psm_1_1_adjusted_female <- rbind(
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = 50, lci = 50, uci = 50, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_egfr_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=433)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=380)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1568)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2298)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1408)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=947)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(egfr_axis_min, egfr_axis_max),
             xticks = seq(egfr_axis_min, egfr_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "eGFR change") %>%
  fp_add_header("Female (n=7034)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


pdf(width = 14, height = 12, "Plots/egfr.pdf")
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
# fourth plot
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

breaks_discontinuation <- c(-5, -3, 0, 3, 5)

group.discontinuation.dataset <- group_values(data = discontinuation.dataset,
                                   variable = "effects",
                                   breaks = breaks_discontinuation) %>%
  drop_na(intervals)


## Adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_discontinuation_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

formula <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) + ", paste(breakdown_adjust[factors], collapse = "+"))

## Fit a bcf for the discontinuation change
if (class(try(
  
  predictions_discontinuation_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_discontinuation_adjusted <- stan_glm(formula, 
                                   data = group.discontinuation.dataset,
                                   family = binomial(link = "logit"),
                                   prior = normal(0, 2),
                                   prior_intercept = normal(0, 2))
  
  saveRDS(models_discontinuation_adjusted, paste0(output_path, "/additional_outcomes/models_discontinuation_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                            filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"),
                          type = "response")
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                  lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                  uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                  sex = "Male",
                                                                                  drugclass = "SGLT2",
                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                           filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"),
                         type = "response")
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                  lci = quantile(male_glp1, probs = c(0.05)),
                                                                                  uci = quantile(male_glp1, probs = c(0.95)),
                                                                                  sex = "Male",
                                                                                  drugclass = "GLP1",
                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                              filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"),
                            type = "response")
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                  lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                  uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                  sex = "Female",
                                                                                  drugclass = "SGLT2",
                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_discontinuation_adjusted, newdata = group.discontinuation.dataset %>%
                             filter(intervals == levels(group.discontinuation.dataset$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"),
                           type = "response")
    
    predictions_discontinuation_stan_adjusted <- rbind(predictions_discontinuation_stan_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                  lci = quantile(female_glp1, probs = c(0.05)),
                                                                                  uci = quantile(female_glp1, probs = c(0.95)),
                                                                                  sex = "Female",
                                                                                  drugclass = "GLP1",
                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_adjusted <- predictions_discontinuation_stan_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_adjusted.rds"))
  
}



## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_discontinuation_stan_psm_1_1 <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

formula <- "stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + sex"

matching_discontinuation <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_discontinuation_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
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
                            filter(drugclass == "SGLT2") %>%
                            slice(1),
                          type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                sex = "Male",
                                                                                drugclass = "SGLT2",
                                                                                intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                           filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1") %>%
                             slice(1),
                         type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                lci = quantile(male_glp1, probs = c(0.05)),
                                                                                uci = quantile(male_glp1, probs = c(0.95)),
                                                                                sex = "Male",
                                                                                drugclass = "GLP1",
                                                                                intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                              filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2") %>%
                                slice(1),
                            type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                sex = "Female",
                                                                                drugclass = "SGLT2",
                                                                                intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- posterior_epred(models_discontinuation_psm_1_1, newdata = group.discontinuation.dataset.matched %>%
                             filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1") %>%
                               slice(1),
                           type = "response")
    
    predictions_discontinuation_stan_psm_1_1 <- rbind(predictions_discontinuation_stan_psm_1_1, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                lci = quantile(female_glp1, probs = c(0.05)),
                                                                                uci = quantile(female_glp1, probs = c(0.95)),
                                                                                sex = "Female",
                                                                                drugclass = "GLP1",
                                                                                intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1 <- predictions_discontinuation_stan_psm_1_1 %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1.rds"))
  
}

## Propensity score matching + adjustment

# maximum number of deciles being tested
quantiles <- length(levels(group.discontinuation.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_discontinuation_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

formula <- paste0("stopdrug_6m_3mFU ~ factor(drugclass) + intervals + factor(drugclass)*intervals + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(breakdown_adjust[factors], collapse = "+"))

matching_discontinuation <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~ agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt +", paste(breakdown_adjust[factors], collapse = "+"))),
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

if (class(try(
  
  predictions_discontinuation_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  models_discontinuation_psm_1_1_adjusted <- stan_glm(formula, 
                                           data = group.discontinuation.dataset.matched,
                                           family = binomial(link = "logit"),
                                           prior = normal(0, 2),
                                           prior_intercept = normal(0, 2))
  
  saveRDS(models_discontinuation_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/models_discontinuation_psm_1_1_adjusted.rds"))
  
  for (i in mnumber) {
    
    male_sglt2 <- predict(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                            filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                            filter(sex == "Male") %>%
                            filter(drugclass == "SGLT2"),
                          type = "response")
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(male_sglt2, na.rm = TRUE),
                                                                                                  lci = quantile(male_sglt2, probs = c(0.05)),
                                                                                                  uci = quantile(male_sglt2, probs = c(0.95)),
                                                                                                  sex = "Male",
                                                                                                  drugclass = "SGLT2",
                                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    male_glp1 <- predict(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                           filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                           filter(sex == "Male") %>%
                           filter(drugclass == "GLP1"),
                         type = "response")
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(male_glp1, na.rm = TRUE),
                                                                                                  lci = quantile(male_glp1, probs = c(0.05)),
                                                                                                  uci = quantile(male_glp1, probs = c(0.95)),
                                                                                                  sex = "Male",
                                                                                                  drugclass = "GLP1",
                                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_sglt2 <- predict(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                              filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                              filter(sex == "Female") %>%
                              filter(drugclass == "SGLT2"),
                            type = "response")
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(female_sglt2, na.rm = TRUE),
                                                                                                  lci = quantile(female_sglt2, probs = c(0.05)),
                                                                                                  uci = quantile(female_sglt2, probs = c(0.95)),
                                                                                                  sex = "Female",
                                                                                                  drugclass = "SGLT2",
                                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
    female_glp1 <- predict(models_discontinuation_psm_1_1_adjusted, newdata = group.discontinuation.dataset.matched %>%
                             filter(intervals == levels(group.discontinuation.dataset.matched$intervals)[i]) %>%
                             filter(sex == "Female") %>%
                             filter(drugclass == "GLP1"),
                           type = "response")
    
    predictions_discontinuation_stan_psm_1_1_adjusted <- rbind(predictions_discontinuation_stan_psm_1_1_adjusted, cbind(mean = mean(female_glp1, na.rm = TRUE),
                                                                                                  lci = quantile(female_glp1, probs = c(0.05)),
                                                                                                  uci = quantile(female_glp1, probs = c(0.95)),
                                                                                                  sex = "Female",
                                                                                                  drugclass = "GLP1",
                                                                                                  intervals = levels(group.discontinuation.dataset$intervals)[i]))
    
  }
  
  predictions_discontinuation_stan_psm_1_1_adjusted <- predictions_discontinuation_stan_psm_1_1_adjusted %>%
    as.data.frame()
  
  saveRDS(predictions_discontinuation_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_discontinuation_stan_psm_1_1_adjusted.rds"))
  
}

#:--------------------
#:---- PLOTS
#:--------------------

## limits

discontinuation_axis_max <- plyr::round_any(ceiling(max(c(predictions_discontinuation_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100, 
                                               predictions_discontinuation_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100,
                                               predictions_discontinuation_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max() * 100))), 2, f = ceiling)


#:---- Adjusted

plot_discontinuation_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=2642)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=3896)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=10916)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=8766)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1844)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=934)", intervals)))))),
         uci = as.numeric(uci) * 100,) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Male (n=28998)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_discontinuation_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=974)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=884)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3744)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=6311)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=4347)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=2502)", intervals)))))),
         uci = as.numeric(uci) * 100) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Adjusted",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Female (n=18762)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1

plot_discontinuation_psm_1_1_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1203)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1665)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3915)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2571)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=560)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=326)", intervals)))))),
         uci = as.numeric(uci) * 100) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Male (n=10240)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_discontinuation_psm_1_1_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=511)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=428)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1741)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2538)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1565)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=1043)", intervals)))))),
         uci = as.numeric(uci) * 100) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Female (n=7826)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


#:---- PSM 1:1 + adjusted


plot_discontinuation_psm_1_1_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1203)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1665)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3915)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2571)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=560)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=326)", intervals)))))),
         uci = as.numeric(uci) * 100) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Male (n=7826)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))


plot_discontinuation_psm_1_1_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on SGLT2i"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:6)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "SGLT2", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", drugclass = "GLP1", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_discontinuation_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:6))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean) * 100,
         lci = as.numeric(lci) * 100,
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=511)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=428)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1741)", 
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2538)", 
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1565)", 
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=1043)", intervals)))))),
         uci = as.numeric(uci) * 100) %>%
  rename("lower" = "lci", "upper" = "uci", "group" = "drugclass", "labeltext" = "intervals") %>%
  group_by(group) %>%
  forestplot(ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             title = "Propensity score matching + adjusted",
             clip = c(0, discontinuation_axis_max),
             xticks = seq(0, discontinuation_axis_max, 2),
             boxsize = .2,
             txt_gp = fpTxtGp(ticks=gpar(cex=0.8), xlab=gpar(cex=1)),
             xlab = "Discontinuation") %>%
  fp_add_header("Female (n=7826)") %>%
  fp_set_style(box = c("#f1a340", "dodgerblue2") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE))



pdf(width = 14, height = 12, "Plots/discontinuation.pdf")
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
# fourth plot
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


breaks_no_co <- c(-5, -3, 0, 3, 5)

group.no_co.dataset <- group_values(data = no_co.dataset,
                                    variable = "effects",
                                    breaks = breaks_no_co) %>%
  drop_na(intervals)

# CVD survival analysis

## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_cvd_stan_psm_1_1 <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.discontinuation.dataset[,breakdown_adjust], is.factor)

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_no_co_cvd_psm_1_1_male <- vector("list", quantiles)

models_no_co_cvd_psm_1_1_female <- vector("list", quantiles)


if (class(try(

  predictions_no_co_cvd_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1.rds"))

  , silent = TRUE)) == "try-error") {

  for (i in mnumber) {

    models_no_co_cvd_psm_1_1_male[[i]] <- coxph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass,
                                            data = group.no_co.dataset.matched %>%
                                              slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
                                              filter(sex == "Male"))

    
    predictions_no_co_cvd_stan_psm_1_1 <- rbind(predictions_no_co_cvd_stan_psm_1_1, cbind(mean = models_no_co_cvd_psm_1_1_male[[i]]$coefficients[1],
                                                                                lci = confint(models_no_co_cvd_psm_1_1_male[[i]])[1],
                                                                                uci = confint(models_no_co_cvd_psm_1_1_male[[i]])[2],
                                                                                sex = "Male",
                                                                                intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    
    models_no_co_cvd_psm_1_1_female[[i]] <- coxph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass,
                                            data = group.no_co.dataset.matched %>%
                                              slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
                                              filter(sex == "Female"))
    
    
    predictions_no_co_cvd_stan_psm_1_1 <- rbind(predictions_no_co_cvd_stan_psm_1_1, cbind(mean = models_no_co_cvd_psm_1_1_female[[i]]$coefficients[1],
                                                                                  lci = confint(models_no_co_cvd_psm_1_1_female[[i]])[1],
                                                                                  uci = confint(models_no_co_cvd_psm_1_1_female[[i]])[2],
                                                                                  sex = "Female",
                                                                                  intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    


  }
  
  
  predictions_no_co_cvd_stan_psm_1_1 <- predictions_no_co_cvd_stan_psm_1_1 %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))

  
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1.rds"))

}


## Propensity score matching + adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_cvd_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE


matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_no_co_cvd_psm_1_1_adjusted_male <- vector("list", quantiles)

models_no_co_cvd_psm_1_1_adjusted_female <- vector("list", quantiles)

if (class(try(
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.no_co.dataset.matched %>%
      slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    
    models_no_co_cvd_psm_1_1_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                data = data.new)
    
    
    predictions_no_co_cvd_stan_psm_1_1_adjusted <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_cvd_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                          uci = confint(models_no_co_cvd_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                          sex = "Male",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    # female data
    data.new <- group.no_co.dataset.matched %>%
      slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    models_no_co_cvd_psm_1_1_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                  data = data.new)
    
    
    predictions_no_co_cvd_stan_psm_1_1_adjusted <- rbind(predictions_no_co_cvd_stan_psm_1_1_adjusted, cbind(mean = models_no_co_cvd_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_cvd_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                          uci = confint(models_no_co_cvd_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                          sex = "Female",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_no_co_cvd_stan_psm_1_1_adjusted <- predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_no_co_cvd_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_psm_1_1_adjusted.rds"))
  
}


## Adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_cvd_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

models_no_co_cvd_adjusted_male <- vector("list", quantiles)

models_no_co_cvd_adjusted_female <- vector("list", quantiles)



if (class(try(

  predictions_no_co_cvd_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted.rds"))

  , silent = TRUE)) == "try-error") {

  for (i in mnumber) {
    
    # male data
    data.new <- group.no_co.dataset %>%
      slice(which(group.no_co.dataset$intervals == levels(group.no_co.dataset$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_no_co_cvd_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                         data = data.new)


    predictions_no_co_cvd_stan_adjusted <- rbind(predictions_no_co_cvd_stan_adjusted, cbind(mean = models_no_co_cvd_adjusted_male[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_cvd_adjusted_male[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_cvd_adjusted_male[[i]])[1,2],
                                                                                                            sex = "Male",
                                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))

    # female data
    data.new <- group.no_co.dataset %>%
      slice(which(group.no_co.dataset$intervals == levels(group.no_co.dataset$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_no_co_cvd_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                           data = data.new)


    predictions_no_co_cvd_stan_adjusted <- rbind(predictions_no_co_cvd_stan_adjusted, cbind(mean = models_no_co_cvd_adjusted_female[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_cvd_adjusted_female[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_cvd_adjusted_female[[i]])[1,2],
                                                                                                            sex = "Female",
                                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))




  }


  predictions_no_co_cvd_stan_adjusted <- predictions_no_co_cvd_stan_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))



  saveRDS(predictions_no_co_cvd_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_cvd_stan_adjusted.rds"))

}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

no_co_axis_min <- plyr::round_any(ceiling(min(c(
  predictions_no_co_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  # predictions_no_co_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  predictions_no_co_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  1))
  ), 2, f = floor)


no_co_axis_max <- plyr::round_any(ceiling(max(c(
  predictions_no_co_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  # predictions_no_co_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  predictions_no_co_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  1))
  ), 2, f = ceiling)

#:---- PSM 1:1

plot_no_co_cvd_psm_1_1_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
  ) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=738)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=931)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2194)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1336)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=168)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=80)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5447)")



plot_no_co_cvd_psm_1_1_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=355)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=274)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1145)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1614)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1074)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=447)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=4909)")


#:---- PSM 1:1 + adjusted

plot_no_co_cvd_psm_1_1_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=738)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=931)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2194)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1336)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=168)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=80)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5447)")



plot_no_co_cvd_psm_1_1_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=355)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=274)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1145)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1614)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1074)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=447)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=4909)")


#:---- Adjusted

plot_no_co_cvd_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1852)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=2709)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=7648)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5774)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=770)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=330)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=19083)")



plot_no_co_cvd_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_cvd_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_cvd_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=740)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=661)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2846)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=4631)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3369)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=1359)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Adjusted",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=13606)")



pdf(width = 14, height = 12, "Plots/no_comorbidities_cvd.pdf")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("CVD outcomes for no comorbidities population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
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
# fourth plot
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

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_hf_stan_psm_1_1 <- vector()

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_no_co_hf_psm_1_1_male <- vector("list", quantiles)

models_no_co_hf_psm_1_1_female <- vector("list", quantiles)


if (class(try(
  
  predictions_no_co_hf_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    models_no_co_hf_psm_1_1_male[[i]] <- coxph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass,
                                                data = group.no_co.dataset.matched %>%
                                                  slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
                                                  filter(sex == "Male"))
    
    
    predictions_no_co_hf_stan_psm_1_1 <- rbind(predictions_no_co_hf_stan_psm_1_1, cbind(mean = models_no_co_hf_psm_1_1_male[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_hf_psm_1_1_male[[i]])[1],
                                                                                          uci = confint(models_no_co_hf_psm_1_1_male[[i]])[2],
                                                                                          sex = "Male",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    
    models_no_co_hf_psm_1_1_female[[i]] <- coxph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass,
                                                  data = group.no_co.dataset.matched %>%
                                                    slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
                                                    filter(sex == "Female"))
    
    
    predictions_no_co_hf_stan_psm_1_1 <- rbind(predictions_no_co_hf_stan_psm_1_1, cbind(mean = models_no_co_hf_psm_1_1_female[[i]]$coefficients[1],
                                                                                          lci = confint(models_no_co_hf_psm_1_1_female[[i]])[1],
                                                                                          uci = confint(models_no_co_hf_psm_1_1_female[[i]])[2],
                                                                                          sex = "Female",
                                                                                          intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_no_co_hf_stan_psm_1_1 <- predictions_no_co_hf_stan_psm_1_1 %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1.rds"))
  
}


## Propensity score matching + adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_hf_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

matching_no_co <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_no_co_hf_psm_1_1_adjusted_male <- vector("list", quantiles)

models_no_co_hf_psm_1_1_adjusted_female <- vector("list", quantiles)

if (class(try(
  
  predictions_no_co_hf_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.no_co.dataset.matched %>%
      slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    
    models_no_co_hf_psm_1_1_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                         data = data.new)
    
    
    predictions_no_co_hf_stan_psm_1_1_adjusted <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted, cbind(mean = models_no_co_hf_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_hf_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_hf_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                            sex = "Male",
                                                                                                            intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    # female data
    data.new <- group.no_co.dataset.matched %>%
      slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    models_no_co_hf_psm_1_1_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                           data = group.no_co.dataset.matched %>%
                                                             slice(which(group.no_co.dataset.matched$intervals == levels(group.no_co.dataset.matched$intervals)[i])) %>%
                                                             filter(sex == "Female"))
    
    
    predictions_no_co_hf_stan_psm_1_1_adjusted <- rbind(predictions_no_co_hf_stan_psm_1_1_adjusted, cbind(mean = models_no_co_hf_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                            lci = confint(models_no_co_hf_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                            uci = confint(models_no_co_hf_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                            sex = "Female",
                                                                                                            intervals = levels(group.no_co.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_no_co_hf_stan_psm_1_1_adjusted <- predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_no_co_hf_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_psm_1_1_adjusted.rds"))
  
}


## Adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.no_co.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_no_co_hf_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.no_co.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

models_no_co_hf_adjusted_male <- vector("list", quantiles)

models_no_co_hf_adjusted_female <- vector("list", quantiles)



if (class(try(
  
  predictions_no_co_hf_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.no_co.dataset %>%
      slice(which(group.no_co.dataset$intervals == levels(group.no_co.dataset$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_no_co_hf_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                 data = data.new)
    
    
    predictions_no_co_hf_stan_adjusted <- rbind(predictions_no_co_hf_stan_adjusted, cbind(mean = models_no_co_hf_adjusted_male[[i]]$coefficients[1],
                                                                                            lci = confint(models_no_co_hf_adjusted_male[[i]])[1,1],
                                                                                            uci = confint(models_no_co_hf_adjusted_male[[i]])[1,2],
                                                                                            sex = "Male",
                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))
    
    # female data
    data.new <- group.no_co.dataset %>%
      slice(which(group.no_co.dataset$intervals == levels(group.no_co.dataset$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_no_co_hf_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                   data = data.new)
    
    
    predictions_no_co_hf_stan_adjusted <- rbind(predictions_no_co_hf_stan_adjusted, cbind(mean = models_no_co_hf_adjusted_female[[i]]$coefficients[1],
                                                                                            lci = confint(models_no_co_hf_adjusted_female[[i]])[1,1],
                                                                                            uci = confint(models_no_co_hf_adjusted_female[[i]])[1,2],
                                                                                            sex = "Female",
                                                                                            intervals = levels(group.no_co.dataset$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_no_co_hf_stan_adjusted <- predictions_no_co_hf_stan_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_no_co_hf_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_no_co_hf_stan_adjusted.rds"))
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

no_co_axis_min <- plyr::round_any(ceiling(min(c(
  predictions_no_co_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  # predictions_no_co_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  predictions_no_co_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
  1))
), 2, f = floor)


no_co_axis_max <- plyr::round_any(ceiling(max(c(
  predictions_no_co_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  # predictions_no_co_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  predictions_no_co_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
  1))
), 2, f = ceiling)

#:---- PSM 1:1

plot_no_co_hf_psm_1_1_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=738)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=931)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2194)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1336)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=168)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=80)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5447)")



plot_no_co_hf_psm_1_1_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=355)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=274)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1145)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1614)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1074)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=447)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=4909)")


#:---- PSM 1:1 + adjusted

plot_no_co_hf_psm_1_1_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=738)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=931)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2194)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1336)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=168)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=80)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5447)")



plot_no_co_hf_psm_1_1_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=355)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=274)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1145)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1614)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1074)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=447)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=4909)")


#:---- Adjusted

plot_no_co_hf_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1852)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=2709)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=7648)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5774)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=770)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=330)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=19083)")



plot_no_co_hf_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_no_co_hf_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_no_co_hf_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=740)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=661)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2846)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=4631)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3369)",
                                                        ifelse(intervals == "(5,28]", ">5 mmol/mol (n=1359)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Adjusted",
             clip = c(no_co_axis_min, no_co_axis_max),
             xticks = seq(no_co_axis_min, no_co_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=13606)")



pdf(width = 14, height = 12, "Plots/no_comorbidities_hf.pdf")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4,
                                           ncol = 2, heights = unit(c(0.5, 5, 5, 5), "null"))))
# title
grid.text("HF outcomes for no comorbidities population (GLP1 control, values correspond to SGLT2)", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
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
# fourth plot
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


breaks_cvd <- c(-5, -3, 0, 3, 5)

group.cvd.dataset <- group_values(data = cvd.dataset,
                                    variable = "effects",
                                    breaks = breaks_cvd) %>%
  drop_na(intervals)





# CVD survival analysis

## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_cvd_cvd_stan_psm_1_1 <- vector()

matching_cvd <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_cvd_cvd_psm_1_1_male <- vector("list", quantiles)

models_cvd_cvd_psm_1_1_female <- vector("list", quantiles)


if (class(try(
  
  predictions_cvd_cvd_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    models_cvd_cvd_psm_1_1_male[[i]] <- coxph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass,
                                                data = group.cvd.dataset.matched %>%
                                                  slice(which(group.cvd.dataset.matched$intervals == levels(group.cvd.dataset.matched$intervals)[i])) %>%
                                                  filter(sex == "Male"))
    
    
    predictions_cvd_cvd_stan_psm_1_1 <- rbind(predictions_cvd_cvd_stan_psm_1_1, cbind(mean = models_cvd_cvd_psm_1_1_male[[i]]$coefficients[1],
                                                                                          lci = confint(models_cvd_cvd_psm_1_1_male[[i]])[1],
                                                                                          uci = confint(models_cvd_cvd_psm_1_1_male[[i]])[2],
                                                                                          sex = "Male",
                                                                                          intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    
    models_cvd_cvd_psm_1_1_female[[i]] <- coxph(Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass,
                                                  data = group.cvd.dataset.matched %>%
                                                    slice(which(group.cvd.dataset.matched$intervals == levels(group.cvd.dataset.matched$intervals)[i])) %>%
                                                    filter(sex == "Female"))
    
    
    predictions_cvd_cvd_stan_psm_1_1 <- rbind(predictions_cvd_cvd_stan_psm_1_1, cbind(mean = models_cvd_cvd_psm_1_1_female[[i]]$coefficients[1],
                                                                                          lci = confint(models_cvd_cvd_psm_1_1_female[[i]])[1],
                                                                                          uci = confint(models_cvd_cvd_psm_1_1_female[[i]])[2],
                                                                                          sex = "Female",
                                                                                          intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_cvd_cvd_stan_psm_1_1 <- predictions_cvd_cvd_stan_psm_1_1 %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1.rds"))
  
}


## Propensity score matching + adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_cvd_cvd_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.cvd.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

matching_cvd <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_cvd_cvd_psm_1_1_adjusted_male <- vector("list", quantiles)

models_cvd_cvd_psm_1_1_adjusted_female <- vector("list", quantiles)

if (class(try(
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.cvd.dataset.matched %>%
      slice(which(group.cvd.dataset.matched$intervals == levels(group.cvd.dataset.matched$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    
    models_cvd_cvd_psm_1_1_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                         data = data.new)
    
    
    predictions_cvd_cvd_stan_psm_1_1_adjusted <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                            lci = confint(models_cvd_cvd_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                            uci = confint(models_cvd_cvd_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                            sex = "Male",
                                                                                                            intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    # female data
    data.new <- group.cvd.dataset.matched %>%
      slice(which(group.cvd.dataset.matched$intervals == levels(group.cvd.dataset.matched$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    models_cvd_cvd_psm_1_1_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                           data = group.cvd.dataset.matched %>%
                                                             slice(which(group.cvd.dataset.matched$intervals == levels(group.cvd.dataset.matched$intervals)[i])) %>%
                                                             filter(sex == "Female"))
    
    
    predictions_cvd_cvd_stan_psm_1_1_adjusted <- rbind(predictions_cvd_cvd_stan_psm_1_1_adjusted, cbind(mean = models_cvd_cvd_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                            lci = confint(models_cvd_cvd_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                            uci = confint(models_cvd_cvd_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                            sex = "Female",
                                                                                                            intervals = levels(group.cvd.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_cvd_cvd_stan_psm_1_1_adjusted <- predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_cvd_cvd_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_psm_1_1_adjusted.rds"))
  
}


## Adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.cvd.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_cvd_cvd_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.cvd.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

models_cvd_cvd_adjusted_male <- vector("list", quantiles)

models_cvd_cvd_adjusted_female <- vector("list", quantiles)



if (class(try(
  
  predictions_cvd_cvd_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.cvd.dataset %>%
      slice(which(group.cvd.dataset$intervals == levels(group.cvd.dataset$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_cvd_cvd_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                 data = data.new)
    
    
    predictions_cvd_cvd_stan_adjusted <- rbind(predictions_cvd_cvd_stan_adjusted, cbind(mean = models_cvd_cvd_adjusted_male[[i]]$coefficients[1],
                                                                                            lci = confint(models_cvd_cvd_adjusted_male[[i]])[1,1],
                                                                                            uci = confint(models_cvd_cvd_adjusted_male[[i]])[1,2],
                                                                                            sex = "Male",
                                                                                            intervals = levels(group.cvd.dataset$intervals)[i]))
    
    # female data
    data.new <- group.cvd.dataset %>%
      slice(which(group.cvd.dataset$intervals == levels(group.cvd.dataset$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_mace_censtime_yrs, postdrug_mace_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_cvd_cvd_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                   data = data.new)
    
    
    predictions_cvd_cvd_stan_adjusted <- rbind(predictions_cvd_cvd_stan_adjusted, cbind(mean = models_cvd_cvd_adjusted_female[[i]]$coefficients[1],
                                                                                            lci = confint(models_cvd_cvd_adjusted_female[[i]])[1,1],
                                                                                            uci = confint(models_cvd_cvd_adjusted_female[[i]])[1,2],
                                                                                            sex = "Female",
                                                                                            intervals = levels(group.cvd.dataset$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_cvd_cvd_stan_adjusted <- predictions_cvd_cvd_stan_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_cvd_cvd_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_cvd_cvd_stan_adjusted.rds"))
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

# cvd_axis_min <- plyr::round_any(ceiling(min(c(
#   predictions_cvd_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
#   predictions_cvd_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
#   predictions_cvd_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()),
#   1)
# ), 2, f = floor)

cvd_axis_min = -4


# cvd_axis_max <- plyr::round_any(ceiling(max(c(
#   predictions_cvd_cvd_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
#   predictions_cvd_cvd_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
#   predictions_cvd_cvd_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()),
#   1)
# ), 2, f = ceiling)

cvd_axis_max = 4

#:---- PSM 1:1

plot_cvd_cvd_psm_1_1_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=733)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=974)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2240)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1442)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=249)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=157)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5795)")



plot_cvd_cvd_psm_1_1_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=348)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=280)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1167)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1686)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1161)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=705)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=5347)")


#:---- PSM 1:1 + adjusted

plot_cvd_cvd_psm_1_1_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=733)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=974)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2240)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1442)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=249)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=157)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=5795)")



plot_cvd_cvd_psm_1_1_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=348)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=280)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1167)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1686)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1161)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=705)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=5347)")


#:---- Adjusted

plot_cvd_cvd_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=1891)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=2741)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=7818)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=6019)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=880)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=478)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=19827)")



plot_cvd_cvd_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_cvd_cvd_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_cvd_cvd_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=750)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=674)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2876)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=4823)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=3541)",
                                                        ifelse(intervals == "(5,31]", ">5 mmol/mol (n=1797)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Adjusted",
             clip = c(cvd_axis_min, cvd_axis_max),
             xticks = seq(cvd_axis_min, cvd_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=14461)")



pdf(width = 14, height = 12, "Plots/cvd_population_cvd_outcome.pdf")
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
# fourth plot
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


breaks_hf <- c(-5, -3, 0, 3, 5)

group.hf.dataset <- group_values(data = hf.dataset,
                                  variable = "effects",
                                  breaks = breaks_hf) %>%
  drop_na(intervals)


# Heart failure survival analysis

## Propensity score matching

# maximum number of deciles being tested
quantiles <- length(levels(group.hf.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_hf_hf_stan_psm_1_1 <- vector()

matching_hf <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_hf_hf_psm_1_1_male <- vector("list", quantiles)

models_hf_hf_psm_1_1_female <- vector("list", quantiles)


if (class(try(
  
  predictions_hf_hf_stan_psm_1_1 <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    models_hf_hf_psm_1_1_male[[i]] <- coxph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass,
                                              data = group.hf.dataset.matched %>%
                                                slice(which(group.hf.dataset.matched$intervals == levels(group.hf.dataset.matched$intervals)[i])) %>%
                                                filter(sex == "Male"))
    
    
    predictions_hf_hf_stan_psm_1_1 <- rbind(predictions_hf_hf_stan_psm_1_1, cbind(mean = models_hf_hf_psm_1_1_male[[i]]$coefficients[1],
                                                                                      lci = confint(models_hf_hf_psm_1_1_male[[i]])[1],
                                                                                      uci = confint(models_hf_hf_psm_1_1_male[[i]])[2],
                                                                                      sex = "Male",
                                                                                      intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
    
    models_hf_hf_psm_1_1_female[[i]] <- coxph(Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass,
                                                data = group.hf.dataset.matched %>%
                                                  slice(which(group.hf.dataset.matched$intervals == levels(group.hf.dataset.matched$intervals)[i])) %>%
                                                  filter(sex == "Female"))
    
    
    predictions_hf_hf_stan_psm_1_1 <- rbind(predictions_hf_hf_stan_psm_1_1, cbind(mean = models_hf_hf_psm_1_1_female[[i]]$coefficients[1],
                                                                                      lci = confint(models_hf_hf_psm_1_1_female[[i]])[1],
                                                                                      uci = confint(models_hf_hf_psm_1_1_female[[i]])[2],
                                                                                      sex = "Female",
                                                                                      intervals = levels(group.hf.dataset.matched$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_hf_hf_stan_psm_1_1 <- predictions_hf_hf_stan_psm_1_1 %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_hf_hf_stan_psm_1_1, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1.rds"))
  
}


## Propensity score matching + adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.hf.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_hf_hf_stan_psm_1_1_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.hf.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

matching_hf <- MatchIt::matchit(
  formula = formula(paste0("drugclass ~  agetx + t2dmduration + hba1cmonth + prehba1c + preegfr + prealt + drugline + ncurrtx + sex + preneuropathy + preretinopathy")),
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

models_hf_hf_psm_1_1_adjusted_male <- vector("list", quantiles)

models_hf_hf_psm_1_1_adjusted_female <- vector("list", quantiles)

if (class(try(

  predictions_hf_hf_stan_psm_1_1_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted.rds"))

  , silent = TRUE)) == "try-error") {

  for (i in mnumber) {

    # male data
    data.new <- group.hf.dataset.matched %>%
      slice(which(group.hf.dataset.matched$intervals == levels(group.hf.dataset.matched$intervals)[i])) %>%
      filter(sex == "Male")

    # vars to use
    vars_to_use <- breakdown_adjust[factors]

    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)

    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }

    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))



    models_hf_hf_psm_1_1_adjusted_male[[i]] <- coxph(as.formula(formula),
                                                       data = data.new)


    predictions_hf_hf_stan_psm_1_1_adjusted <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted, cbind(mean = models_hf_hf_psm_1_1_adjusted_male[[i]]$coefficients[1],
                                                                                                        lci = confint(models_hf_hf_psm_1_1_adjusted_male[[i]])[1,1],
                                                                                                        uci = confint(models_hf_hf_psm_1_1_adjusted_male[[i]])[1,2],
                                                                                                        sex = "Male",
                                                                                                        intervals = levels(group.hf.dataset.matched$intervals)[i]))

    # female data
    data.new <- group.hf.dataset.matched %>%
      slice(which(group.hf.dataset.matched$intervals == levels(group.hf.dataset.matched$intervals)[i])) %>%
      filter(sex == "Female")

    # vars to use
    vars_to_use <- breakdown_adjust[factors]

    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)

    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }

    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))

    models_hf_hf_psm_1_1_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                         data = data.new)


    predictions_hf_hf_stan_psm_1_1_adjusted <- rbind(predictions_hf_hf_stan_psm_1_1_adjusted, cbind(mean = models_hf_hf_psm_1_1_adjusted_female[[i]]$coefficients[1],
                                                                                                        lci = confint(models_hf_hf_psm_1_1_adjusted_female[[i]])[1,1],
                                                                                                        uci = confint(models_hf_hf_psm_1_1_adjusted_female[[i]])[1,2],
                                                                                                        sex = "Female",
                                                                                                        intervals = levels(group.hf.dataset.matched$intervals)[i]))




  }


  predictions_hf_hf_stan_psm_1_1_adjusted <- predictions_hf_hf_stan_psm_1_1_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))



  saveRDS(predictions_hf_hf_stan_psm_1_1_adjusted, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_psm_1_1_adjusted.rds"))

}


## Adjusted

# maximum number of deciles being tested
quantiles <- length(levels(group.hf.dataset[,"intervals"]))
# create lists with results
mnumber = c(1:quantiles)
predictions_hf_hf_stan_adjusted <- vector()

breakdown_adjust <- unique(c(variables_mu, variables_tau))
# categorical variables in breakdown
factors <- sapply(group.hf.dataset[,breakdown_adjust], is.factor)
factors["sex"] = FALSE

models_hf_hf_adjusted_male <- vector("list", quantiles)

models_hf_hf_adjusted_female <- vector("list", quantiles)



if (class(try(
  
  predictions_hf_hf_stan_adjusted <- readRDS(paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  for (i in mnumber) {
    
    # male data
    data.new <- group.hf.dataset %>%
      slice(which(group.hf.dataset$intervals == levels(group.hf.dataset$intervals)[i])) %>%
      filter(sex == "Male")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_hf_hf_adjusted_male[[i]] <- coxph(as.formula(formula),
                                               data = data.new)
    
    
    predictions_hf_hf_stan_adjusted <- rbind(predictions_hf_hf_stan_adjusted, cbind(mean = models_hf_hf_adjusted_male[[i]]$coefficients[1],
                                                                                        lci = confint(models_hf_hf_adjusted_male[[i]])[1,1],
                                                                                        uci = confint(models_hf_hf_adjusted_male[[i]])[1,2],
                                                                                        sex = "Male",
                                                                                        intervals = levels(group.hf.dataset$intervals)[i]))
    
    # female data
    data.new <- group.hf.dataset %>%
      slice(which(group.hf.dataset$intervals == levels(group.hf.dataset$intervals)[i])) %>%
      filter(sex == "Female")
    
    # vars to use
    vars_to_use <- breakdown_adjust[factors]
    
    one_category_factors <- sapply(data.new[, vars_to_use], function(x) length(unique(x)) == 1)
    
    if (!is_empty(vars_to_use[one_category_factors])) {
      # which variables in breakdown_adjust only have one category represented
      for (l in 1:length(breakdown_adjust[factors][one_category_factors])) {
        # position
        position <- which(vars_to_use == breakdown_adjust[factors][one_category_factors][l])
        # remove var from breakdown_adjust
        vars_to_use <- vars_to_use[-position]
      }
    }
    
    formula <- paste0("Surv(postdrug_hf_censtime_yrs, postdrug_hf_censvar) ~ drugclass + rcs(agetx, 3) + rcs(t2dmduration, 3) + rcs(hba1cmonth, 3) + rcs(prehba1c, 3) + rcs(preegfr, 3) + rcs(prealt, 3) + rcs(prebmi, 3) +", paste(vars_to_use, collapse = "+"))
    
    
    models_hf_hf_adjusted_female[[i]] <- coxph(as.formula(formula),
                                                 data = data.new)
    
    
    predictions_hf_hf_stan_adjusted <- rbind(predictions_hf_hf_stan_adjusted, cbind(mean = models_hf_hf_adjusted_female[[i]]$coefficients[1],
                                                                                        lci = confint(models_hf_hf_adjusted_female[[i]])[1,1],
                                                                                        uci = confint(models_hf_hf_adjusted_female[[i]])[1,2],
                                                                                        sex = "Female",
                                                                                        intervals = levels(group.hf.dataset$intervals)[i]))
    
    
    
    
  }
  
  
  predictions_hf_hf_stan_adjusted <- predictions_hf_hf_stan_adjusted %>%
    as.data.frame() %>%
    mutate(mean = as.numeric(mean),
           lci = as.numeric(lci),
           uci = as.numeric(uci))
  
  
  
  saveRDS(predictions_hf_hf_stan_adjusted, paste0(output_path, "/additional_outcomes/predictions_hf_hf_stan_adjusted.rds"))
  
}


#:--------------------
#:---- PLOTS
#:--------------------

## limits

# hf_axis_min <- plyr::round_any(ceiling(min(c(
#   predictions_hf_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
#   predictions_hf_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min(),
#   predictions_hf_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% min()),
#   1)
# ), 2, f = floor)

hf_axis_min = -4


# hf_axis_max <- plyr::round_any(ceiling(max(c(
#   predictions_hf_hf_stan_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
#   predictions_hf_hf_stan_psm_1_1_adjusted %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max(),
#   predictions_hf_hf_stan_psm_1_1 %>% select(c("mean","lci","uci")) %>% as.data.frame() %>% gather() %>% select(value) %>% unlist() %>% as.numeric() %>% max()),
#   1)
# ), 2, f = ceiling)

hf_axis_max = 4

#:---- PSM 1:1

plot_hf_hf_psm_1_1_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_psm_1_1 %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=850)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1145)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2756)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2020)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=425)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=276)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=7472)")



plot_hf_hf_psm_1_1_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_psm_1_1 %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=378)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=343)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1316)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1963)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1395)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=843)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=6238)")


#:---- PSM 1:1 + adjusted

plot_hf_hf_psm_1_1_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=850)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=1145)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=2756)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=2020)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=425)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=276)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=7472)")



plot_hf_hf_psm_1_1_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_psm_1_1_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=378)", 
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=343)", 
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=1316)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=1963)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1395)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=843)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Propensity score matching + adjusted",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=6238)")


#:---- Adjusted

plot_hf_hf_adjusted_male <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Male", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_adjusted %>%
    filter(sex == "Male") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=2105)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=3154)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=9281)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=7924)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=1574)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=833)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Male (n=24871)")



plot_hf_hf_adjusted_female <- rbind(
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on SGLT2i"),
  predictions_hf_hf_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(c(1:3)),
  cbind(mean = -10, lci = -10, uci = -10, sex = "Female", intervals = "Predicted HbA1c benefit on GLP1-RA"),
  predictions_hf_hf_stan_adjusted %>%
    filter(sex == "Female") %>%
    slice(-c(1:3))
) %>%
  as.data.frame() %>%
  mutate(mean = as.numeric(mean),
         lci = as.numeric(lci),
         intervals = ifelse(intervals == "(-20,-5]", ">5 mmol/mol (n=801)",
                            ifelse(intervals == "(-5,-3]", "3-5 mmol/mol (n=747)",
                                   ifelse(intervals == "(-3,0]", "0-3 mmol/mol (n=3183)",
                                          ifelse(intervals == "(0,3]", "0-3 mmol/mol (n=5481)",
                                                 ifelse(intervals == "(3,5]", "3-5 mmol/mol (n=4081)",
                                                        ifelse(intervals == "(5,40]", ">5 mmol/mol (n=2187)", intervals)))))),
         uci = as.numeric(uci)) %>%
  rename("lower" = "lci", "upper" = "uci", "labeltext" = "intervals") %>%
  forestplot(labeltext = labeltext,
             ci.vertices = TRUE,
             title = "Adjusted",
             clip = c(hf_axis_min, hf_axis_max),
             xticks = seq(hf_axis_min, hf_axis_max, 2),
             ci.vertices.height = 0.05,
             # xlog = TRUE,
             boxsize = .1,
             xlab = "log(Hazards Ratio)") %>%
  fp_add_header("Female (n=16480)")



pdf(width = 14, height = 12, "Plots/hf_population_hf_outcome.pdf")
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
# fourth plot
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


















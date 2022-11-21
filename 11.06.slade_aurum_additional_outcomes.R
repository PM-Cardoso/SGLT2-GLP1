####################
## Description:
##  - In this file we check:
##    - Weight change
##    - eGFR change
##    - Discontinuation
####################

library(tidyverse)
library(bcf)


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
dir.create(paste0(output_path, "/additional_outcomes"))

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
treatment_effects <- readRDS(paste0(output_path, "/response_model/patient_effects.rds"))



#:--------------------------------------------------------
# Weight change

## Read in data for weight
weight.dataset <- set_up_data_sglt2_glp1(dataset.type = "weight.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(w.change = postweight - preweight)


breaks_weight <- c(-8, -5, -3, 0, 3, 5, 8)

group.weight.dataset <- group_values(data = weight.dataset,
                                     variable = "effects",
                                     breaks = breaks_weight) %>%
  select(postweight, preweight, w.change, drugclass, intervals, prop.score) %>%
  drop_na()

# keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
prop_score <- 1 - group.weight.dataset$prop.score

# split predicted treatment effects into deciles
predicted_treatment_effect <- group.weight.dataset %>%
  plyr::ddply("intervals", dplyr::summarise,
              N = length(postweight))

# maximum number of deciles being tested
quantiles <- length(levels(group.weight.dataset[,"intervals"]))

# create lists with results
mnumber = c(1:quantiles)
models  <- as.list(1:quantiles)
obs <- vector(); lci <- vector(); uci <- vector(); intercept <- vector(); predictions <- data.frame(NULL)

# weights for SGLT2 Z = 1
sglt2.data <- group.weight.dataset %>%
  filter(drugclass == "SGLT2") %>%
  mutate(calc_prop = 1/(prop.score))

# weights for GLP1 Z = 0
glp1.data <- group.weight.dataset %>%
  filter(drugclass == "GLP1") %>%
  mutate(calc_prop = 1/(1-prop.score))

data.new <- rbind(sglt2.data, glp1.data)

# formula
formula <- "w.change ~ factor(drugclass) + preweight"

# iterate through deciles
for (i in mnumber) {
  # fit linear regression for decile
  models[[i]] <- lm(as.formula(formula),data=data.new,subset=data.new[,"intervals"]==levels(group.weight.dataset[,"intervals"])[i], weights = calc_prop)
  
  # collect treatment effect from regression
  obs <- append(obs,models[[i]]$coefficients[2])
  
  # collect intercept from regression
  intercept <- append(intercept,models[[i]]$coefficients[1])
  
  # calculate confidence intervals
  confint_all <- confint(models[[i]], levels=0.95)
  
  # collect lower bound CI
  lci <- append(lci,confint_all[2,1])
  
  # collect upper bound CI
  uci <- append(uci,confint_all[2,2])
  
  # predictions in SGLT2
  values <- predict(models[[i]], data.frame(drugclass = factor("SGLT2", levels = c("SGLT2", "GLP1")), preweight = mean(data.new$preweight, na.rm = TRUE)), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "SGLT2", intervals = levels(group.weight.dataset[,"intervals"])[i]))
  
  # predictions in GLP1
  values <- predict(models[[i]], data.frame(drugclass = factor("GLP1", levels = c("SGLT2", "GLP1")), preweight = mean(data.new$preweight, na.rm = TRUE)), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "GLP1", intervals = levels(group.weight.dataset[,"intervals"])[i]))
  
}

# join treatment effects for deciles in a data.frame
effects <- data.frame(predicted_treatment_effect,cbind(obs,lci,uci))

plot_weight_benefit <- effects %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = obs, ymin = lci, ymax = uci)) +
  coord_flip()

plot_weight <- predictions %>%
  mutate(fit = as.numeric(fit),
         lwr = as.numeric(lwr),
         upr = as.numeric(upr),
         intervals = factor(intervals, levels = levels(group.weight.dataset[,"intervals"]))) %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = fit, ymin = lwr, ymax = upr, colour = drugclass), position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = c("GLP1 benefit >8 mmol/mol",
                              "GLP1 benefit 5-8 mmol/mol",
                              "GLP1 benefit 3-5 mmol/mol",
                              "GLP1 benefit 0-3 mmol/mol",
                              "SGLT2 benefit 0-3 mmol/mol",
                              "SGLT2 benefit 3-5 mmol/mol",
                              "SGLT2 benefit 5-8 mmol/mol",
                              "SGLT2 benefit >8 mmol/mol")) +
  ylab("Observed weight change") +
  xlab("Predicted treatment effect") +
  ggtitle("Weight change in CPRD") +
  coord_flip() +
  scale_colour_manual(values = c("red", "#f1a340")) +
  theme(legend.position = "bottom",
        legend.title = element_blank())



#:--------------------------------------------------------
# eGFR change

## Read in data for eGFR
egfr.dataset <- set_up_data_sglt2_glp1(dataset.type = "egfr.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated")) %>%
  mutate(egfr.change = postegfr - preegfr)


breaks_egfr <- c(-8, -5, -3, 0, 3, 5, 8)

group.egfr.dataset <- group_values(data = egfr.dataset,
                                     variable = "effects",
                                     breaks = breaks_egfr)
select(postegfr, preegfr, w.change, drugclass, intervals, prop.score) %>%
  drop_na()

# keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
prop_score <- 1 - group.egfr.dataset$prop.score

# split predicted treatment effects into deciles
predicted_treatment_effect <- group.egfr.dataset %>%
  plyr::ddply("intervals", dplyr::summarise,
              N = length(postegfr)) %>%
  drop_na()

# maximum number of deciles being tested
quantiles <- length(levels(group.egfr.dataset[,"intervals"]))

# create lists with results
mnumber = c(1:quantiles)
models  <- as.list(1:quantiles)
obs <- vector(); lci <- vector(); uci <- vector(); intercept <- vector(); predictions <- data.frame(NULL)

# weights for SGLT2 Z = 1
sglt2.data <- group.egfr.dataset %>%
  filter(drugclass == "SGLT2") %>%
  mutate(calc_prop = 1/(prop.score))

# weights for GLP1 Z = 0
glp1.data <- group.egfr.dataset %>%
  filter(drugclass == "GLP1") %>%
  mutate(calc_prop = 1/(1-prop.score))

data.new <- rbind(sglt2.data, glp1.data)

# formula
formula <- "egfr.change ~ factor(drugclass) + preegfr"

# iterate through deciles
for (i in mnumber) {
  # fit linear regression for decile
  models[[i]] <- lm(as.formula(formula),data=data.new,subset=data.new[,"intervals"]==levels(group.egfr.dataset[,"intervals"])[i], weights = calc_prop)
  
  # collect treatment effect from regression
  obs <- append(obs,models[[i]]$coefficients[2])
  
  # collect intercept from regression
  intercept <- append(intercept,models[[i]]$coefficients[1])
  
  # calculate confidence intervals
  confint_all <- confint(models[[i]], levels=0.95)
  
  # collect lower bound CI
  lci <- append(lci,confint_all[2,1])
  
  # collect upper bound CI
  uci <- append(uci,confint_all[2,2])
  
  # predictions in SGLT2
  values <- predict(models[[i]], data.frame(drugclass = factor("SGLT2", levels = c("SGLT2", "GLP1")), preegfr = mean(data.new$preegfr, na.rm = TRUE)), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "SGLT2", intervals = levels(group.egfr.dataset[,"intervals"])[i]))
  
  # predictions in GLP1
  values <- predict(models[[i]], data.frame(drugclass = factor("GLP1", levels = c("SGLT2", "GLP1")), preegfr = mean(data.new$preegfr, na.rm = TRUE)), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "GLP1", intervals = levels(group.egfr.dataset[,"intervals"])[i]))
  
}

# join treatment effects for deciles in a data.frame
effects <- data.frame(predicted_treatment_effect,cbind(obs,lci,uci))

plot_egfr_benefit <- effects %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = obs, ymin = lci, ymax = uci)) +
  coord_flip()

plot_egfr <- predictions %>%
  mutate(fit = as.numeric(fit),
         lwr = as.numeric(lwr),
         upr = as.numeric(upr),
         intervals = factor(intervals, levels = levels(group.egfr.dataset[,"intervals"]))) %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = fit, ymin = lwr, ymax = upr, colour = drugclass), position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = c("GLP1 benefit >8 mmol/mol",
                              "GLP1 benefit 5-8 mmol/mol",
                              "GLP1 benefit 3-5 mmol/mol",
                              "GLP1 benefit 0-3 mmol/mol",
                              "SGLT2 benefit 0-3 mmol/mol",
                              "SGLT2 benefit 3-5 mmol/mol",
                              "SGLT2 benefit 5-8 mmol/mol",
                              "SGLT2 benefit >8 mmol/mol")) +
  ylab("Observed eGFR change") +
  xlab("Predicted treatment effect") +
  ggtitle("eGFR change in CPRD") +
  coord_flip() +
  scale_colour_manual(values = c("red", "#f1a340")) +
  theme(legend.position = "bottom",
        legend.title = element_blank())




#:--------------------------------------------------------
# Discontinuation
## Read in data for discontinuation
discontinuation.dataset <- set_up_data_sglt2_glp1(dataset.type = "discontinuation.dataset") %>%
  left_join(patient_prop_scores, by = c("patid", "pated")) %>%
  left_join(treatment_effects, by = c("patid", "pated"))












pdf(width = 10, file = "Plots/11.06.additional_outcomes.pdf")

patchwork::wrap_plots(list(plot_weight, 
                           plot_egfr +
                             theme(axis.title.y = element_blank(),
                                   axis.text.y = element_blank()))) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom")) # center title of full plot

dev.off()




















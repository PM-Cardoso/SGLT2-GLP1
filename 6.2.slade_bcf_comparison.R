####################
## Description:
##  - In this file we make a collected of plots comparing assessment
##      measurements of several models: 
##        - 3.2. bcf vs 3.1. grf
##        - 3.2. bcf vs 5.1. Complete/Routine model/No propensity score
##        - 3.2. bcf vs 5.1. Complete/Routine model/Propensity score
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")


library(tidyverse)
library(bcf)
library(grf)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)


## make directory for outputs
dir.create(paste0(output_path, "/Comparison"))

## make directory for outputs
dir.create("Plots")

###############################################################################
###############################################################################
######################### Read Data / Model In ################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))
# name: final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))

###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################

source("0.1.slade_functions.R")

############################# BCF
### Complete model of only routine data, no propensity score (n: 9866))
#############################

data_complete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass, ncurrtx, drugline, yrdrugstart, t2dmduration, agetx, malesex, Category, hba1cmonth, prebmi, prealt, egfr_ckdepi, prehba1cmmol
  ) %>%
  drop_na() # removed 1302


data_complete_routine_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass, ncurrtx, drugline, yrdrugstart, t2dmduration, agetx, malesex, Category, hba1cmonth, prebmi, prealt, egfr_ckdepi, prehba1cmmol
  ) %>%
  drop_na() # removed 804


# Produce a model matrix for fitting bcf

dataset_full <- rbind(data_complete_routine_dev, data_complete_routine_val)

dataset_model.matrix <- model.matrix(~posthba1c_final + drugclass + ncurrtx + drugline + yrdrugstart + t2dmduration + agetx +
                                       malesex + Category + hba1cmonth + prebmi + prealt + egfr_ckdepi + prehba1cmmol, dataset_full) %>%
  as.data.frame() %>%
  select(-`(Intercept)`) %>%
  mutate(drugclass = drugclassSGLT2) %>%
  select(-drugclassSGLT2)

# calculate a propensity score for the model

prop.score <- glm(drugclass ~ ncurrtx + drugline + t2dmduration + agetx + 
                    malesex + Category + hba1cmonth + prebmi + prealt + egfr_ckdepi + prehba1cmmol, family = binomial(link = "logit"), data = dataset_full[1:nrow(data_complete_routine_dev),])


dataset_full_bcf <- dataset_model.matrix %>%
  mutate_all(function(x) as.numeric(x)) %>%
  as.matrix()

# fit the model

post <- bcf::bcf(y = dataset_full_bcf[1:nrow(data_complete_routine_dev),1],
                 z = dataset_full_bcf[1:nrow(data_complete_routine_dev),19],
                 x_control = dataset_full_bcf[1:nrow(data_complete_routine_dev),-c(1,19)],
                 pihat = prop.score$fitted.values,
                 nburn = 1000,
                 nsim = 1000)

# collect average treatment effect

bcf.effects.dev <- cbind(mean = post$tau %>% colMeans()) %>%
  data.frame() %>%
  set_names(c("mean"))


##########




grf_model <- grf::causal_forest(X = dataset_model.matrix %>%
                                  slice(1:nrow(data_complete_routine_dev)) %>%
                                  select(-posthba1c_final, -drugclass),
                                Y = dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"],
                                W = dataset_model.matrix[1:nrow(data_complete_routine_dev), "drugclass"],
                                W.hat = prop.score$fitted.values)

#Dev
grf.effects.dev <- cbind(mean = grf_model$predictions) %>%
  data.frame() %>%
  set_names(c("mean"))


### BCF vs GRF

plot_effect_comparison_1 <- cbind(`BCF Effect` = bcf.effects.dev[,"mean"],
                                  `GRF Effect` = grf.effects.dev[,"mean"]) %>%
  as.data.frame() %>%
  ggplot() +
    theme_bw() +
    geom_point(aes(x = `BCF Effect`, y = `GRF Effect`)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    geom_smooth(aes(x = `BCF Effect`, y = `GRF Effect`), colour = "red", method = "lm") +
    xlab("BCF: Predicted CATE") +
    ylab("GRF: Predicted CATE") +
    ggtitle("BCF vs GRF")




### Complete model of only routine data, no propensity score

comp_routine_no_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_no_prop_effects_summary_dev.rds"))
    
plot_effect_comparison_2 <- as.data.frame(comp_routine_no_prop_effects_summary_dev) %>%
  cbind(`BCF Effect` = bcf.effects.dev[,"mean"]) %>%
  ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`, x = `BCF Effect`), colour = "grey") +
    geom_point(aes(x = `BCF Effect`, y = `mean`)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    geom_smooth(aes(x = `BCF Effect`, y = `mean`), colour = "red", method = "lm") +
    xlab("BCF: Predicted CATE") +
    ylab("BART: Predicted treatment Heterogeneity") +
    ggtitle("BCF vs BART (no propensity score)")
    

bart_comp_routine_no_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_no_prop.rds"))

# Dev
data_complete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_no_prop$X)
  ) %>% 
  drop_na()


# calculate effects
if (class(try(
  
  comp_routine_no_prop_effects_dev <- readRDS(paste0(output_path, "/Comparison/comp_routine_no_prop_effects_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_no_prop_effects_dev <- calc_effect(bart_comp_routine_no_prop, data_complete_routine_dev)
  
  saveRDS(comp_routine_no_prop_effects_dev, paste0(output_path, "/Comparison/comp_routine_no_prop_effects_dev.rds"))
  
}

# To calculate error at each iteration, we do BART effects - BCF effects. 
#   - A positive value: BART is over-estimating
#   - A negative value: BART is under-estimating

error_comparison_2 <- comp_routine_no_prop_effects_dev - t(post$tau)

plot_error_comparison_2 <- cbind(
    `5%` = apply(error_comparison_2, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(error_comparison_2, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(error_comparison_2, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(error_comparison_2, MARGIN = 1, function(x) mean(c(x))),
    `BCF Effect` = bcf.effects.dev[,"mean"]
  ) %>%
  as.data.frame() %>%
  ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`, x = `BCF Effect`), colour = "grey") +
    geom_point(aes(x = `BCF Effect`, y = `mean`)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = viridis::viridis(1, begin = 0.6), lwd = 0.75) +
    geom_smooth(aes(x = `BCF Effect`, y = `mean`), colour = "red", method = "lm") +
    xlab("BCF: Predicted CATE") +
    ylab("Treatment Effect difference") +
    ggtitle("BCF vs BART (no propensity score)")
  



### Complete model of only routine data, propensity score

comp_routine_prop_effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/Assessment/comp_routine_prop_effects_summary_dev.rds"))

plot_effect_comparison_3 <- as.data.frame(comp_routine_prop_effects_summary_dev) %>%
  cbind(`BCF Effect` = bcf.effects.dev[,"mean"]) %>%
  ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`, x = `BCF Effect`), colour = "grey") +
    geom_point(aes(x = `BCF Effect`, y = `mean`)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    geom_smooth(aes(x = `BCF Effect`, y = `mean`), colour = "red", method = "lm") +
    xlab("BCF: Predicted CATE") +
    ylab("BART: Predicted treatment Heterogeneity") +
    ggtitle("BCF vs BART (propensity score)")


bart_comp_routine_prop <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop.rds"))

bart_comp_routine_prop_model <- readRDS(paste0(output_path, "/Model_fit/bart_comp_routine_prop_model.rds"))


# Dev
data_complete_routine_prop_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    colnames(bart_comp_routine_prop_model$X)[which(colnames(bart_comp_routine_prop_model$X) != "prop_score")]
  ) %>% 
  drop_na() %>%
  cbind(prop_score = bart_comp_routine_prop$p_hat_train)


# calculate effects
if (class(try(
  
  comp_routine_prop_effects_dev <- readRDS(paste0(output_path, "/Comparison/comp_routine_prop_effects_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  comp_routine_prop_effects_dev <- calc_effect(bart_comp_routine_prop_model, data_complete_routine_prop_dev)
  
  saveRDS(comp_routine_prop_effects_dev, paste0(output_path, "/Comparison/comp_routine_prop_effects_dev.rds"))
  
}

# To calculate error at each iteration, we do BART effects - BCF effects. 
#   - A positive value: BART is over-estimating
#   - A negative value: BART is under-estimating

error_comparison_3 <- comp_routine_prop_effects_dev - t(post$tau)

plot_error_comparison_3 <- cbind(
    `5%` = apply(error_comparison_3, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(error_comparison_3, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(error_comparison_3, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(error_comparison_3, MARGIN = 1, function(x) mean(c(x))),
    `BCF Effect` = bcf.effects.dev[,"mean"]
  ) %>%
  as.data.frame() %>%
  ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`, x = `BCF Effect`), colour = "grey") +
    geom_point(aes(x = `BCF Effect`, y = `mean`)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = viridis::viridis(1, begin = 0.6), lwd = 0.75) +
    geom_smooth(aes(x = `BCF Effect`, y = `mean`), colour = "red", method = "lm") +
    xlab("BCF: Predicted CATE") +
    ylab("Treatment Effect difference") +
    ggtitle("BCF vs BART (propensity score)")



pdf(file = "Plots/6.2.comparison_bcf_bart.pdf")
plot_effect_comparison_1
plot_effect_comparison_2
plot_error_comparison_2
plot_effect_comparison_3
plot_error_comparison_3
dev.off()




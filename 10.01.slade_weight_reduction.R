####################
## Description:
##  - In this file we find a way to know weight change expected:
##      1. fit a BART model for predicting weight reduction in a therapy.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")

library(tidyverse)
library(bartMachine)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path, "/Weight_reduction"))

## make directory for outputs
dir.create("Plots")



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.all.extra.vars
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_allcohort.Rda"))

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

# name:final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))


###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################

source("0.1.slade_functions.R")


resid_plot <- function(pred_dev, pred_val, title) {
  ##### Imput variables
  # pred_dev - predicted/observed values for development dataset
  # pred_val - predicted/observed values for validation dataset
  # title - plot title
  
  # Full grid plot
  cowplot::plot_grid(
    
    # title
    cowplot::ggdraw() +
      cowplot::draw_label(title)
    
    ,
    
    cowplot::plot_grid(
      
      # Plot of predicted vs observed for development dataset
      pred_dev %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
        geom_point(aes(x = orig, y = mean)) +
        geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        xlab("Observed Post Weight") +
        ylab("Predicted Post Weight")
      
      ,
      
      # Plot of predicted vs observed for validation dataset
      pred_val %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
        geom_point(aes(x = orig, y = mean)) +
        geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        xlab("Observed Post Weight") +
        ylab("Predicted Post Weight")
      
      ,
      
      # Plot of standardised residuals for development dataset
      pred_dev %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
        geom_point(aes(x = mean, y = std.resid)) +
        geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        stat_smooth(aes(x = mean, y = std.resid)) +
        xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
        ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
        xlab("Average Predicted  Post Weight") +
        ylab("Standardised Residuals")
      
      ,
      
      # Plot of standardised residuals for validation dataset
      pred_val %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
        geom_point(aes(x = mean, y = std.resid)) +
        geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        stat_smooth(aes(x = mean, y = std.resid)) +
        xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
        ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
        xlab("Average Predicted Post Weight") +
        ylab("Standardised Residuals")
      
      , ncol = 2, nrow = 2, labels = c("A", "B", "", "")
      
    )
    
    , ncol = 1, nrow = 2, rel_heights = c(0.1,1)
    
  )
  
}


###############################################################################
###############################################################################
########################### Model Fitting #####################################
###############################################################################
###############################################################################
##
## Initially we start with 9866 individuals
##

dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi, postweight6m)) %>%
  drop_na(postweight6m) # drop 2846


## Fit model using all the available in treatment selection BART model variables to estimate weight change
if (class(try(
  
  bart_model_weight <- readRDS(paste0(output_path, "/Weight_reduction/bart_model_weight.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model_weight <- bartMachine::bartMachine(X = dataset.dev %>%
                                           select(drugclass,
                                                  preweight,
                                                  # below is vars from both variable selections
                                                  egfr_ckdepi,
                                                  hba1cmonth,
                                                  prealt,
                                                  prehba1cmmol,
                                                  score.excl.mi,
                                                  # below is vars from BART variable selection
                                                  Category,
                                                  drugline,
                                                  ncurrtx,
                                                  yrdrugstart,
                                                  # below is the vars from grf variable selection
                                                  agetx,
                                                  malesex,
                                                  prehdl,
                                                  # prebmi, # used preweight instead of prebmi
                                                  prebil,
                                                  preplatelets,
                                                  t2dmduration,
                                                  prealb,
                                                  presys,
                                                  preast),
                                         # y = dataset.dev[,"weight_change"],
                                         y = dataset.dev[,"postweight6m"],
                                         use_missing_data = TRUE,
                                         impute_missingness_with_rf_impute = FALSE,
                                         impute_missingness_with_x_j_bar_for_lm = TRUE,
                                         num_trees = 200,
                                         num_burn_in = 3000,
                                         num_iterations_after_burn_in = 1000,
                                         serialize = TRUE)
  
  saveRDS(bart_model_weight, paste0(output_path, "/Weight_reduction/bart_model_weight.rds"))
  
}


########
### Validation of model
########

# Dev

data_dev <- dataset.dev %>% 
  select(c(patid, pateddrug, postweight6m, 
           colnames(bart_model_weight$X)))


## Get posteriors
if (class(try(
  
  posteriors_dev <- readRDS(paste0(output_path, "/Weight_reduction/posteriors_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_dev <- bartMachine::bart_machine_get_posterior(bart_model_weight, data_dev %>%
                                                              select(
                                                                colnames(bart_model_weight$X)
                                                              ))
  saveRDS(posteriors_dev, paste0(output_path, "/Weight_reduction/posteriors_dev.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  cred_pred_dev <- readRDS(paste0(output_path, "/Weight_reduction/cred_pred_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  cred_pred_dev <- calc_resid(data_dev, posteriors_dev, "postweight6m")
  
  saveRDS(cred_pred_dev, paste0(output_path, "/Weight_reduction/cred_pred_dev.rds"))
  
}


# Val
data_val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi, postweight6m)) %>%
  drop_na(postweight6m) %>%  # drop 1934
  select(c(patid, pateddrug, postweight6m, 
           colnames(bart_model_weight$X)))


### Get posteriors
if (class(try(
  
  posteriors_val <- readRDS(paste0(output_path, "/Weight_reduction/posteriors_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  posteriors_val <- bartMachine::bart_machine_get_posterior(bart_model_weight, data_val %>%
                                                              select(
                                                                colnames(bart_model_weight$X)
                                                              ))
  saveRDS(posteriors_val, paste0(output_path, "/Weight_reduction/posteriors_val.rds"))
  
  
}

### residuals calculation
if (class(try(
  
  cred_pred_val <- readRDS(paste0(output_path, "/Weight_reduction/cred_pred_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  interim_dataset <- data_val %>%
    mutate(posthba1c_final = postweight6m) %>%
    select(-postweight6m)
  
  cred_pred_val <- calc_resid(interim_dataset, posteriors_val)
  
  saveRDS(cred_pred_val, paste0(output_path, "/Weight_reduction/cred_pred_val.rds"))
  
}



plot_residuals <- resid_plot(cred_pred_dev, cred_pred_val, "Residuals of Model 4")

## calculate best drug for weight reduction

# calculate best weight drug
if (class(try(
  
  effects_summary_dev <- readRDS(paste0(output_path, "/Weight_reduction/effects_summary_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_dev <- calc_effect_summary(bart_model_weight, data_dev)
  
  saveRDS(effects_summary_dev, paste0(output_path, "/Weight_reduction/effects_summary_dev.rds"))
  
}


# calculate best weight drug
if (class(try(
  
  effects_summary_val <- readRDS(paste0(output_path, "/Weight_reduction/effects_summary_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  effects_summary_val <- calc_effect_summary(bart_model_weight, data_val)
  
  saveRDS(effects_summary_val, paste0(output_path, "/Weight_reduction/effects_summary_val.rds"))
  
}


## plot effects validation

predicted_observed_dev <- data_dev %>%
  cbind(hba1c_diff = effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

predicted_observed_val <- data_val %>%
  cbind(hba1c_diff = effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))


##############
# Validating ATE
if (class(try(
  
  ATE_validation_dev <- readRDS(paste0(output_path, "/Weight_reduction/ATE_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_dev <- calc_ATE_validation(predicted_observed_dev, "postweight6m")
  
  saveRDS(ATE_validation_dev, paste0(output_path, "/Weight_reduction/ATE_validation_dev.rds"))
  
}

plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -8, 10)


if (class(try(
  
  ATE_validation_val <- readRDS(paste0(output_path, "/Weight_reduction/ATE_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_validation_val <- calc_ATE_validation(predicted_observed_val, "postweight6m")
  
  saveRDS(ATE_validation_val, paste0(output_path, "/Weight_reduction/ATE_validation_val.rds"))
  
}

plot_ATE_val <- ATE_plot(ATE_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 8)

plot_ATE <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation lm(weght~drugclass+prop_score)")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev, plot_ATE_val, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score matching
if (class(try(
  
  ATE_matching_validation_dev <- readRDS(paste0(output_path, "/Weight_reduction/ATE_matching_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_dev, "postweight6m")
  
  saveRDS(ATE_matching_validation_dev, paste0(output_path, "/Weight_reduction/ATE_matching_validation_dev.rds"))
  
}

plot_ATE_dev_prop_score <- ATE_plot(ATE_matching_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10)

if (class(try(
  
  ATE_matching_validation_val <- readRDS(paste0(output_path, "/Weight_reduction/ATE_matching_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_matching_validation_val <- calc_ATE_validation_prop_matching(predicted_observed_val, "postweight6m")
  
  saveRDS(ATE_matching_validation_val, paste0(output_path, "/Weight_reduction/ATE_matching_validation_val.rds"))
  
}

plot_ATE_val_prop_score <- ATE_plot(ATE_matching_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_ATE_prop_score_matching <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score matching")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score, plot_ATE_val_prop_score, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))


# Validation ATE prop score inverse weighting
if (class(try(
  
  ATE_weighting_validation_dev <- readRDS(paste0(output_path, "/Weight_reduction/ATE_weighting_validation_dev.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_dev <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_dev, "postweight6m")
  
  saveRDS(ATE_weighting_validation_dev, paste0(output_path, "/Weight_reduction/ATE_weighting_validation_dev.rds"))
  
}

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

if (class(try(
  
  ATE_weighting_validation_val <- readRDS(paste0(output_path, "/Weight_reduction/ATE_weighting_validation_val.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "postweight6m")
  
  saveRDS(ATE_weighting_validation_val, paste0(output_path, "/Weight_reduction/ATE_weighting_validation_val.rds"))
  
}

plot_ATE_val_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 8)

plot_ATE_prop_score_weighting <- cowplot::plot_grid(
  
  cowplot::ggdraw() +
    cowplot::draw_label("Effects validation prop score inverse weighting")
  
  ,
  
  cowplot::plot_grid(plot_ATE_dev_prop_score_weighting, plot_ATE_val_prop_score_weighting, ncol = 2, nrow = 1, labels = c("A", "B"))
  
  , nrow = 2, ncol = 1, rel_heights = c(0.1, 1))




## Plot preference of therapy for weight change

plot_effect_1 <- -1*effects_summary_dev %>%
  as.data.frame() %>%
  select(mean)
plot_effect_1 <- plot_effect_1 %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())


plot_effect_2 <- -1*effects_summary_val %>%
  as.data.frame() %>%
  select(mean)
plot_effect_2 <- plot_effect_2 %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# Plot preference of therapy for weight change by gender

effects_summary_dev_male <- effects_summary_dev %>%
  cbind(malesex = data_dev$malesex) %>%
  filter(malesex == 1)

effects_summary_dev_female <- effects_summary_dev %>%
  cbind(malesex = data_dev$malesex) %>%
  filter(malesex == 0)


plot_effect_1_male <- -1*effects_summary_dev_male %>%
  as.data.frame() %>%
  select(mean)
plot_effect_1_male <- plot_effect_1_male %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Male",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

plot_effect_1_female <- -1*effects_summary_dev_female %>%
  as.data.frame() %>%
  select(mean)
plot_effect_1_female <- plot_effect_1_female %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Female",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())


effects_summary_val_male <- effects_summary_val %>%
  cbind(malesex = data_val$malesex) %>%
  filter(malesex == 1)

effects_summary_val_female <- effects_summary_val %>%
  cbind(malesex = data_val$malesex) %>%
  filter(malesex == 0)

plot_effect_2_male <- -1*effects_summary_val_male %>%
  as.data.frame() %>%
  select(mean)
plot_effect_2_male <- plot_effect_2_male %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Male",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

plot_effect_2_female <- -1*effects_summary_val_female %>%
  as.data.frame() %>%
  select(mean)
plot_effect_2_female <- plot_effect_2_female %>%
  mutate(above=ifelse(mean > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-3.5,7.5,by=0.25)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Female",x="Post Weight", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())



###############################
####### Description of weight change for deciles of model

bart_model_final <- readRDS(paste0(output_path, "/Final_model/model_5/bart_model_final.rds"))

data_dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X))) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, preweight, postweight6m))


data_val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X))) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, preweight, postweight6m))


effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/model_5/Assessment/effects_summary_dev.rds"))

effects_summary_val <- readRDS(paste0(output_path, "/Final_model/model_5/Assessment/effects_summary_val.rds"))


#####################

group_values <- function(data, variable, breaks) {
  ### Input variables
  # data: dataset used in splitting
  # variable: variable with values to be split
  # breaks: break points between values
  
  # stop in case 'variable' is not included in 'data'
  if (is.null(data[, variable])) {stop("'variable' not included in 'data'")}
  
  # include extra values so that extremes are included
  breaks.full <- c(breaks, floor(min(data[,variable])), ceiling(max(data[,variable])))
  
  new.data <- data %>%
    cbind(intervals = cut(data[, variable], breaks = breaks.full))
  
  return(new.data)
}

# description of weight change by levels:
## SGLT2:
### < -8 mmol
### -5 - -8 mmol
### -3 - -5 mmol
### 0 - -3 mmol
## GLP1:
### 0 - 3 mmol
### 3 - 5 mmol
### 5 - 8 mmol
### > 8 mmol

# breaks
breaks = c(-5, -3, 0, 3, 5)

######
####

dataset_breakdown_dev <- data_dev %>%
  mutate(weight.change = postweight6m - preweight) %>%
  select(weight.change, drugclass) %>%
  cbind(hba1c_diff = effects_summary_dev$mean)

#########
######## # Unadjusted for indication bias
#########

dataset_intervals_dev <- group_values(data = dataset_breakdown_dev, 
                                  variable = "hba1c_diff", 
                                  breaks = breaks) %>%
  select(-hba1c_diff) %>%
  group_by(intervals, drugclass) %>%
  mutate(mean.values = quantile(weight.change, probs = c(0.5), na.rm = TRUE),
         low.values = quantile(weight.change, probs = c(0.25), na.rm = TRUE),
         high.values = quantile(weight.change, probs = c(0.75), na.rm = TRUE)) %>%
  select(-weight.change) %>%
  unique()


plot_weight_change_dev <- dataset_intervals_dev %>%
  ggplot() +
  geom_pointrange(aes(x = intervals,
                      y = mean.values,
                      ymin = low.values, 
                      ymax = high.values, 
                      colour = drugclass), 
                  position=position_dodge(width=0.5)) +
  ylab("Weight change kg (median [IQR])") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")

##

dataset_breakdown_val <- data_val %>%
  mutate(weight.change = postweight6m - preweight) %>%
  select(weight.change, drugclass) %>%
  cbind(hba1c_diff = effects_summary_val$mean)


dataset_intervals_val <- group_values(data = dataset_breakdown_val, 
                                      variable = "hba1c_diff", 
                                      breaks = breaks) %>%
  select(-hba1c_diff) %>%
  group_by(intervals, drugclass) %>%
  mutate(mean.values = quantile(weight.change, probs = c(0.5), na.rm = TRUE),
         low.values = quantile(weight.change, probs = c(0.25), na.rm = TRUE),
         high.values = quantile(weight.change, probs = c(0.75), na.rm = TRUE)) %>%
  select(-weight.change) %>%
  unique()


plot_weight_change_val <- dataset_intervals_val %>%
  ggplot() +
  geom_pointrange(aes(x = intervals,
                      y = mean.values,
                      ymin = low.values, 
                      ymax = high.values, 
                      colour = drugclass), 
                  position=position_dodge(width=0.5)) +
  ylab("Weight change kg (median [IQR])") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")


plot_weight_change_unadjusted <- patchwork::wrap_plots(
  # Development
  plot_weight_change_dev,
  # Validation
  plot_weight_change_val
) + patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')

#########
######## # Adjusted for indication bias
#########

# dev

# prop score model for indication bias
# load all data for range of variable values; name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# extracting selected variables for individuals in dataset
data.new <- data_dev %>%
  select(patid, pateddrug) %>%
  left_join(final.all.extra.vars %>%
              select(patid, 
                     pateddrug,
                     drugclass,
                     # prebmi,
                     preweight,
                     t2dmduration,
                     prealb,
                     egfr_ckdepi,
                     drugline,
                     prehba1cmmol,
                     ncurrtx,
                     score.excl.mi,
                     Category), by = c("patid", "pateddrug"))

# fit propensity model with the variables that influence therapy indication
prop_model <- bartMachine::bartMachine(X = data.new %>%
                                         select(# prebmi,
                                                preweight,
                                                t2dmduration,
                                                prealb,
                                                egfr_ckdepi,
                                                drugline,
                                                prehba1cmmol,
                                                ncurrtx,
                                                score.excl.mi,
                                                Category),
                                       y = data.new[,"drugclass"],
                                       use_missing_data = TRUE,
                                       impute_missingness_with_rf_impute = FALSE,
                                       impute_missingness_with_x_j_bar_for_lm = TRUE,
                                       num_trees = 200,
                                       num_burn_in = 1000,
                                       num_iterations_after_burn_in = 200)

# keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
prop_score <- 1 - prop_model$p_hat_train


dataset_intervals_dev <- group_values(data = dataset_breakdown_dev, 
                                      variable = "hba1c_diff", 
                                      breaks = breaks)

# split predicted treatment effects into deciles
predicted_treatment_effect <- dataset_intervals_dev %>%
  plyr::ddply("intervals", dplyr::summarise,
              N = length(hba1c_diff),
              hba1c_diff.pred = mean(hba1c_diff),
              weight.change.pred = mean(weight.change, na.rm = TRUE))

# maximum number of deciles being tested
quantiles <- levels(dataset_intervals_dev$intervals)

# create lists with results
mnumber = c(1:length(quantiles))
models  <- as.list(1:length(quantiles))
hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector(); intercept.obs <- vector(); predictions <- data.frame(NULL)

# join dataset and propensity score
data.new <- dataset_intervals_dev %>%
  cbind(calc_prop = prop_score)

# weights for SGLT2 Z = 1
sglt2.data <- data.new %>%
  filter(drugclass == "SGLT2") %>%
  mutate(calc_prop = 1/(calc_prop))

# weights for GLP1 Z = 0
glp1.data <- data.new %>%
  filter(drugclass == "GLP1") %>%
  mutate(calc_prop = 1/(1-calc_prop))

data.new <- rbind(sglt2.data, glp1.data)

variable = "weight.change"
# formula
formula <- paste0(variable, " ~ factor(drugclass)")


# iterate through deciles
for (i in mnumber) {
  # fit linear regression for decile
  models[[i]] <- lm(as.formula(formula),data=data.new,subset=intervals==quantiles[i], weights = calc_prop)
  
  # collect treatment effect from regression
  hba1c_diff.obs <- append(hba1c_diff.obs,models[[i]]$coefficients[2])
  
  # collect intercept from regression
  intercept.obs <- append(intercept.obs,models[[i]]$coefficients[1])
  
  # calculate confidence intervals
  confint_all <- confint(models[[i]], levels=0.95)
  
  # collect lower bound CI
  lower.obs <- append(lower.obs,confint_all[2,1])
  
  # collect upper bound CI
  upper.obs <- append(upper.obs,confint_all[2,2])
  
  # predictions in SGLT2
  values <- predict(models[[i]], data.frame(drugclass = factor("SGLT2", levels = c("SGLT2", "GLP1"))), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "SGLT2", intervals = quantiles[i]))
  
  # predictions in GLP1
  values <- predict(models[[i]], data.frame(drugclass = factor("GLP1", levels = c("SGLT2", "GLP1"))), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "GLP1", intervals = quantiles[i]))
  
}

# join treatment effects for deciles in a data.frame
effects <- data.frame(predicted_treatment_effect,cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>%
  dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)

# returned list with fitted propensity model + decile treatment effects  
t <- list(prop_model = prop_model, effects = effects)

plot_weight_benefit_dev <- t[["effects"]] %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = hba1c_diff.obs, ymin = lower.obs, ymax = upper.obs)) +
  coord_flip()

plot_weight_dev <- predictions %>%
  mutate(fit = as.numeric(fit),
         lwr = as.numeric(lwr),
         upr = as.numeric(upr),
         intervals = factor(intervals, levels = quantiles)) %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = fit, ymin = lwr, ymax = upr, colour = drugclass), position = position_dodge(width = 0.5)) +
  ylim(-6, 0) +
  coord_flip() 
  
# Val

# prop score model for indication bias
# load all data for range of variable values; name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# extracting selected variables for individuals in dataset
data.new <- data_val %>%
  select(patid, pateddrug) %>%
  left_join(final.all.extra.vars %>%
              select(patid, 
                     pateddrug,
                     drugclass,
                     # prebmi,
                     preweight,
                     t2dmduration,
                     prealb,
                     egfr_ckdepi,
                     drugline,
                     prehba1cmmol,
                     ncurrtx,
                     score.excl.mi,
                     Category), by = c("patid", "pateddrug"))

# fit propensity model with the variables that influence therapy indication
prop_model <- bartMachine::bartMachine(X = data.new %>%
                                         select(# prebmi,
                                           preweight,
                                           t2dmduration,
                                           prealb,
                                           egfr_ckdepi,
                                           drugline,
                                           prehba1cmmol,
                                           ncurrtx,
                                           score.excl.mi,
                                           Category),
                                       y = data.new[,"drugclass"],
                                       use_missing_data = TRUE,
                                       impute_missingness_with_rf_impute = FALSE,
                                       impute_missingness_with_x_j_bar_for_lm = TRUE,
                                       num_trees = 200,
                                       num_burn_in = 1000,
                                       num_iterations_after_burn_in = 200)

# keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
prop_score <- 1 - prop_model$p_hat_train


dataset_intervals_val <- group_values(data = dataset_breakdown_val, 
                                      variable = "hba1c_diff", 
                                      breaks = breaks)

# split predicted treatment effects into deciles
predicted_treatment_effect <- dataset_intervals_val %>%
  plyr::ddply("intervals", dplyr::summarise,
              N = length(hba1c_diff),
              hba1c_diff.pred = mean(hba1c_diff),
              weight.change.pred = mean(weight.change, na.rm = TRUE))

# maximum number of deciles being tested
quantiles <- levels(dataset_intervals_val$intervals)

# create lists with results
mnumber = c(1:length(quantiles))
models  <- as.list(1:length(quantiles))
hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector(); intercept.obs <- vector(); predictions <- data.frame(NULL)

# join dataset and propensity score
data.new <- dataset_intervals_val %>%
  cbind(calc_prop = prop_score)

# weights for SGLT2 Z = 1
sglt2.data <- data.new %>%
  filter(drugclass == "SGLT2") %>%
  mutate(calc_prop = 1/(calc_prop))

# weights for GLP1 Z = 0
glp1.data <- data.new %>%
  filter(drugclass == "GLP1") %>%
  mutate(calc_prop = 1/(1-calc_prop))

data.new <- rbind(sglt2.data, glp1.data)

variable = "weight.change"
# formula
formula <- paste0(variable, " ~ factor(drugclass)")


# iterate through deciles
for (i in mnumber) {
  # fit linear regression for decile
  models[[i]] <- lm(as.formula(formula),data=data.new,subset=intervals==quantiles[i], weights = calc_prop)
  
  # collect treatment effect from regression
  hba1c_diff.obs <- append(hba1c_diff.obs,models[[i]]$coefficients[2])
  
  # collect intercept from regression
  intercept.obs <- append(intercept.obs,models[[i]]$coefficients[1])
  
  # calculate confidence intervals
  confint_all <- confint(models[[i]], levels=0.95)
  
  # collect lower bound CI
  lower.obs <- append(lower.obs,confint_all[2,1])
  
  # collect upper bound CI
  upper.obs <- append(upper.obs,confint_all[2,2])
  
  # predictions in SGLT2
  values <- predict(models[[i]], data.frame(drugclass = factor("SGLT2", levels = c("SGLT2", "GLP1"))), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "SGLT2", intervals = quantiles[i]))
  
  # predictions in GLP1
  values <- predict(models[[i]], data.frame(drugclass = factor("GLP1", levels = c("SGLT2", "GLP1"))), interval = "confidence")
  
  predictions <- rbind(predictions, cbind(values, drugclass = "GLP1", intervals = quantiles[i]))
  
}

# join treatment effects for deciles in a data.frame
effects <- data.frame(predicted_treatment_effect,cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>%
  dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)

# returned list with fitted propensity model + decile treatment effects  
t <- list(prop_model = prop_model, effects = effects)

plot_weight_benefit_val <- t[["effects"]] %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = hba1c_diff.obs, ymin = lower.obs, ymax = upper.obs)) +
  coord_flip()

plot_weight_val <- predictions %>%
  mutate(fit = as.numeric(fit),
         lwr = as.numeric(lwr),
         upr = as.numeric(upr),
         intervals = factor(intervals, levels = quantiles)) %>%
  ggplot() +
  geom_pointrange(aes(x = intervals, y = fit, ymin = lwr, ymax = upr, colour = drugclass), position = position_dodge(width = 0.5)) +
  ylim(-6, 0) +
  coord_flip() 


## combining plots

plot_weight_change_adjusted_dev <- patchwork::wrap_plots(
  # Development
  plot_weight_benefit_dev,
  # Validation
  plot_weight_dev
) + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')


plot_weight_change_adjusted_val <- patchwork::wrap_plots(
  # Development
  plot_weight_benefit_val,
  # Validation
  plot_weight_val
) + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')





pdf(file = "Plots/10.1.weight_change_plots.pdf")
plot_residuals
cowplot::plot_grid(
  plot_effect_1,
  plot_effect_2,
  ncol = 2, nrow = 1, labels = c("A", "B")
)

cowplot::plot_grid(
  
  cowplot::plot_grid(plot_effect_1_male, plot_effect_1_female, ncol = 2, nrow = 1)
  
  ,
  
  cowplot::plot_grid(plot_effect_2_male, plot_effect_2_female, ncol = 2, nrow = 1)
  
  , nrow = 2, ncol = 1, labels = c("A", "B")
)
plot_ATE
plot_ATE_prop_score_matching
plot_ATE_prop_score_weighting
plot_weight_change_unadjusted
plot_weight_change_adjusted_dev
plot_weight_change_adjusted_val
dev.off()


###############################
###### Variable selection

# Perform variable selection in the weight mdoel

if (class(try(
  
  bart_var_selection <- readRDS(paste0(output_path, "/Weight_reduction/bart_var_selection.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_var_selection <- bartMachine::var_selection_by_permute(bart_machine = bart_model_weight,
                                                              num_reps_for_avg = 15,
                                                              num_permute_samples = 50,
                                                              num_trees_for_permute = 20)
  
  saveRDS(bart_var_selection, paste0(output_path, "/Weight_reduction/bart_var_selection.rds"))
  
}



###############################
####### Description of weight change for deciles of model

bart_model_final <- readRDS(paste0(output_path, "/Final_model/model_5/bart_model_final.rds"))

data_dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X))) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, preweight, postweight6m))


data_val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X))) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, preweight, postweight6m))


effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/model_5/Assessment/effects_summary_dev.rds"))

effects_summary_val <- readRDS(paste0(output_path, "/Final_model/model_5/Assessment/effects_summary_val.rds"))


#####################

group_values <- function(data, variable, breaks) {
  ### Input variables
  # data: dataset used in splitting
  # variable: variable with values to be split
  # breaks: break points between values
  
  # stop in case 'variable' is not included in 'data'
  if (is.null(data[, variable])) {stop("'variable' not included in 'data'")}
  
  # include extra values so that extremes are included
  breaks.full <- c(breaks, floor(min(data[,variable])), ceiling(max(data[,variable])))
  
  new.data <- data %>%
    cbind(intervals = cut(data[, variable], breaks = breaks.full))
  
  return(new.data)
}

# description of weight change by levels:
## SGLT2:
### < -8 mmol
### -5 - -8 mmol
### -3 - -5 mmol
### 0 - -3 mmol
## GLP1:
### 0 - 3 mmol
### 3 - 5 mmol
### 5 - 8 mmol
### > 8 mmol
#
# breaks
# breaks = c(-8, -5, -3, 0, 3, 5, 8)
#
# dataset_breakdown <- data_dev %>%
#   mutate(weight.change = postweight6m - preweight) %>%
#   select(weight.change) %>%
#   cbind(hba1c_diff = effects_summary_dev$mean)
# 
# 
# dataset_intervals <- group_values(data = dataset_breakdown, 
#                                   variable = "hba1c_diff", 
#                                   breaks = breaks) %>%
#   select(-hba1c_diff) %>%
#   group_by(intervals) %>%
#   mutate(mean.values = mean(weight.change, na.rm = TRUE),
#          low.values = quantile(weight.change, probs = c(0.05), na.rm = TRUE),
#          high.values = quantile(weight.change, probs = c(0.95), na.rm = TRUE)) %>%
#   select(-weight.change) %>%
#   unique()
#
######
####
# 
# ######
# ####
# 
# dataset_new <- data_dev %>%
#   mutate(weight.change = postweight6m - preweight) %>%
#   cbind(hba1c_diff = effects_summary_dev$mean)
# 
# 
# prop_score_model <- calc_ATE_prop_score(dataset_breakdown)
# 
# dataset_breakdown <- dataset_new %>%
#   cbind(prop_score = 1 - prop_score_model$p_hat_train) %>%
#   select(weight.change, hba1c_diff, prop_score, drugclass)
# 
# 
# 
# 
# dataset_intervals <- group_values(data = dataset_breakdown, 
#                                   variable = "hba1c_diff", 
#                                   breaks = breaks) %>%
#   select(-hba1c_diff) %>%
#   filter(intervals == "(5,8]") %>%
#   drop_na()
# 
# # weights for SGLT2 Z = 1
# sglt2.data <- dataset_intervals %>%
#   filter(drugclass == "SGLT2") %>%
#   mutate(prop_score = 1/(prop_score))
# 
# # weights for GLP1 Z = 0
# glp1.data <- dataset_intervals %>%
#   filter(drugclass == "GLP1") %>%
#   mutate(prop_score = 1/(1-prop_score))
# 
# new.data <- rbind(sglt2.data, glp1.data)
# 
# lm(weight.change~drugclass, new.data, weights = prop_score)
# 
# 
# 
# #####
# 
# dataset_new <- data_dev %>%
#   mutate(weight.change = postweight6m - preweight) %>%
#   cbind(hba1c_diff = effects_summary_dev$mean)
# 
# 
# data <- group_values(data = dataset_new, 
#                                   variable = "hba1c_diff", 
#                                   breaks = breaks) %>%
#   select(patid, pateddrug, drugclass, postweight6m, weight.change, hba1c_diff, intervals) %>%
#   drop_na() %>%
#   select(-postweight6m) %>%
#   rename("postweight6m" = "weight.change") %>%
#   mutate(hba1c_diff.q = as.numeric(intervals))
# 
# 
# # calculate propensity score
# prop_model <- calc_ATE_prop_score(data)
# 
# # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
# prop_score <- 1 - prop_model$p_hat_train
# 
# # split predicted treatment effects into deciles
# predicted_treatment_effect <- data %>%
#   plyr::ddply("hba1c_diff.q", dplyr::summarise,
#               N = length(hba1c_diff),
#               hba1c_diff.pred = mean(hba1c_diff))
# 
# # maximum number of deciles being tested
# quantiles <- max(data$hba1c_diff.q)
# 
# # create lists with results
# mnumber = c(1:quantiles)
# models  <- as.list(1:quantiles)
# hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector();
# 
# # join dataset and propensity score
# data.new <- data %>%
#   cbind(prop_score)
# 
# # weights for SGLT2 Z = 1
# sglt2.data <- data.new %>%
#   filter(drugclass == "SGLT2") %>%
#   mutate(prop_score = 1/(prop_score))
# 
# # weights for GLP1 Z = 0
# glp1.data <- data.new %>%
#   filter(drugclass == "GLP1") %>%
#   mutate(prop_score = 1/(1-prop_score))
# 
# data.new <- rbind(sglt2.data, glp1.data)
# 
# 
# # iterate through deciles
# for (i in mnumber) {
#   # fit linear regression for decile
#   models[[i]] <- lm(as.formula(postweight6m ~ factor(drugclass)),data=data.new,subset=hba1c_diff.q==i, weights = prop_score)
#   
# }
# 
#







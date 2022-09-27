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
  
  interim_dataset <- data_dev %>%
    mutate(posthba1c_final = postweight6m) %>%
    select(-postweight6m)
  
  cred_pred_dev <- calc_resid(interim_dataset, posteriors_dev)
  
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





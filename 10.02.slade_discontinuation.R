####################
## Description:
##  - In this file we find a way to know weight change expected:
##      1. fit a BART model for predicting discontinuation in a therapy.
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
dir.create(paste0(output_path, "/Discontinuation"))

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
              select(patid, pateddrug, score.excl.mi, stopdrug6m_3mFU)) %>%
  mutate(stopdrug6m_3mFU = factor(stopdrug6m_3mFU)) %>%
  drop_na(stopdrug6m_3mFU) # drop 699


## Fit model using all the available in treatment selection BART model variables to estimate weight change
if (class(try(
  
  bart_model_discontinuation <- readRDS(paste0(output_path, "/Weight_reduction/bart_model_discontinuation.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  bart_model_discontinuation <- bartMachine::bartMachine(X = dataset.dev %>%
                                                  select(drugclass,
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
                                                         prebmi,
                                                         prebil,
                                                         preplatelets,
                                                         t2dmduration,
                                                         prealb,
                                                         presys,
                                                         preast),
                                                y = dataset.dev[,"stopdrug6m_3mFU"],
                                                use_missing_data = TRUE,
                                                impute_missingness_with_rf_impute = FALSE,
                                                impute_missingness_with_x_j_bar_for_lm = TRUE,
                                                num_trees = 200,
                                                num_burn_in = 3000,
                                                num_iterations_after_burn_in = 1000,
                                                serialize = TRUE)
  
  saveRDS(bart_model_discontinuation, paste0(output_path, "/Weight_reduction/bart_model_discontinuation.rds"))
  
}


# bart_model_discontinuation chose everyone on not discontinuing


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
              select(patid, pateddrug, stopdrug6m_3mFU))


data_val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X))) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, stopdrug6m_3mFU))


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
breaks = c(-8, -5, -3, 0, 3, 5, 8)

######
####

dataset_breakdown_dev <- data_dev %>%
  select(stopdrug6m_3mFU, drugclass) %>%
  cbind(hba1c_diff = effects_summary_dev$mean)


dataset_intervals_dev <- group_values(data = dataset_breakdown_dev, 
                                      variable = "hba1c_diff", 
                                      breaks = breaks) %>%
  select(-hba1c_diff) %>%
  group_by(intervals, drugclass) %>%
  mutate(mean.values = mean(stopdrug6m_3mFU, na.rm = TRUE)) %>%
  select(-stopdrug6m_3mFU) %>%
  unique()


plot_discontinuation_dev <- dataset_intervals_dev %>%
  ggplot() +
  geom_point(aes(x = intervals,
                 y = mean.values,
                 colour = drugclass),
             position = position_dodge(width = 0.5)) +
  ylab("% discontinuing treatment within 6 months") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")
  

####

dataset_breakdown_val <- data_val %>%
  select(stopdrug6m_3mFU, drugclass) %>%
  cbind(hba1c_diff = effects_summary_val$mean)


dataset_intervals_val <- group_values(data = dataset_breakdown_val, 
                                      variable = "hba1c_diff", 
                                      breaks = breaks) %>%
  select(-hba1c_diff) %>%
  group_by(intervals, drugclass) %>%
  mutate(mean.values = mean(stopdrug6m_3mFU, na.rm = TRUE)) %>%
  select(-stopdrug6m_3mFU) %>%
  unique()


plot_discontinuation_val <- dataset_intervals_val %>%
  ggplot() +
  geom_point(aes(x = intervals,
                 y = mean.values,
                 colour = drugclass),
             position = position_dodge(width = 0.5)) +
  ylab("% discontinuing treatment within 6 months") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom")



plot_discontinuation <- patchwork::wrap_plots(
  # Development
  plot_discontinuation_dev,
  # Validation
  plot_discontinuation_val
) + patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')



pdf(file = "Plots/10.2.discontinuation_plots.pdf")
plot_discontinuation
dev.off()







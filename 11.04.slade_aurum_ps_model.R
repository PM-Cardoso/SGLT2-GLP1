####################
## Description:
##  - In this file we:
##    - Fit a BART PS model to all variables.
##    - Perform variable selection on PS model and refit model.
##    - Refit BART PS model
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

output_path <- "Samples/SGLT2-GLP1/Aurum"

## make directory for outputs
dir.create(output_path)

## male directory for outputs
dir.create(paste0(output_path, "/ps_model"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")



ps.model.train <- set_up_data_sglt2_glp1(dataset.type = "ps.model.train")

###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################


########
### Fit a propensity score model to all variables in dataset
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(

  bart_ps_model <- readRDS(paste0(output_path, "/ps_model/bart_ps_model.rds"))

  , silent = TRUE)) == "try-error") {

  set.seed(123)
  bart_ps_model <- bartMachine::bartMachine(X = ps.model.train %>%
                                              select(-patid,
                                                     -pated,
                                                     -drugclass),
                                            y = ps.model.train[,"drugclass"] %>%
                                              mutate(drugclass = factor(drugclass)) %>%
                                              unlist(),
                                            use_missing_data = TRUE,
                                            num_trees = 50,
                                            num_burn_in = 2000,
                                            num_iterations_after_burn_in = 1000,
                                            serialize = TRUE,
                                            seed = 123)

  saveRDS(bart_ps_model, paste0(output_path, "/ps_model/bart_ps_model.rds"))

}


# Check the optimal number of trees
## Does not work for classification


# Cross validation of model
if (class(try(
  
  bart_ps_model_cv <- readRDS(paste0(output_path, "/ps_model/bart_ps_model_cv.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  set.seed(123)
  bart_ps_model_cv <- bartMachine::bartMachineCV(X = ps.model.train %>%
                                                   select(-patid,
                                                          -pated,
                                                          -drugclass),
                                                 y = ps.model.train[,"drugclass"] %>%
                                                   mutate(drugclass = factor(drugclass)) %>%
                                                   unlist(),
                                                 num_tree_cvs = c(50, 100, 150, 200),
                                                 k_cvs = c(2, 3),
                                                 use_missing_data = TRUE,
                                                 num_burn_in = 2000,
                                                 num_iterations_after_burn_in = 1000,
                                                 serialize = TRUE,
                                                 seed = 123)
  
  saveRDS(bart_ps_model_cv, paste0(output_path, "/ps_model/bart_ps_model_cv.rds"))
  
}


########
### Variable selection of propensity score model
########

if (class(try(

  vs_bart_ps_model <- readRDS(paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))

  , silent = TRUE)) == "try-error") {

  pdf(file = "Plots/11.04.prop_model_vs.pdf", width = 18, height = 11)
  # error with cv
  set.seed(123)
  vs_bart_ps_model <- var_selection_by_permute(bart_ps_model)
  dev.off()

  ## Variables selected


  saveRDS(vs_bart_ps_model, paste0(output_path, "/ps_model/vs_bart_ps_model.rds"))

}


########
### Refit PS model with selected vars
########

## Fit initial model using all the available variables to estimate HbA1c outcome
if (class(try(

  bart_ps_model_final <- readRDS(paste0(output_path, "/ps_model/bart_ps_model_final.rds"))

  , silent = TRUE)) == "try-error") {

  set.seed(123)
  bart_ps_model_final <- bartMachine::bartMachineCV(X = ps.model.train %>%
                                                    select(
                                                      vs_bart_ps_model$important_vars_local_names
                                                      ),
                                                  y = ps.model.train[,"drugclass"] %>%
                                                    mutate(drugclass = factor(drugclass)) %>%
                                                    unlist(),
                                                  num_tree_cvs = c(50, 100, 150, 200),
                                                  k_cvs = c(2, 3),
                                                  use_missing_data = TRUE,
                                                  num_burn_in = 2000,
                                                  num_iterations_after_burn_in = 1000,
                                                  serialize = TRUE,
                                                  seed = 123)

  saveRDS(bart_ps_model_final, paste0(output_path, "/ps_model/bart_ps_model_final.rds"))

}


# 
# predictions <- bart_ps_model_final$p_hat_train
# 
# df <- as.data.frame(cbind(predictions, ps.model.train[,"drugclass"] %>%
#                             mutate(drugclass = ifelse(drugclass == 2, 0, 1)))) %>%
#   set_names(c("predictions", "labels"))
# 
# library(ROCR)
# pred <- prediction(df$predictions, df$labels)
# perf <- performance(pred,"tpr","fpr")
# pdf()
# plot(perf,colorize=TRUE)
# dev.off()
# 
# library(pROC)
# pROC_obj <- roc(df$labels,df$predictions,
#                 smoothed = TRUE,
#                 # arguments for ci
#                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE)
# 
# 
# 
# values <- coords(pROC_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
#                                          "precision", "recall"), transpose = FALSE)
# 
# 
# as.data.frame(cbind(predictions, ps.model.train[,"drugclass"])) %>%
#   set_names(c("predictions", "labels")) %>%
#   mutate(predictions = ifelse(predictions > values$threshold, 1, 2)) %>% table()



###############################################################################
###############################################################################
###################### Propensity scores for patients #########################
###############################################################################
###############################################################################


# Prop scores for train dataset
patient_prop_scores <- ps.model.train %>%
  select(patid, pated) %>%
  cbind(prop.score = bart_ps_model_final$p_hat_train)


# Prop scores for test dataset

ps.model.test <- set_up_data_sglt2_glp1(dataset.type = "ps.model.test")

# calculate prop score
if (class(try(

  prop_score_testing_data <- readRDS(paste0(output_path, "/ps_model/prop_score_testing_data.rds"))

  , silent = TRUE)) == "try-error") {

  set.seed(123)
  prop_score_testing_data <- predict(bart_ps_model_final, ps.model.test %>%
                        select(
                          colnames(bart_ps_model_final$X)
                        ))

  saveRDS(prop_score_testing_data, paste0(output_path, "/ps_model/prop_score_testing_data.rds"))

}

predictions <- prop_score_testing_data

df <- as.data.frame(cbind(predictions, ps.model.test[,"drugclass"] %>%
                            mutate(drugclass = ifelse(drugclass == 2, 0, 1)))) %>%
  set_names(c("predictions", "labels"))

library(ROCR)
pred <- prediction(df$predictions, df$labels)
perf <- performance(pred,"tpr","fpr")
pdf()
plot(perf,colorize=TRUE)
dev.off()

library(pROC)
pROC_obj <- roc(df$labels,df$predictions,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)



values <- coords(pROC_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                         "precision", "recall"), transpose = FALSE)


as.data.frame(cbind(predictions, ps.model.test[,"drugclass"])) %>%
  set_names(c("predictions", "labels")) %>%
  mutate(predictions = ifelse(predictions > values$threshold, 1, 2)) %>% table()


patient_prop_scores <- patient_prop_scores %>%
  rbind(
    ps.model.test %>%
      select(patid, pated) %>%
      cbind(prop.score = prop_score_testing_data)
  ) %>%
  as.data.frame()

saveRDS(patient_prop_scores, paste0(output_path, "/ps_model/patient_prop_scores.rds"))



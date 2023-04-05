####################
## Description:    (this file uses BCF v.2.0.1)
##  - In this file we:
##    - Fit a similar version of the model but with reduced size to enable its use in a shiny app.
####################


## Load libraries
library(tidyverse)
require(bcf)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

dir.create(paste0(output_path, "/response_model_bcf"))

dir.create(paste0(output_path, "/response_model_bcf/shiny_app"))

dir.create(paste0(output_path, "/response_model_bcf/shiny_app/trees_no_prop"))

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("11.01.slade_aurum_functions.R")
source("11.02.slade_aurum_set_data.R")

## Load dataset
hba1c.train <- set_up_data_sglt2_glp1(dataset.type = "hba1c.train")

## Collect propensity score values
patient_prop_scores <- readRDS(paste0(output_path, "/ps_model/patient_prop_scores.rds"))

## Append values
hba1c.train <- hba1c.train %>%
  left_join(patient_prop_scores, by = c("patid", "pated"))

## Prepare dataset to be used in SparseBCF
hba1c.train.complete <- hba1c.train %>%
  # drop the variables with the most missingness (>40%)
  select(-preacr, -preast, -prehaematocrit, -prehaemoglobin, -pretriglyceride) %>%
  # only complete cases
  drop_na() %>%
  as.data.frame()


variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))


hba1c.train.complete.vs <- hba1c.train.complete %>%
  select(patid, pated, posthba1cfinal, drugclass, unique(c(variables_mu, variables_tau)), prop.score)


# Fit BCF model without propensity score included

if (class(try(
  
  bcf_model <- readRDS(paste0(output_path, "/response_model_bcf/shiny_app/bcf_model.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  
  bcf_model = bcf::bcf(y = hba1c.train.complete.vs$posthba1cfinal,
                       z = hba1c.train.complete.vs %>%
                         select(drugclass) %>%
                         mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                         unlist(),
                       x_control = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_mu)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       x_moderate = hba1c.train.complete.vs %>%
                         select(
                           all_of(variables_tau)
                         ) %>%
                         mutate_all(funs(as.numeric(.))) %>%
                         as.matrix(),
                       pihat = 1-hba1c.train.complete.vs$prop.score,
                       nburn = 200000,
                       nsim = 150,
                       nthin = 149,
                       n_chains = 2,
                       # n_threads was max((RcppParallel::defaultNumThreads() - 2)/n_cores, 1) (this uses all of the server)
                       n_threads = 4,
                       update_interval = 500,
                       ntree_control = 200,
                       sd_control = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_control = 0.95,
                       power_control = 2,
                       ntree_moderate = 200,
                       sd_moderate = 2 * sd(hba1c.train.complete.vs$posthba1cfinal),
                       base_moderate = 0.95,
                       power_moderate = 2,
                       use_muscale = FALSE,
                       use_tauscale = FALSE,
                       include_pi = "none",
                       save_tree_directory = paste0(output_path, "/response_model_bcf/shiny_app/trees_no_prop"))
  
  saveRDS(bcf_model, paste0(output_path, "/response_model_bcf/shiny_app/bcf_model.rds"))
  
}

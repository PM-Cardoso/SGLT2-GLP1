####################
## Description:
##  - This file includes a table descriptive analyses of development
##      and validation datasets.
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])

library(tidyverse)
library(psych)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path,"/Data_Analysis"))

## make directory for outputs
dir.create("Plots")



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

# name: final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))


###############################################################################
###############################################################################
################################# Analysis ####################################
###############################################################################
###############################################################################

dataset.dev <- final.dev %>%
  select(posthba1c_final,
         drugclass,
         egfr_ckdepi,
         hba1cmonth,
         prealt,
         prehba1cmmol,
         score,
         Category,
         drugline,
         ncurrtx,
         yrdrugstart,
         agetx,
         malesex,
         prehdl,
         prebmi,
         prebil,
         preplatelets,
         t2dmduration,
         prealb,
         presys,
         preast
  )


## fit table
describeBy(dataset.dev, group = dataset.dev$drugclass, fast = TRUE)

table(dataset.dev$drugline, dataset.dev$drugclass)

table(dataset.dev$ncurrtx, dataset.dev$drugclass)

table(dataset.dev$Category, dataset.dev$drugclass)

table(dataset.dev$malesex, dataset.dev$drugclass)

####


dataset.val <- final.val %>%
  select(posthba1c_final,
         drugclass,
         egfr_ckdepi,
         hba1cmonth,
         prealt,
         prehba1cmmol,
         score,
         Category,
         drugline,
         ncurrtx,
         yrdrugstart,
         agetx,
         malesex,
         prehdl,
         prebmi,
         prebil,
         preplatelets,
         t2dmduration,
         prealb,
         presys,
         preast
  )

## fit table
describeBy(dataset.val, group = dataset.val$drugclass, fast = TRUE)

table(dataset.val$drugline, dataset.val$drugclass)

table(dataset.val$ncurrtx, dataset.val$drugclass)

table(dataset.val$Category, dataset.val$drugclass)

table(dataset.val$malesex, dataset.val$drugclass)


###########################################
###########################################
###########################################

dataset.dev.post <- dataset.dev %>%
  filter(yrdrugstart > 2012)

## fit table
describeBy(dataset.dev.post, group = dataset.dev.post$drugclass, fast = TRUE)

table(dataset.dev.post$drugline, dataset.dev.post$drugclass)

table(dataset.dev.post$ncurrtx, dataset.dev.post$drugclass)

table(dataset.dev.post$Category, dataset.dev.post$drugclass)

table(dataset.dev.post$malesex, dataset.dev.post$drugclass)


####


dataset.val.post <- dataset.val %>%
  filter(yrdrugstart > 2012)

## fit table
describeBy(dataset.val.post, group = dataset.val.post$drugclass, fast = TRUE)

table(dataset.val.post$drugline, dataset.val.post$drugclass)

table(dataset.val.post$ncurrtx, dataset.val.post$drugclass)

table(dataset.val.post$Category, dataset.val.post$drugclass)

table(dataset.val.post$malesex, dataset.val.post$drugclass)























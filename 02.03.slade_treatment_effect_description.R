####################
## Description:
##  - This file includes a table/plot descrition of varying treatment effect
##      for all variables
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

# name: final.all.extra.vars
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_allcohort.Rda"))

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

# name: final.val
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))




###############################################################################
###############################################################################
################################# Analysis ####################################
###############################################################################
###############################################################################

# bart model

bart_model_final <- readRDS(paste0(output_path, "/Final_model/cvd_new/bart_model_final.rds"))


# Dev dataset
dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

# Val dataset
dataset.val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

# treatment effects dev

effects_summary_dev <- readRDS(paste0(output_path, "/Final_model/cvd_new/Assessment/effects_summary_dev.rds"))

# treatment effects val

effects_summary_val <- readRDS(paste0(output_path, "/Final_model/cvd_new/Assessment/effects_summary_val.rds"))

#################
############ Extreme deciles characteristics
#################

full.dataset <- rbind( # join dev + val datasets
  dataset.dev %>%
    # join hba1c_diff - treatment effect
    cbind(hba1c_diff = effects_summary_dev$mean),
  dataset.val %>%
    cbind(hba1c_diff = effects_summary_val$mean)
) %>%
  # select extremes by above 7.5 mmol/mol
  filter(hba1c_diff > 7.5 | hba1c_diff < -7.5) %>%
  mutate(favor = ifelse(hba1c_diff < 0, "SGLT2", "GLP1")) %>%
  # select extremes by deciles
  # mutate(
  #   hba1c_diff.q = ntile(hba1c_diff, 10)
  # ) %>%
  # filter(hba1c_diff.q %in% c(1,10)) %>%
  select(-patid, -pateddrug, -posthba1c_final, -hba1c_diff)


describeBy(full.dataset, group = as.factor(full.dataset$favor))

# table of drug taken vs drug preferred
table(full.dataset$drugclass, full.dataset$favor)

# table of sex vs drug preferred
table(full.dataset$malesex, full.dataset$favor)


#################
############ Plot variable characteristics and colour of therapy
#################

library(GGally)

dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

dat1 <- effects_summary_dev %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


combination <- cbind(dataset.dev %>% select(-patid, -pateddrug, -posthba1c_final, -hba1cmonth, -drugline, -yrdrugstart, -ncurrtx), Benefit = dat1$above)


my_dens_lower <- function(data, mapping, ...) {
  ggplot(mapping=mapping) +
    geom_density2d(data = filter(data, Benefit == "Favours GLP1"), size = 1.2, alpha = 0.7, colour = "red") +
    geom_density2d(data = filter(data, Benefit == "Favours SGLT2"), size = 1.2, alpha = 0.7, colour = '#f1a340')
}

my_dens_diagonal <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_density(aes(fill = Benefit), alpha = 0.7) +
    scale_fill_manual(values = c("red", "#f1a340")) +
    scale_colour_manual(values = c("red", "#f1a340"))
}


plot_characteristics <- ggpairs(combination, columns = 1:(ncol(combination)-1), 
        aes(color = Benefit),
        showStrips = TRUE,
        lower = list(continuous = my_dens_lower, discrete = wrap(ggally_facetbar, position = "dodge", alpha = 0.7), combo = wrap(ggally_facetdensity,alpha=0.7)),
        diag = list(continuous = my_dens_diagonal, discrete = wrap(ggally_barDiag, position = "dodge", alpha = 0.7)),
        upper = NULL,
        legend = 1,
        title = "Characteristics") +
  theme(legend.position = 'bottom',
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank()) +
  scale_fill_manual(values = c("red", "#f1a340")) +
  scale_colour_manual(values = c("red", "#f1a340"))







#### PDF with all the plots

pdf(file = "Plots/2.3.effect_characteristics.pdf", width = 20, height = 20)
plot_characteristics
dev.off()





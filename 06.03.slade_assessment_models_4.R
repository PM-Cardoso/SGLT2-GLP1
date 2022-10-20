####################
## Description:
##  - In this file we make a collected of plots comparing assessment
##      measurements of several models from 4.1 to 4.4.
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")


library(tidyverse)
# library(bartMachine)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Plots #########################################
###############################################################################
###############################################################################


assessment <- rbind(
  cbind(readRDS(paste0(output_path, "/Final_model/model_1/Assessment/assessment.rds")), Model = "Model 1"),
  cbind(readRDS(paste0(output_path, "/Final_model/model_2/Assessment/assessment.rds")), Model = "Model 2"),
  cbind(readRDS(paste0(output_path, "/Final_model/model_3/Assessment/assessment.rds")), Model = "Model 3"),
  cbind(readRDS(paste0(output_path, "/Final_model/model_4/Assessment/assessment.rds")), Model = "Model 4")
  ) %>%
  as.data.frame() %>%
  mutate(`5%` = as.numeric(`5%`),
         `50%` = as.numeric(`50%`),
         `95%` = as.numeric(`95%`),
         Model = factor(Model, levels = c("Model 4", "Model 3", "Model 2", "Model 1")))


plot_assessment <- assessment %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(y = Model, xmin = `5%`, xmax = `95%`, colour = Model), width = 0.2) +
  geom_point(aes(x = `50%`, y = Model, shape = Dataset), size = 2, colour = "black") +
  facet_wrap(~statistic, ncol = 1, scales = "free") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  guides(colour = "none")




#### PDF with all the plots


pdf(file = "Plots/6.3.assessment.pdf")
plot_assessment
dev.off()


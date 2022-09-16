####################
## Description:
##  - In this file we investigate variable importance
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")

library(tidyverse)
library(bartMachine)
library(forcats)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path, "/Final_model"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/With_grf_no_prop"))


## make directory for outputs
dir.create(paste0(output_path, "/Final_model/With_grf_no_prop/Assessment"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################## Read Data in ###################################
###############################################################################
###############################################################################

# bart model
bart_model_final <- readRDS(paste0(output_path, "/Final_model/With_grf_no_prop/bart_model_final.rds"))



# variable importance

if (class(try(
  
  variable_importance <- readRDS(paste0(output_path, "/Final_model/With_grf_no_prop/Assessment/variable_importance.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  pdf(file = paste0(output_path, "/Final_model/With_grf_no_prop/Assessment/variable_importance.pdf"))
  
  variable_importance <- bartMachine::investigate_var_importance(bart_model_final)
  
  dev.off()
  
  saveRDS(variable_importance, paste0(output_path, "/Final_model/With_grf_no_prop/Assessment/variable_importance.rds"))
  
}



plot_var_importance <- variable_importance$avg_var_props %>%
  t() %>%
  as.data.frame() %>%
  mutate(drugclass = drugclass_SGLT2 + drugclass_GLP1,
         drugline = drugline_2 + drugline_3 + drugline_4 + drugline_5,
         ncurrtx = ncurrtx_0 + ncurrtx_1 + ncurrtx_2 + ncurrtx_3,
         malesex = malesex_0 + malesex_1,
         Category = `Category_Active smoker` + `Category_Ex-smoker` + `Category_Non-smoker`) %>%
  select(-drugclass_SGLT2, -drugclass_GLP1, -drugline_2, -drugline_3, -drugline_4, -drugline_5, -ncurrtx_0, -ncurrtx_1, -ncurrtx_2, -ncurrtx_3, -malesex_0, -malesex_1, -`Category_Active smoker`, -`Category_Ex-smoker`, -`Category_Non-smoker`) %>%
  rename("Therapy" = "drugclass",
         "HbA1c" = "prehba1cmmol",
         "Year Drug Start" = "yrdrugstart",
         "eGFR" = "egfr_ckdepi",
         "Outcome month" = "hba1cmonth",
         "Systolic" = "presys",
         "ALT" = "prealt",
         "Bilirubin" = "prebil",
         "AST" = "preast",
         "HDL" = "prehdl",
         "Age" = "agetx",
         "CVD score" = "score",
         "Time to prescription" = "t2dmduration",
         "BMI" = "prebmi",
         "Albuminuria" = "prealb",
         "Platelets" = "preplatelets",
         "Number of past therapies" = "drugline",
         "Number of Current therapies" = "ncurrtx",
         "Sex" = "malesex",
         "Smoking status" = "Category") %>%
  gather() %>%
  mutate(value = value * 100) %>%
  ggplot(aes(x = fct_rev(fct_reorder(key, value)), y = value)) +
  theme_classic() +
  geom_col() +
  ylab("Inclusion proportions (%)") +
  scale_y_continuous(limits = c(0, 25), expand = (c(0,0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate(x = 0, xend=0, y=0, yend=25, colour="black", lwd=0.75, geom="segment")


pdf(file = paste0(output_path, "/Final_model/With_grf_no_prop/Assessment/variable_importance_2.pdf"))
plot_var_importance
dev.off()











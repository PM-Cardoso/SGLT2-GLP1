####################
## Description:
##  - In this file we compare the predictions of SGLT2/GLP1 from the BCF model
##    to the predictions form the SGLT2/DPP4 linear model.
####################


library(rms)
library(tidyverse)

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

## make directory for outputs
# dir.create(paste0(output_path, "/linear_model_comparison"))

## make directory for outputs
dir.create("Plots")


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

source("11.02.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################

# read in predictions from BCF model
patient_predicted_outcomes <- readRDS("Samples/SGLT2-GLP1/Aurum/response_model_bcf/patient_predicted_outcomes.rds")

# load in linear regression model
load("m1_hba1cmodel_SGLT2_DPP4.Rdata")


# read in full cohort with modifications for linear model
full.cohort.updated <- set_up_data_sglt2_glp1(dataset.type = "full.cohort") %>%
  select(patid, pated, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth) %>%
  rename("prealtlog"="prealt",
         "prehba1cmmol"="prehba1c",
         "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog),
         ncurrtx = factor(ncurrtx, levels = c("1","2","3","4","5+"), labels = c("0","1","2","3","3")),
         drugline = factor(drugline, levels = c("2","3","4","5+"), labels = c("2","3","4","5")),
         hba1cmonth = ifelse(is.na(hba1cmonth), 12, hba1cmonth)) %>%
  drop_na() %>%
  filter(prehba1cmmol < 120) %>%
  filter(egfr_ckdepi > 45)


# predict for SGLT2
dataset.sglt2 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("SGLT2", levels = c("SGLT2", "DPP4")))
predictions.sglt2 <- predict(m1, dataset.sglt2)

# predict for DPP4
dataset.dpp4 <- full.cohort.updated %>%
  select(-patid, -pated) %>%
  mutate(drugclass = factor("DPP4", levels = c("SGLT2", "DPP4")))
predictions.dpp4 <- predict(m1, dataset.dpp4)

# treatment effects
effects <- predictions.sglt2-predictions.dpp4



## Compare SGLT2 BCF vs Linear regression
interim.dataset <- full.cohort.updated %>%
  select(patid, pated) %>%
  cbind(pred.SGLT2.lm = predictions.sglt2) %>%
  left_join(patient_predicted_outcomes %>%
              select(patid, pated, pred.SGLT2) %>%
              rename("pred.SGLT2.bcf" = "pred.SGLT2"), by = c("patid", "pated")) %>%
  select(-patid, -pated)

plot_comparison <- interim.dataset %>%
  ggplot(aes(y = pred.SGLT2.bcf, x = pred.SGLT2.lm)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red", linetype = "dashed") +
  stat_smooth() +
  ylab("SGLT2 predictions with BCF") +
  xlab("SGLT2 predictions with Linear Regression") +
  ylim(min(interim.dataset), max(interim.dataset)) +
  xlim(min(interim.dataset), max(interim.dataset)) +
  ggtitle(paste0("Comparison of SGLT2 predictions using Linear Regression and BCF (n=", interim.dataset %>% nrow(), ")")) +
  theme_bw()

  
pdf(width = 7, height = 7, "Plots/11.09.plot_1.pdf")
plot_comparison
dev.off()
  
## What drug is best:
#-----------------
# Using SGLT2 BCF
interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "DPP4i", NA))),
         best_drug = factor(best_drug))

plot_bar <- interim.dataset %>%
  select(best_drug) %>%
  table() %>%
  as.data.frame() %>%
  rename("best_drug" = ".") %>%
  ggplot(aes(x = best_drug, y = Freq, fill = best_drug)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = Freq), vjust = -0.5) +
  xlab("Optimal predicted therapy") +
  ylab("Number of patients") +
  ggtitle(paste0("Predicted optimal therapy (n=", interim.dataset%>%nrow(), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_bw() +
  theme(legend.position = "none")

pdf(width = 7, height = 7, "Plots/11.09.plot_2.pdf")
plot_bar
dev.off()

# # Using SGLT2 Linear
# plot_bar <- patient_predicted_outcomes %>%
#   select(-pred.SGLT2) %>%
#   left_join(full.cohort.updated %>%
#               select(patid, pated) %>%
#               cbind(pred.DPP4 = predictions.dpp4,
#                     pred.SGLT2 = predictions.sglt2), by = c("patid", "pated")) %>%
#   drop_na() %>%
#   mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "SGLT2i",
#                             ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "GLP1-RA",
#                                    ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "DPP4i", NA))),
#          best_drug = factor(best_drug)) %>%
#   select(best_drug) %>%
#   count(best_drug) %>%
#   ggplot(aes(x = best_drug, y = n, fill = best_drug)) +
#   geom_bar(stat="identity") +
#   xlab("Optimal predicted therapy") +
#   ylab("Number of patients") +
#   scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
#   theme_bw() +
#   theme(legend.position = "none")
  

#---------------
interim.dataset <- patient_predicted_outcomes %>%
  left_join(full.cohort.updated %>%
              select(patid, pated) %>%
              cbind(pred.DPP4 = predictions.dpp4), by = c("patid", "pated")) %>%
  drop_na() %>%
  mutate(best_drug = ifelse(pred.SGLT2 < pred.GLP1 & pred.SGLT2 < pred.DPP4, "Favours SGLT2i",
                            ifelse(pred.GLP1 < pred.SGLT2 & pred.GLP1 < pred.DPP4, "Favours GLP1-RA",
                                   ifelse(pred.DPP4 < pred.SGLT2 & pred.DPP4 < pred.GLP1, "Favours DPP4i", NA))),
         best_drug = factor(best_drug),
         effect = ifelse(best_drug == "Favours SGLT2i" & pred.GLP1 < pred.DPP4, pred.SGLT2 - pred.GLP1,
                         ifelse(best_drug == "Favours SGLT2i" & pred.DPP4 < pred.GLP1, pred.SGLT2 - pred.DPP4,
                                ifelse(best_drug == "Favours GLP1-RA" & pred.SGLT2 < pred.DPP4, pred.GLP1 - pred.SGLT2,
                                       ifelse(best_drug == "Favours GLP1-RA" & pred.DPP4 < pred.SGLT2, pred.GLP1 - pred.DPP4,
                                              ifelse(best_drug == "Favours DPP4i" & pred.SGLT2 < pred.GLP1, pred.DPP4 - pred.SGLT2,
                                                     ifelse(best_drug == "Favours DPP4i" & pred.GLP1 < pred.SGLT2, pred.DPP4 - pred.GLP1, NA))))))) 

plot_histogram <- interim.dataset %>%
  select(best_drug, effect) %>%
  ggplot(aes(x = effect, fill = best_drug)) +
  geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(-15, 0, by = 1)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  xlab("Predicted HbA1c benefit (mmol/mol)") +
  ggtitle(paste0("Predicted HbA1c benefit against the next best therapy (n=", interim.dataset %>% nrow(), ")")) +
  scale_fill_manual(values = c("red", "dodgerblue2", "#f1a340")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~best_drug, nrow = 1)
  
  
pdf(width = 9, height = 4, "Plots/11.09.plot_3.pdf")
plot_histogram
dev.off()



## merge pdfs

qpdf::pdf_combine(input = c("Plots/11.09.plot_1.pdf",
                            "Plots/11.09.plot_2.pdf",
                            "Plots/11.09.plot_3.pdf"),
                  output = "Plots/11.09.comparison_SGLT2_GLP1_DPP4.pdf")

file.remove(c("Plots/11.09.plot_1.pdf", "Plots/11.09.plot_2.pdf", "Plots/11.09.plot_3.pdf"))



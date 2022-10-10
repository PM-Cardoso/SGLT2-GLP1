####################
## Description:
##  - In this file we plot all the figures used in the SGLT2 vs GLP1 paper
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])

## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")

library(tidyverse)
library(bartMachine)

## make directory for outputs
dir.create("Samples")

## make directory for outputs
dir.create("Samples/SGLT2-GLP1")


#####

# Plot 1 

# 1. Histogram of treatment benefit, overall, female, male


# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/effects_summary_dev.rds")

## hist_plot function
dat1 <- effects_summary_dev %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))



# Plot overall treatment effect benefit
plot_effects_1 <- ggplot(data=dat1, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Overall Population",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# Stratify by sex
effects_summary_dev_male <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 1)

effects_summary_dev_female <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 0)

dat2 <- effects_summary_dev_male %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


# Plot treatment effect of male sex
plot_effects_2 <- ggplot(data=dat2, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Stratified: Male sex",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

dat3 <- effects_summary_dev_female %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


# Plot treatment effect of female sex
plot_effects_3 <- ggplot(data=dat3, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Stratified: Female sex",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

plot_1 <- patchwork::wrap_plots(list(plot_effects_1, plot_effects_2, plot_effects_3), ncol = 3) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pdf("Plot1.pdf", width = 12, height = 4)
plot_1
dev.off()



# 2. Differential treatment effect, only important vars (eGFR, age, sex)

# 3. A validation of outcome, plot: y-defined outcomes, x-ATE
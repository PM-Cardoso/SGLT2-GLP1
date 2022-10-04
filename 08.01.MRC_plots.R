####################
## Description:
##  - In this file we plot all the figures used in the presentation
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

## make directory for outputs
dir.create("Samples/SGLT2-GLP1/Plots")

## make directory for outputs
dir.create("Samples/SGLT2-GLP1/Plots/MRC_plots")

#####

# Plot 1 

## Averate treatment effect for population + sex split

# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/effects_summary_dev.rds")

## hist_plot function
dat1 <- effects_summary_dev %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


plot_effects_1 <- ggplot(data=dat1, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Overall Population",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# breaks by quantile
# breaks <- c(
#   -10,
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 1) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 2) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 2) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 3) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 3) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 4) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 4) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 5) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 5) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 6) %>% select(hba1c_diff) %>% min()),
#   0,
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 6) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 7) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 7) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 8) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 8) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 9) %>% select(hba1c_diff) %>% min()),
#   mean(predicted_observed_dev %>% filter(hba1c_diff.q == 9) %>% select(hba1c_diff) %>% max(), predicted_observed_dev %>% filter(hba1c_diff.q == 10) %>% select(hba1c_diff) %>% min()),
#   10
# )
# plot_effects_1 <- ggplot(data=dat1, aes(x=mean,fill=above)) +
#   geom_histogram(position="identity", alpha=0.5,color="black", breaks = breaks) +
#   geom_vline(aes(xintercept=0), linetype="dashed")+
#   labs(title="Overall Population",x="HbA1c difference (mmol/mol)", y = "Number of people") +
#   scale_fill_manual(values=c("red","#f1a340"))+
#   theme_classic() +
#   theme(legend.position = c(0.80, 0.97)) +
#   theme(legend.title = element_blank())


effects_summary_dev_male <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 1)

effects_summary_dev_female <- effects_summary_dev %>%
  cbind(malesex = final.dev$malesex) %>%
  filter(malesex == 0)

dat2 <- effects_summary_dev_male %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


plot_effects_2 <- ggplot(data=dat2, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Male",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# breaks by quantile
# plot_effects_2 <- ggplot(data=dat2, aes(x=mean,fill=above)) +
#   geom_histogram(position="identity", alpha=0.5,color="black",breaks=breaks) +
#   geom_vline(aes(xintercept=0), linetype="dashed")+
#   labs(title="Male",x="HbA1c difference (mmol/mol)", y = "Number of people") +
#   scale_fill_manual(values=c("red","#f1a340"))+
#   theme_classic() +
#   theme(legend.position = c(0.80, 0.97)) +
#   theme(legend.title = element_blank())


dat3 <- effects_summary_dev_female %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


plot_effects_3 <- ggplot(data=dat3, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Female",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# breaks by quantile
# plot_effects_3 <- ggplot(data=dat3, aes(x=mean,fill=above)) +
#   geom_histogram(position="identity", alpha=0.5,color="black",breaks=breaks) +
#   geom_vline(aes(xintercept=0), linetype="dashed")+
#   labs(title="Female",x="HbA1c difference (mmol/mol)", y = "Number of people") +
#   scale_fill_manual(values=c("red","#f1a340"))+
#   theme_classic() +
#   theme(legend.position = c(0.80, 0.97)) +
#   theme(legend.title = element_blank())


plot_1 <- cowplot::plot_grid(
  plot_effects_1,
  plot_effects_2,
  plot_effects_3,
  ncol = 3, nrow = 1
)

pdf("Plot1.pdf", width = 12, height = 4)
plot_1
dev.off()


#####

# Plot 2

## Differential treatment effect for variables

effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Comparison/Model_4_4/effects_summary_dev.rds")

# load all data for range of variable values; name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# HbA1c

# plot treatment effect + histogram marginal
effects <- effects_summary_dev[["prehba1cmmol"]]
variable <- "prehba1cmmol"
xtitle <- "HbA1c"

plot_hist <- ggplot() +
  theme_classic() +
  geom_histogram(aes(x = final.all.extra.vars[, variable])) +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "white"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())


plot_diff <- effects %>%
  group_by(ntile, ntile.value) %>%
  mutate(`5%`= quantile(mean, probs = c(0.05)),
         `50%` = quantile(mean, probs = c(0.50)),
         `95%` = quantile(mean, probs = c(0.95)),
         mean = mean(mean)) %>%
  ungroup() %>%
  unique() %>%
  ggplot() +
  geom_line(aes(x = ntile.value, y = mean), col = "red") +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`, x = ntile.value), alpha = 0.1) +
  xlab(xtitle) + ylab("Treatment Effect")


plot.hba1c.diff.marg <- cowplot::plot_grid(
  
  plot_diff
  
  ,
  
  cowplot::plot_grid(
    
    ggplot() +
      theme_void()
    
    ,
    
    plot_hist
    
    , ncol = 2, nrow = 1, rel_widths = c(0.06, 1)
    
  )
  
  , ncol = 1, nrow = 2, rel_heights = c(0.85, 0.15)
  
)


# eGFR

# plot treatment effect + histogram marginal
effects <- effects_summary_dev[["egfr_ckdepi"]]
variable <- "egfr_ckdepi"
xtitle <- "eGFR"

plot_hist <- ggplot() +
  theme_classic() +
  geom_histogram(aes(x = final.all.extra.vars[, variable])) +
  xlim(45, 155) +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "white"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())
  
  
plot_diff <- effects %>%
  group_by(ntile, ntile.value) %>%
  mutate(`5%`= quantile(mean, probs = c(0.05)),
         `50%` = quantile(mean, probs = c(0.50)),
         `95%` = quantile(mean, probs = c(0.95)),
         mean = mean(mean)) %>%
  ungroup() %>%
  unique() %>%
  ggplot() +
  geom_line(aes(x = ntile.value, y = mean), col = "red") +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`, x = ntile.value), alpha = 0.1) +
  xlab(xtitle) + ylab("Treatment Effect") +
  xlim(45, 155)


plot.egfr.diff.marg <- cowplot::plot_grid(
  
  plot_diff
  
  ,
  
  cowplot::plot_grid(
    
    ggplot() +
      theme_void()
    
    ,
    
    plot_hist
    
    , ncol = 2, nrow = 1, rel_widths = c(0.06, 1)
    
  )
  
  , ncol = 1, nrow = 2, rel_heights = c(0.85, 0.15)
  
)

# Sex

# plot treatment effect + histogram marginal
effects <- effects_summary_dev[["malesex"]]
variable <- "malesex"
xtitle <- "Sex"

plot_hist <- ggplot() +
  theme_classic() +
  geom_bar(aes(x = final.all.extra.vars[, variable])) +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(color = "white"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())

plot_diff <- effects %>%
  group_by(ntile, ntile.value) %>%
  mutate(`5%`= quantile(mean, probs = c(0.05)),
         `50%` = quantile(mean, probs = c(0.50)),
         `95%` = quantile(mean, probs = c(0.95)),
         mean = mean(mean)) %>%
  ungroup() %>%
  unique() %>%
  mutate(ntile = as.double(ntile)) %>%
  ggplot() +
  geom_point(aes(x = ntile, y = mean), col = "red") +
  geom_errorbar(aes(ymin = `5%`, ymax = `95%`, x = ntile), alpha = 0.1, width = 0.1) +
  xlab(xtitle) + ylab("Treatment Effect") +
  scale_x_continuous(labels = c("Female", "Male"), breaks = 1:length(levels(final.all.extra.vars[,variable])), limits = c(0.5, 2.5))


plot.malesex.diff.marg <- cowplot::plot_grid(
  
  plot_diff
  
  ,
  
  cowplot::plot_grid(
    
    ggplot() +
      theme_void()
    
    ,
    
    plot_hist
    
    , ncol = 2, nrow = 1, rel_widths = c(0.06, 1)
    
  )
  
  , ncol = 1, nrow = 2, rel_heights = c(0.85, 0.15)
  
)


plot_2 <- cowplot::plot_grid(
  plot.egfr.diff.marg,
  plot.hba1c.diff.marg,
  plot.malesex.diff.marg,
  ncol = 3, nrow = 1
)

pdf("Plot2.pdf", width = 9, height = 4)
plot_2
dev.off()


#####

# Plot 3

cred_pred_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/cred_pred_dev.rds")

plot_3 <- cowplot::plot_grid(
  
  cred_pred_dev %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
    geom_point(aes(x = orig, y = mean)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    xlab("Observed HbA1c (mmol/mol)") +
    ylab("Predicted HbA1c (mmol/mol)")
  
  ,
  
  cred_pred_dev %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
    geom_point(aes(x = mean, y = std.resid)) +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    stat_smooth(aes(x = mean, y = std.resid)) +
    xlab("Average Predicted HbA1c (mmol/mol)") +
    ylab("Standardised Residuals")
  
  , ncol = 1, nrow = 2
  
)

pdf("Plot3.pdf", width = 4, height = 7)
plot_3
dev.off()


#####

# Plot 4

source("0.1.slade_functions.R")

# name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

# name: final.val
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_valcohort.Rda")

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/bart_model_final.rds")

data_dev <- final.dev %>%
  # select(-score) %>%
  # left_join(final.all.extra.vars %>%
  #             select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


data_val <- final.val %>%
  # select(-score) %>%
  # left_join(final.all.extra.vars %>%
  #             select(patid, pateddrug, score.excl.mi)) %>% 
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/effects_summary_dev.rds")

effects_summary_val <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/effects_summary_val.rds")



predicted_observed_dev <- data_dev %>%
  cbind(hba1c_diff = effects_summary_dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

predicted_observed_val <- data_val %>%
  cbind(hba1c_diff = effects_summary_val$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10))

ATE_matching_validation_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/ATE_matching_validation_dev.rds")

plot_ATE_dev_prop_score <- ATE_matching_validation_dev[["effects"]] %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = hba1c_diff.pred, y = obs), alpha = 1) + 
  theme_bw() +
  geom_errorbar(aes(x = hba1c_diff.pred, y = obs, ymin = lci, ymax = uci), colour = "black", width = 0.1) +
  ylab("Decile average treatment effect") + 
  xlab("Predicted conditional average treatment effects") +
  ggtitle("Development cohort") +
  scale_x_continuous(limits = c(-15, 15), breaks = c(seq(-15, 15, by = 2))) +
  scale_y_continuous(limits = c(-15, 15), breaks = c(seq(-15, 15, by = 2))) +
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") 

ATE_matching_validation_val <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/ATE_matching_validation_val.rds")

plot_ATE_val_prop_score <- ATE_matching_validation_val[["effects"]] %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = hba1c_diff.pred, y = obs), alpha = 1) + 
  theme_bw() +
  geom_errorbar(aes(x = hba1c_diff.pred, y = obs, ymin = lci, ymax = uci), colour = "black", width = 0.1) +
  ylab("Decile average treatment effect") + 
  xlab("Predicted conditional average treatment effects") +
  ggtitle("Validation cohort") +
  scale_x_continuous(limits = c(-15, 15), breaks = c(seq(-15, 15, by = 2))) +
  scale_y_continuous(limits = c(-15, 15), breaks = c(seq(-15, 15, by = 2))) +
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") 


plot_4 <- cowplot::plot_grid(
  plot_ATE_dev_prop_score,
  plot_ATE_val_prop_score,
  ncol = 2, nrow = 1
)


pdf("Plot4.pdf", width = 8, height = 4)
plot_4
dev.off()


#####

# Plot 5
# How we get predictions for patients with both drugs
# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/bart_model_final.rds")

# patient 39
patient <- rbind(
  final.dev %>% 
    select(c(patid, pateddrug, posthba1c_final, 
             colnames(bart_model_final$X))) %>%
    mutate(drugclass = factor("SGLT2", levels = levels(final.dev$drugclass))) %>%
    slice(39)
  ,
  final.dev %>% 
    select(c(patid, pateddrug, posthba1c_final, 
             colnames(bart_model_final$X))) %>%
    mutate(drugclass = factor("GLP1", levels = levels(final.dev$drugclass))) %>%
    slice(39)
)
  

posteriors <- bartMachine::bart_machine_get_posterior(bart_model_final, 
                                                      patient %>%
                                                        select(colnames(bart_model_final$X)))$y_hat_posterior_samples %>%
  t() %>%
  as.data.frame()
colnames(posteriors) <- c("SGLT2", "GLP1")

plot_5 <- posteriors %>%
  gather() %>%
  ggplot() +
  geom_density(aes(x = value, fill = key), alpha = 0.5, colour = "black") +
  labs(x = "Average HbA1c (mmol/mol)", y = "Density") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.87)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlim(54, 68)


pdf("Plot5.pdf", width = 4.5, height = 4.5)
plot_5
dev.off()

#####

# Plot 6
# Probability of therapy being better than the other
plot_6 <- posteriors %>%
  mutate(hba1c.diff = SGLT2-GLP1) %>%
  select(hba1c.diff) %>%
  mutate(above=ifelse(hba1c.diff > 0, "Favours GLP1", "Favours SGLT2")) %>%
  ggplot(aes(x=hba1c.diff,fill=above)) +
  geom_histogram(position="identity", alpha = 0.5, color="black",breaks=seq(-6,6,by=0.5)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(x="HbA1c difference (mmol/mol)", y = "Density") +
  scale_fill_manual(values=c("red","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.87)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


pdf("Plot6.pdf", width = 4, height = 4)
plot_6
dev.off()


#####

# Plot 7

source("0.1.slade_functions.R")


plot_diff_treatment_response <- function(response, post_hba1c, variable, xtitle, k = 1, ymin = NULL, ymax = NULL, title = "") {
  ##### Input variables
  # response: response summary calculated from calc_diff_treatment_effect function
  # post_hba1c: patient hba1c value post therapy
  # variable: variable being investigated
  # xtitle: title of x axis
  # ymin, ymax: limits of y axis in plot
  
  response$`5%` <- response$`5%` - post_hba1c
  response$`50%` <- response$`50%` - post_hba1c
  response$`95%` <- response$`95%` - post_hba1c
  response$mean <- response$mean - post_hba1c
  
  # load all data for range of variable values; name: final.all.extra.vars
  load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  
  # different approaches whether the variable is continuous or categorical
  if (is.numeric(final.all.extra.vars[, variable])) {
    # if variable is continuous
    
    # plot histogram of all values in variable
    plot_hist <- ggplot() +
      theme_void() +
      geom_histogram(aes(x = final.all.extra.vars[, variable]))
    
    
    
    plot_diff <- response %>%
      ggplot() +
      stat_smooth(aes(x = ntile.value, y = mean, colour = key)) +
      scale_colour_manual(values=c("red","#f1a340")) +
      xlab(xtitle) + ylab("HbA1c response (mmol/mol)") + theme(legend.position = "none")
    
    
    # some variables require logging the x-axis due to extreme values
    if (variable == "preast" | variable == "prebil" | variable == "prealt") {
      
      # log scale of histogram plot
      plot_hist <- plot_hist +
        scale_x_log10()
      
      # log scale of differential response plot
      plot_diff  <- plot_diff +
        scale_x_log10() +
        xlab(paste0(xtitle, " (log)"))
      
    }
    
  } else {
    # if variable is categorical
    
    # plot histogram of all values in variable
    plot_hist <- ggplot() +
      theme_void() +
      geom_bar(aes(x = final.all.extra.vars[, variable]))
    
    
    plot_diff <- response %>%
      ggplot() +
      geom_pointrange(aes(y = mean, ymin = `5%`, ymax = `95%`, x = ntile, colour = key), position = position_dodge2(width=0.5)) +
      scale_colour_manual(values=c("red","#f1a340")) +
      xlab(xtitle) + ylab("HbA1c response (mmol/mol)") + 
      theme(legend.position = "none") +
      scale_x_continuous(labels = levels(final.all.extra.vars[, variable]), breaks = 1:length(levels(final.all.extra.vars[,variable])))
    
    
  }
  
  if (!is.null(ymin) & !is.null(ymax)) {
    plot_diff <- plot_diff +
      ylim(ymin, ymax)
  }
  
  plot_diff <- plot_diff +
    ggtitle(title)
  
  # plot of combined histogram + differential response
  plot.diff.marg <- cowplot::plot_grid(
    
    # plot of differential response
    plot_diff
    
    ,
    
    cowplot::plot_grid(
      
      # spacing of plots
      ggplot() +
        theme_void()
      
      ,
      
      # plot of histogram
      plot_hist
      
      , ncol = 2, nrow = 1, rel_widths = c(0.06, 1)
      
    )
    
    , ncol = 1, nrow = 2, rel_heights = c(0.85, 0.15)
    
  )
  
  return(plot.diff.marg)
  
}

# Differential treatment response for patient: eGFR change

# name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/bart_model_final.rds")

dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))

# patient 39, MRC_plots.R
# patient final hba1c = dataset.dev[9241,"posthba1c_final"] = 54
specific.patient <- cbind(
  drugclass = "SGLT2",
  egfr_ckdepi = as.numeric(dataset.dev[9433, "egfr_ckdepi"]),
  hba1cmonth = 12,
  prealt = as.numeric(dataset.dev[9433, "prealt"]),
  prehba1cmmol = as.numeric(dataset.dev[9433, "prehba1cmmol"]),
  score.excl.mi = as.numeric(dataset.dev[9433, "score.excl.mi"]),
  Category = "Non-smoker",
  drugline = "2",
  ncurrtx = "1",
  yrdrugstart = 2016,
  agetx = as.numeric(dataset.dev[9433, "agetx"]),
  malesex = "0",
  prehdl = as.numeric(dataset.dev[9433, "prehdl"]),
  prebmi = as.numeric(dataset.dev[9433, "prebmi"]),
  prebil = as.numeric(dataset.dev[9433, "prebil"]),
  preplatelets = as.numeric(dataset.dev[9433, "preplatelets"]),
  t2dmduration = as.numeric(dataset.dev[9433, "t2dmduration"]),
  prealb = as.numeric(dataset.dev[9433, "prealb"]),
  presys = as.numeric(dataset.dev[9433, "presys"]),
  preast = as.numeric(dataset.dev[9433, "preast"])
) %>%
  as.data.frame() %>%
  mutate(drugclass = factor(drugclass, levels = levels(dataset.dev$drugclass)),
         egfr_ckdepi = as.numeric(egfr_ckdepi),
         hba1cmonth = as.numeric(hba1cmonth),
         prealt = as.numeric(prealt),
         prehba1cmmol = as.numeric(prehba1cmmol),
         score.excl.mi = as.numeric(score.excl.mi),
         Category = factor(Category, levels = levels(dataset.dev$Category)),
         drugline = factor(drugline, levels = levels(dataset.dev$drugline)),
         ncurrtx = factor(ncurrtx, levels = levels(dataset.dev$ncurrtx)),
         yrdrugstart = as.numeric(yrdrugstart),
         agetx = as.numeric(agetx),
         malesex = factor(malesex, levels = levels(dataset.dev$malesex)),
         prehdl = as.numeric(prehdl),
         prebmi = as.numeric(prebmi),
         prebil = as.numeric(prebil),
         preplatelets = as.numeric(preplatelets),
         t2dmduration = as.numeric(t2dmduration),
         prealb = as.numeric(prealb),
         presys = as.numeric(presys),
         preast = as.numeric(preast)
  )


response_summary_patient <- diff_treatment_response(bart_model_final, specific.patient, 25)


## egfr_ckdepi

# plot treatment response + histogram marginal
plot.egfr.diff.marg.female <- plot_diff_treatment_response(response = response_summary_patient[["egfr_ckdepi"]], 
                                                    post_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                    variable = "egfr_ckdepi", xtitle = "eGFR",
                                                    ymin = -25, ymax = -13, title = "Female")


specific.patient$malesex <- factor("1", levels = levels(dataset.dev$malesex))

response_summary_patient <- diff_treatment_response(bart_model_final, specific.patient, 25)


## egfr_ckdepi

# plot treatment response + histogram marginal
plot.egfr.diff.marg.male <- plot_diff_treatment_response(response = response_summary_patient[["egfr_ckdepi"]], 
                                                           post_hba1c = as.numeric(final.dev[9433, "prehba1cmmol"]),
                                                           variable = "egfr_ckdepi", xtitle = "eGFR",
                                                         ymin = -25, ymax = -13, title = "Male")


plot_7 <- cowplot::plot_grid(
  plot.egfr.diff.marg.female,
  plot.egfr.diff.marg.male
)


pdf("Plot7.pdf", width = 8, height = 5)
plot_7
dev.off()



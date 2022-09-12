####################
## Description:
##  - In this file we plot all the figures used in the presentation
####################

# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])

library(tidyverse)


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
  scale_fill_manual(values=c("#998ec3","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())


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
  scale_fill_manual(values=c("#998ec3","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())


dat3 <- effects_summary_dev_female %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


plot_effects_3 <- ggplot(data=dat3, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-12,18,by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Female",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(values=c("#998ec3","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())


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

# load all data for range of variable values; name: final.all
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# HbA1c

# plot treatment effect + histogram marginal
effects <- effects_summary_dev[["prehba1cmmol"]]
variable <- "prehba1cmmol"
xtitle <- "HbA1c"

plot_hist <- ggplot() +
  theme_classic() +
  geom_histogram(aes(x = final.all[, variable])) +
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
  geom_histogram(aes(x = final.all[, variable])) +
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
  geom_bar(aes(x = final.all[, variable])) +
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
  scale_x_continuous(labels = c("Female", "Male"), breaks = 1:length(levels(final.all[,variable])), limits = c(0.5, 2.5))


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

ATE_validation_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/ATE_validation_dev.rds")

ATE_validation_val <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/Assessment/ATE_validation_val.rds")


#Function to output HTE by subgroup
ATE_plot <- function(data,pred,obs,obslowerci,obsupperci, ymin, ymax, title) {
  ###
  # data: dataset used in fitting,
  # pred: column with predicted values
  # obs: observed values
  # obslowerci: lower bound of CI for prediction
  # obsupperci: upper bound of CI for prediction
  # dataset.type: type of dataset to choose axis: "Development" or "Validation"
  
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_point(alpha=1) + theme_bw() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    xlab("Predicted HbA1c difference (mmol/mol)") + ylab("Observed HbA1c difference (mmol/mol") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle(title) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") 
}



plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -13, 13, "Development cohort")

plot_ATE_val <- ATE_plot(ATE_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -13, 13, "Validation cohort")


plot_4 <- cowplot::plot_grid(
  plot_ATE_dev,
  plot_ATE_val,
  ncol = 2, nrow = 1
)


pdf("Plot4.pdf", width = 7, height = 4)
plot_4
dev.off()


#####

# Plot 5
# How we get predictions for patients with both drugs
# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/bart_model_final.rds")

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
  scale_fill_manual(values=c("#998ec3","#f1a340"))+
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
  scale_fill_manual(values=c("#998ec3","#f1a340"))+
  theme_classic() +
  theme(legend.position = c(0.80, 0.87)) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


pdf("Plot6.pdf", width = 4, height = 4)
plot_6
dev.off()


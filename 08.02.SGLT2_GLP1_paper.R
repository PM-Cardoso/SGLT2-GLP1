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


# name: final.val
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_valcohort.Rda")

effects_summary_val <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/Assessment/effects_summary_val.rds")

## hist_plot function
dat1 <- effects_summary_val %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))



# Plot overall treatment effect benefit
plot_effects_1 <- ggplot(data=dat1, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(plyr::round_any(min(effects_summary_val$mean), 2, f = floor),
                                                                         plyr::round_any(max(effects_summary_val$mean), 2, f = ceiling),
                                                                         by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Overall Population",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(breaks = c("Favours SGLT2", "Favours GLP1"), values = c("#f1a340", "red"))+
  ylim(0, 1300) +
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

# Stratify by sex
effects_summary_val_male <- effects_summary_val %>%
  cbind(malesex = final.val$malesex) %>%
  filter(malesex == 1)

effects_summary_val_female <- effects_summary_val %>%
  cbind(malesex = final.val$malesex) %>%
  filter(malesex == 0)

dat2 <- effects_summary_val_male %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


# Plot treatment effect of male sex
plot_effects_2 <- ggplot(data=dat2, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(plyr::round_any(min(effects_summary_val$mean), 2, f = floor),
                                                                         plyr::round_any(max(effects_summary_val$mean), 2, f = ceiling),
                                                                         by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Strata: Male",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(breaks = c("Favours SGLT2", "Favours GLP1"), values = c("#f1a340", "red"))+
  ylim(0, 1300) +
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

dat3 <- effects_summary_val_female %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))


# Plot treatment effect of female sex
plot_effects_3 <- ggplot(data=dat3, aes(x=mean,fill=above)) +
  geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(plyr::round_any(min(effects_summary_val$mean), 2, f = floor),
                                                                         plyr::round_any(max(effects_summary_val$mean), 2, f = ceiling),
                                                                         by=2)) +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  labs(title="Strata: Female",x="HbA1c difference (mmol/mol)", y = "Number of people") +
  scale_fill_manual(breaks = c("Favours SGLT2", "Favours GLP1"), values = c("#f1a340", "red"))+
  ylim(0, 1300) +
  theme_classic() +
  # theme(legend.position = c(0.80, 0.97)) +
  theme(legend.title = element_blank())

plot_1 <- patchwork::wrap_plots(list(plot_effects_1, plot_effects_2, plot_effects_3), ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(caption = "Figure 1:\nDistribution of predicted individualised treatment effect for SGLT2 treatment compared to GLP1 treatment in the CPRD validation cohort. Negative values reflect a predicted glucose-lowering\ntreatment benefit on SGLT2 treatment, positive values reflect a predicted treatment benefit on GLP1 treatment.") & theme(legend.position = 'bottom',
                                                                                                                                                                                                                                                                                                                                                                         plot.caption = element_text(hjust = 0))

# plot_1


pdf("Plot1.pdf", width = 11, height = 5)
plot_1
dev.off()



# 2. Differential treatment effect, only important vars (eGFR, age, sex)

# name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
# name: final.dev
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")
# name: final.val
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_valcohort.Rda")

source("0.1.slade_functions.R")


plot_diff_treatment_response <- function(response, pre_hba1c, variable, xtitle, k = 1, ymin = NULL, ymax = NULL, title = NULL) {
  ##### Input variables
  # response: response summary calculated from calc_diff_treatment_effect function
  # pre_hba1c: patient hba1c value post therapy
  # variable: variable being investigated
  # xtitle: title of x axis
  # ymin, ymax: limits of y axis in plot
  
  if (variable == "prehba1cmmol") {
    response$`5%` <- response$`5%` - response$ntile.value
    response$`50%` <- response$`50%` - response$ntile.value
    response$`95%` <- response$`95%` - response$ntile.value
    response$mean <- response$mean - response$ntile.value 
  } else {
    response$`5%` <- response$`5%` - pre_hba1c
    response$`50%` <- response$`50%` - pre_hba1c
    response$`95%` <- response$`95%` - pre_hba1c
    response$mean <- response$mean - pre_hba1c
  }
  
  # load all data for range of variable values; name: final.all.extra.vars
  load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  
  # different approaches whether the variable is continuous or categorical
  if (is.numeric(final.all.extra.vars[, variable])) {
    # if variable is continuous
    
    # plot histogram of all values in variable
    plot_hist <- ggplot() +
      theme_void() +
      geom_histogram(aes(x = final.all.extra.vars[, variable]))
    
    if (is.null(title)) {
      
      plot_diff <- response %>%
        ggplot() +
        stat_smooth(aes(x = ntile.value, y = mean, colour = key)) +
        scale_colour_manual(values=c("red","#f1a340")) +
        xlab(xtitle) + ylab("HbA1c response (mmol/mol)") + theme(legend.position = "none",
                                                                 axis.title.x = element_text(face = "bold"))
      
    } else {
      
      plot_diff <- response %>%
        ggplot() +
        stat_smooth(aes(x = ntile.value, y = mean, colour = key)) +
        scale_colour_manual(values=c("red","#f1a340")) +
        xlab(xtitle) + ylab("HbA1c response (mmol/mol)") + theme(legend.position = "none",
                                                                 axis.title.x = element_text(face = "bold")) +
        ggtitle(title)
      
    }
      
    
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
      geom_bar(aes(x = final.all.extra.vars[, variable]), na.rm = TRUE) +
      scale_x_discrete(na.translate = FALSE)
    
    
    plot_diff <- response %>%
      ggplot() +
      geom_pointrange(aes(y = mean, ymin = `5%`, ymax = `95%`, x = ntile, colour = key), position = position_dodge2(width=0.5)) +
      scale_colour_manual(values=c("red","#f1a340")) +
      xlab(xtitle) + ylab("HbA1c response (mmol/mol)") + 
      theme(legend.position = "none",
            axis.title.x = element_text(face = "bold")) +
      scale_x_continuous(labels = c("Female", "Male"), breaks = 1:length(levels(final.all.extra.vars[,variable])))
    
    
  }
  
  if (!is.null(ymin) & !is.null(ymax)) {
    plot_diff <- plot_diff +
      ylim(ymin, ymax)
  }
  
  # plot of combined histogram + differential effects
  plot.diff.marg <- patchwork::wrap_plots(list(plot_diff, plot_hist), ncol = 1, heights = c(0.85, 0.15))
  
  return(plot.diff.marg)
  
}



bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/bart_model_final.rds")

dataset.dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


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
  yrdrugstart = 2015,
  agetx = as.numeric(dataset.dev[9433, "agetx"]),
  malesex = "1",
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


response_summary_patient_male <- diff_treatment_response(bart_model_final, specific.patient, 25)

specific.patient$malesex <- factor("0", levels = levels(dataset.dev$malesex))

response_summary_patient_female <- diff_treatment_response(bart_model_final, specific.patient, 25)


# plot treatment response + histogram marginal
plot.egfr.diff.marg.male <- plot_diff_treatment_response(response = response_summary_patient_male[["egfr_ckdepi"]], 
                                                         pre_hba1c = as.numeric(specific.patient["prehba1cmmol"]),
                                                         variable = "egfr_ckdepi", xtitle = "eGFR",
                                                         ymin = -25, ymax = -9,
                                                         title = "Strata: Male")


# plot treatment response + histogram marginal
plot.egfr.diff.marg.female <- plot_diff_treatment_response(response = response_summary_patient_female[["egfr_ckdepi"]], 
                                                           pre_hba1c = as.numeric(specific.patient["prehba1cmmol"]),
                                                           variable = "egfr_ckdepi", xtitle = "eGFR",
                                                           ymin = -25, ymax = -9,
                                                           title = "Strata: Female")


# plot treatment response + histogram marginal
plot.age.diff.marg.male <- plot_diff_treatment_response(response = response_summary_patient_male[["agetx"]],
                                                        pre_hba1c = as.numeric(specific.patient["prehba1cmmol"]),
                                                        variable = "agetx", xtitle = "Age",
                                                        ymin = -25, ymax = -9,
                                                        title = "Strata: Male")


# plot treatment response + histogram marginal
plot.age.diff.marg.female <- plot_diff_treatment_response(response = response_summary_patient_female[["agetx"]], 
                                                           pre_hba1c = as.numeric(specific.patient["prehba1cmmol"]),
                                                           variable = "agetx", xtitle = "Age",
                                                           ymin = -25, ymax = -9,
                                                          title = "Strata: Female")


plot.malesex.diff.marg <- plot_diff_treatment_response(response = response_summary_patient_male[["malesex"]],
                                                       pre_hba1c = as.numeric(specific.patient["prehba1cmmol"]),
                                                       variable = "malesex", xtitle = "Sex",
                                                       ymin = -25, ymax = -9)


plot_2 <- patchwork::wrap_plots(((plot.egfr.diff.marg.female|plot.egfr.diff.marg.male)/(plot.age.diff.marg.female|plot.age.diff.marg.male))|plot.malesex.diff.marg) +
  patchwork::plot_annotation(caption = "Figure 2:\nDifferential treatment response for a range of values in baseline clinical features. Treatment response estimates for a patient: HbA1c - 76, ALT - 16, eGFR - 98.9, CVD risk - 1.3%, Year of\ndrug start - 2018, Age - 66, HDL - 0.92, BMI - 25.8, Bilirubin - 11, Platelets - 251, Albumin - 43, Systolic - 124, AST - 17, Drug line treatment - 2, number of current drugs - 1,\nSmoking status - Non-smoker, time between diabetes diagnosis and drug indication = 12 years") & theme(plot.caption = element_text(hjust = 0))

plot_2


pdf("Plot2.pdf", width = 11, height = 6)
plot_2
dev.off()


# 3. A validation of outcome, plot: y-defined outcomes, x-ATE

source("0.1.slade_functions.R")

# name: final.all.extra.vars
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")

# # name: final.dev
# load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_devcohort.Rda")

# name: final.val
load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_valcohort.Rda")

bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/bart_model_final.rds")

data_dev <- final.dev %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


data_val <- final.val %>%
  select(-score) %>%
  left_join(final.all.extra.vars %>%
              select(patid, pateddrug, score.excl.mi)) %>%
  select(c(patid, pateddrug, posthba1c_final, 
           colnames(bart_model_final$X)))


# effects_summary_dev <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/Assessment/effects_summary_dev.rds")

effects_summary_val <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/Assessment/effects_summary_val.rds")


predicted_observed_val_initial <- data_val %>%
  cbind(hba1c_diff = effects_summary_val$mean)
  # mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
  # hba1c_diff.q = ntile(hba1c_diff, 10))

group_values <- function(data, variable, breaks) {
  ### Input variables
  # data: dataset used in splitting
  # variable: variable with values to be split
  # breaks: break points between values
  
  # stop in case 'variable' is not included in 'data'
  if (is.null(data[, variable])) {stop("'variable' not included in 'data'")}
  
  # include extra values so that extremes are included
  breaks.full <- c(breaks, floor(min(data[,variable])), ceiling(max(data[,variable])))
  
  new.data <- data %>%
    cbind(hba1c_diff.q = cut(data[, variable], breaks = breaks.full))
  
  return(new.data)
}

breaks = c(-5, 0, 5)

predicted_observed_val <- group_values(data = predicted_observed_val_initial, 
                                       variable = "hba1c_diff", 
                                       breaks = breaks)
vector.test <- predicted_observed_val$hba1c_diff.q

predicted_observed_val$hba1c_diff.q <- as.numeric(predicted_observed_val$hba1c_diff.q)

ATE_weighting_validation_val <- calc_ATE_validation_inverse_prop_weighting(predicted_observed_val, "posthba1c_final")

plot_ATE_dev_prop_score_weighting  <- ATE_plot(ATE_weighting_validation_val[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)

plot_3_SGLT2 <- ATE_weighting_validation_val[["effects"]] %>%
  select(hba1c_diff.q, hba1c_diff.obs, lower.obs, upper.obs) %>%
  mutate(hba1c_diff.q = as.factor(hba1c_diff.q)) %>%
  cbind(best_drug = factor(c("SGLT2", "SGLT2", "GLP1", "GLP1"), levels = c("SGLT2", "GLP1")),
        labels = factor(c("SGLT2 benefit >5 mmol/mol (15.9% CPRD validation set)",
                          "SGLT2 benefit 0-5 mmol/mol (39.0%)",
                          "GLP1 benefit 0-5 mmol/mol (30.0%)",
                          "GLP1 benefit >5 mmol/mol (15.1%)"),
                   levels = c("GLP1 benefit >5 mmol/mol (15.1%)",
                              "GLP1 benefit 0-5 mmol/mol (30.0%)",
                              "SGLT2 benefit 0-5 mmol/mol (39.0%)",
                              "SGLT2 benefit >5 mmol/mol (15.9% CPRD validation set)"))) %>%
  filter(best_drug == "SGLT2") %>%
  ggplot(aes(x = labels, y = hba1c_diff.obs, ymin = lower.obs, ymax = upper.obs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey") +
  scale_y_continuous(limits = c(
    plyr::round_any(min(ATE_weighting_validation_val[["effects"]]$lower.obs), 2, f = floor),
    plyr::round_any(max(ATE_weighting_validation_val[["effects"]]$upper.obs), 2, f = ceiling)),
    breaks = seq(
      plyr::round_any(min(ATE_weighting_validation_val[["effects"]]$lower.obs), 2, f = floor),
      plyr::round_any(max(ATE_weighting_validation_val[["effects"]]$upper.obs), 2, f = ceiling),
      by = 2
    )
  ) +
  geom_pointrange() +
  geom_errorbar(width = 0.2) +
  ylab("Average treatment difference (mmol/mol, negative favours SGLT2i)") +
  ggtitle("Predicted SGLT2i benefit") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_blank()) +
  coord_flip()

plot_3_GLP1 <- ATE_weighting_validation_val[["effects"]] %>%
  select(hba1c_diff.q, hba1c_diff.obs, lower.obs, upper.obs) %>%
  mutate(hba1c_diff.q = as.factor(hba1c_diff.q)) %>%
  cbind(best_drug = factor(c("SGLT2", "SGLT2", "GLP1", "GLP1"), levels = c("SGLT2", "GLP1")),
        labels = factor(c("SGLT2 benefit >5 mmol/mol (15.9% CPRD validation set)",
                          "SGLT2 benefit 0-5 mmol/mol (39.0%)",
                          "GLP1 benefit 0-5 mmol/mol (30.0%)",
                          "GLP1 benefit >5 mmol/mol (15.1%)"),
                        levels = c("GLP1 benefit >5 mmol/mol (15.1%)",
                                   "GLP1 benefit 0-5 mmol/mol (30.0%)",
                                   "SGLT2 benefit 0-5 mmol/mol (39.0%)",
                                   "SGLT2 benefit >5 mmol/mol (15.9% CPRD validation set)"))) %>%
  filter(best_drug == "GLP1") %>%
  ggplot(aes(x = labels, y = hba1c_diff.obs, ymin = lower.obs, ymax = upper.obs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey") +
  scale_y_continuous(limits = c(
    plyr::round_any(min(ATE_weighting_validation_val[["effects"]]$lower.obs), 2, f = floor),
    plyr::round_any(max(ATE_weighting_validation_val[["effects"]]$upper.obs), 2, f = ceiling)),
    breaks = seq(
      plyr::round_any(min(ATE_weighting_validation_val[["effects"]]$lower.obs), 2, f = floor),
      plyr::round_any(max(ATE_weighting_validation_val[["effects"]]$upper.obs), 2, f = ceiling),
      by = 2
    )
  ) +
  geom_pointrange() +
  geom_errorbar(width = 0.2) +
  ylab("Average treatment difference (mmol/mol, negative favours SGLT2i)") +
  ggtitle("Predicted GLP1 benefit") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  coord_flip()

plot_3 <- patchwork::wrap_plots(list(plot_3_SGLT2, plot_3_GLP1), ncol = 1) +
  patchwork::plot_annotation(caption = "Figure 3:\nAverage treatment effects across subgroups defined by clinical cut-offs of predicted treatment benefit. In CPRD, estimates of average\ntreatment effects are calculated using inverse propensity score weighting (IPSW). Figure 1 contains the full distribution of predicted\ntreatment difference estimates.") & theme(plot.caption = element_text(hjust = 0))

plot_3


pdf("Plot3.pdf", width = 8, height = 4)
plot_3
dev.off()






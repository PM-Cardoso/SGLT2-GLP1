####################
## Description:
##  - In this file we use generalised random forests (grf), to model 
##      conditional average treatment effect in a causal model.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])



## increase memery usage to 50gb of RAM
options(java.parameters = "-Xmx50g")

library(tidyverse)
library(bartMachine)


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
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.dev
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_devcohort.Rda"))

load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_valcohort.Rda"))



###############################################################################
###############################################################################
################################ FUNCTIONS ####################################
###############################################################################
###############################################################################

source("0.1.slade_functions.R")


############################# GRF
### Complete model of only routine data, no propensity score (n: 9866))
#############################

data_complete_routine_dev <- final.dev %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass,
    ncurrtx,
    drugline,
    yrdrugstart,
    t2dmduration,
    agetx,
    malesex,
    Category,
    hba1cmonth,
    prebmi,
    prealt,
    egfr_ckdepi,
    prehba1cmmol
  ) %>%
  drop_na() # removed 1302


data_complete_routine_val <- final.val %>%
  select(
    patid,
    pateddrug,
    posthba1c_final,
    drugclass,
    ncurrtx,
    drugline,
    yrdrugstart,
    t2dmduration,
    agetx,
    malesex,
    Category,
    hba1cmonth,
    prebmi,
    prealt,
    egfr_ckdepi,
    prehba1cmmol
  ) %>%
  drop_na() # removed 804


dataset_full <- rbind(data_complete_routine_dev, data_complete_routine_val)

dataset_model.matrix <- model.matrix(~posthba1c_final + drugclass + ncurrtx + drugline + yrdrugstart + t2dmduration + agetx +
                                       malesex + Category + hba1cmonth + prebmi + prealt + egfr_ckdepi + prehba1cmmol, dataset_full) %>%
  as.data.frame() %>%
  select(-`(Intercept)`) %>%
  mutate(drugclass = drugclassSGLT2) %>%
  select(-drugclassSGLT2)


prop.score <- glm(drugclass ~ ncurrtx + drugline + t2dmduration + agetx + 
                    malesex + Category + hba1cmonth + prebmi + prealt + egfr_ckdepi + prehba1cmmol, family = binomial(link = "logit"), data = dataset_full[1:nrow(data_complete_routine_dev),])



grf_model <- grf::causal_forest(X = dataset_model.matrix %>%
                             slice(1:nrow(data_complete_routine_dev)) %>%
                             select(-posthba1c_final, -drugclass),
                           Y = dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"],
                           W = dataset_model.matrix[1:nrow(data_complete_routine_dev), "drugclass"],
                           W.hat = prop.score$fitted.values)


# Calibration of the model
grf.calibration <- grf::test_calibration(grf_model)
# Best linear fit using forest predictions (on held-out data)
# as well as the mean forest prediction as regressors, along
# with one-sided heteroskedasticity-robust (HC3) SEs:
#   
#                                 Estimate Std. Error t value    Pr(>t)
# mean.forest.prediction          0.73197    1.04672  0.6993    0.2422
# differential.forest.prediction  1.31753    0.16636  7.9198 1.339e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Dev
effects.dev <- cbind(mean = grf_model$predictions) %>%
  data.frame() %>%
  set_names(c("mean"))


priority.cate.dev <- 1 * grf_model$predictions

rate.dev <- toc_function(dataset_model.matrix[1:nrow(data_complete_routine_dev),],
                         priority.cate.dev, 
                         prop.score$fitted.values, 
                         grf_model$Y.hat,
                         q = seq(0.1,1,by = 0.05),
                         target = "AUTOC")

# rate.dev$TOC %>%
#   ggplot() +
#   geom_line(aes(x = q, y = estimate)) +
#   geom_line(aes(x = q, y = estimate-1.95*std.err), linetype = "dashed") +
#   geom_line(aes(x = q, y = estimate+1.95*std.err), linetype = "dashed") +
#   ggtitle(paste0("Dev GRF: TOC - ",signif(rate.dev$estimate, 3)," [sd:", signif(rate.dev$std.err, 3),"]"))


#Val
prop.score_val <- predict(prop.score, dataset_full[-c(1:nrow(data_complete_routine_dev)),])
cf.eval <- grf::causal_forest(X = dataset_model.matrix %>%
                           slice(-c(1:nrow(data_complete_routine_dev))) %>%
                           select(-posthba1c_final, -drugclass),
                         dataset_model.matrix[-c(1:nrow(data_complete_routine_dev)), "posthba1c_final"],
                         dataset_model.matrix[-c(1:nrow(data_complete_routine_dev)), "drugclass"],
                         W.hat = prop.score_val)

priority.cate.val <- 1 * cf.eval$predictions


rate.val <- toc_function(dataset_model.matrix[-c(1:nrow(data_complete_routine_dev)),],
                         priority.cate.val, 
                         prop.score_val, 
                         predict(grf_model, dataset_model.matrix %>%
                                   slice(-c(1:nrow(data_complete_routine_dev))) %>%
                                   select(-posthba1c_final, -drugclass)),
                         q = seq(0.1,1,by = 0.05),
                         target = "AUTOC")

# rate.val$TOC %>%
#   ggplot() +
#   geom_line(aes(x = q, y = estimate)) +
#   geom_line(aes(x = q, y = estimate-1.95*std.err), linetype = "dashed") +
#   geom_line(aes(x = q, y = estimate+1.95*std.err), linetype = "dashed") +
#   ggtitle(paste0("Val GRF: TOC - ",signif(rate.val$estimate, 3)," [sd:", signif(rate.val$std.err, 3),"]"))


#######


predicted_observed_complete_routine_dev <- dataset_model.matrix %>%
  slice(1:nrow(data_complete_routine_dev)) %>%
  cbind(hba1c_diff = effects.dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10)) 


ATE_validation_dev <- calc_ATE_validation(predicted_observed_complete_routine_dev %>%
                                            cbind(data_complete_routine_dev[,c("patid", "pateddrug")]))
    
plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -10, 10)



###
# Plot resid

plot_resid_dev <- ggplot() +
  theme_bw() +
  # geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
  geom_point(aes(x = dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"], y = grf_model$Y.hat)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
  xlim(min(dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"], grf_model$Y.hat), max(dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"], grf_model$Y.hat)) +
  ylim(min(dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"], grf_model$Y.hat), max(dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"], grf_model$Y.hat)) +
  xlab("Observed HbA1c (mmol/mol)") +
  ylab("Predicted HbA1c (mmol/mol)")





############


pdf(file = "Plots/3.1.grf_effects.pdf")
prop.score$fitted.values %>%
  as.data.frame() %>%
  set_names(c("value")) %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  ggtitle("Propensity scores")

hist_plot(effects.dev, "Dev GRF: treatment effect", -15, 20)


rate.dev$TOC %>%
  ggplot() +
  geom_line(aes(x = q, y = estimate)) +
  geom_line(aes(x = q, y = estimate-1.95*std.err), linetype = "dashed") +
  geom_line(aes(x = q, y = estimate+1.95*std.err), linetype = "dashed") +
  ggtitle(paste0("Dev GRF: TOC - ",signif(rate.dev$estimate, 3)," [sd:", signif(rate.dev$std.err, 3),"]"))


plot_ATE_dev


plot_resid_dev

dev.off()







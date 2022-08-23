####################
## Description:
##  - In this file we use generalised random forests (grf), to model 
##      conditional average treatment effect in a causal model.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


library(tidyverse)
library(grf)
library(plyr)
library(rms)


## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)



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

source("0.0.slade_functions.R")


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



grf_model <- causal_forest(X = dataset_model.matrix %>%
                             slice(1:nrow(data_complete_routine_dev)) %>%
                             select(-posthba1c_final, -drugclass),
                           Y = dataset_model.matrix[1:nrow(data_complete_routine_dev), "posthba1c_final"],
                           W = dataset_model.matrix[1:nrow(data_complete_routine_dev), "drugclass"],
                           W.hat = prop.score$fitted.values)

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
cf.eval <- causal_forest(X = dataset_model.matrix %>%
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


############


# predicted_observed_complete_routine_dev <- data_complete_routine_dev %>%
#   select(posthba1c_final, prehba1cmmol, drugclass, drugline, ncurrtx, hba1cmonth, egfr_ckdepi, prealt, agetx, prebmi) %>%
#   mutate(prealtlog = log(prealt)) %>%
#   select(-prealt) %>%
#   cbind(hba1c_diff = effects.dev$mean) %>%
#   mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
#          hba1c_diff.q = ntile(hba1c_diff, 10)) 
# 
# predicted_observed_complete_routine_dev_t1 <- predicted_observed_complete_routine_dev %>%
#   ddply("hba1c_diff.q", dplyr::summarise,
#         N = length(hba1c_diff),
#         hba1c_diff.pred = mean(hba1c_diff))
# 
# 
# 
# mnumber = c(1:10)
# models  <- as.list(1:10)
# 
# hba1c_diff.obs.unadj <- vector();lower.unadj <- vector();upper.unadj <- vector();hba1c_diff.obs.sim <- vector()
# lower.sim <- vector();upper.sim <- vector();hba1c_diff.obs.adj <- vector();lower.adj <- vector();upper.adj <- vector() 
# 
# formula1 <- "posthba1c_final~factor(drugclass)"
# formula2 <- "posthba1c_final~factor(drugclass)+prehba1cmmol+ncurrtx+drugline+rcs(hba1cmonth,3)+egfr_ckdepi+prealtlog"
# formula3 <- "posthba1c_final~factor(drugclass)+rcs(prehba1cmmol,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+rcs(egfr_ckdepi,3)+rcs(prealtlog,3)+rcs(agetx,3)+rcs(prebmi,3)"
# f <- as.list(c(formula1,formula2,formula3))
# 
# # #Unadj
# for(i in mnumber) {
#   models[[i]] <- lm(as.formula(formula1),data=predicted_observed_complete_routine_dev,subset=hba1c_diff.q==i)
#   hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,models[[i]]$coefficients[2])
#   confint_all <- confint(models[[i]], levels=0.95)
#   lower.unadj <- append(lower.unadj,confint_all[2,1])
#   upper.unadj <- append(upper.unadj,confint_all[2,2])
# }
# # #Simple 
# for(i in mnumber) {
#   models[[i]] <- lm(as.formula(formula2),data=predicted_observed_complete_routine_dev,subset=hba1c_diff.q==i)
#   hba1c_diff.obs.sim <- append(hba1c_diff.obs.sim,models[[i]]$coefficients[2])
#   confint_all <- confint(models[[i]], levels=0.95)
#   lower.sim <- append(lower.sim,confint_all[2,1])
#   upper.sim <- append(upper.sim,confint_all[2,2])
# }
# #Full
# for(i in mnumber) {
#   models[[i]] <- lm(as.formula(formula3),data=predicted_observed_complete_routine_dev,subset=hba1c_diff.q==i)
#   hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
#   confint_all <- confint(models[[i]], levels=0.95)
#   lower.adj <- append(lower.adj,confint_all[2,1])
#   upper.adj <- append(upper.adj,confint_all[2,2])
# }
# 
# #Final data.frame  
# t1 <- data.frame(predicted_observed_complete_routine_dev_t1,cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj,
#                                                                   hba1c_diff.obs.sim,lower.sim,upper.sim,
#                                                                   hba1c_diff.obs.adj,lower.adj,upper.adj))


predicted_observed_complete_routine_dev <- data_complete_routine_dev %>%
  # select(posthba1c_final, prehba1cmmol, drugclass, drugline, ncurrtx, hba1cmonth, egfr_ckdepi, prealt, agetx, prebmi) %>%
  # mutate(prealtlog = log(prealt)) %>%
  # select(-prealt) %>%
  cbind(hba1c_diff = effects.dev$mean) %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10)) 


t1 <- effects_calibration(predicted_observed_complete_routine_dev)



#simple adj
plotdata_1 <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.unadj,lci=lower.unadj,uci=upper.unadj)
# plotdata_2 <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.sim,lci=lower.sim,uci=upper.sim)
# plotdata_3 <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)


plot_predicted_observed <- hte_plot(plotdata_1,"hba1c_diff.pred","obs","lci","uci") 






############


pdf(paste0(output_path, "/grf_effects.pdf"))
prop.score$fitted.values %>%
  as.data.frame() %>%
  set_names(c("value")) %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  ggtitle("Propensity scores")

hist_plot(effects.dev,-2.5,2.3,1100, "Dev GRF: treatment effect", -15, 20)


rate.dev$TOC %>%
  ggplot() +
  geom_line(aes(x = q, y = estimate)) +
  geom_line(aes(x = q, y = estimate-1.95*std.err), linetype = "dashed") +
  geom_line(aes(x = q, y = estimate+1.95*std.err), linetype = "dashed") +
  ggtitle(paste0("Dev GRF: TOC - ",signif(rate.dev$estimate, 3)," [sd:", signif(rate.dev$std.err, 3),"]"))


plot_predicted_observed


dev.off()







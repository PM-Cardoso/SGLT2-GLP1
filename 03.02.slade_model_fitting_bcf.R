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

############################# BCF
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


dataset_full_bcf <- dataset_model.matrix %>%
  mutate_all(function(x) as.numeric(x)) %>%
  as.matrix()


post <- bcf::bcf(y = dataset_full_bcf[1:nrow(data_complete_routine_dev),1],
            z = dataset_full_bcf[1:nrow(data_complete_routine_dev),19],
            x_control = dataset_full_bcf[1:nrow(data_complete_routine_dev),-c(1,19)],
            pihat = prop.score$fitted.values,
            nburn = 1000,
            nsim = 1000)


effects.dev <- cbind(mean = post$tau %>% colMeans()) %>%
  data.frame() %>%
  set_names(c("mean"))


#########


predicted_observed_complete_routine_dev <- dataset_full_bcf[1:nrow(data_complete_routine_dev), ] %>%
  cbind(hba1c_diff = effects.dev$mean) %>%
  as.data.frame() %>%
  mutate(bestdrug = ifelse(hba1c_diff < 0, "SGLT2", "GLP1"),
         hba1c_diff.q = ntile(hba1c_diff, 10)) 


ATE_validation_dev <- calc_ATE_validation(predicted_observed_complete_routine_dev %>%
                                            cbind(data_complete_routine_dev[,c("patid", "pateddrug")]))

plot_ATE_dev <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -13, 13)



predicted_observed_complete_routine_dev[,"drugclass"][predicted_observed_complete_routine_dev[,"drugclass"] == 1] <- "SGLT2"
predicted_observed_complete_routine_dev[,"drugclass"][predicted_observed_complete_routine_dev[,"drugclass"] == 0] <- "GLP1"

ATE_validation_dev <- calc_ATE_validation_prop_matching(predicted_observed_complete_routine_dev %>%
                                                          cbind(data_complete_routine_dev[,c("patid", "pateddrug")]))

plot_ATE_dev_prop_score <- ATE_plot(ATE_validation_dev[["effects"]], "hba1c_diff.pred", "obs", "lci", "uci", -14, 14)



# Plot
plot_effects_validation <- cowplot::plot_grid(
  
  #title
  cowplot::ggdraw() +
    cowplot::draw_label("Effects Validation")
  
  ,
  
  #effects plot
  plot_ATE_dev
  
  , ncol = 1, nrow = 2, rel_heights = c(0.1, 1)
)


###
# plot residuals


resid_dev <- cbind(lower_bd = apply(post$yhat, MARGIN = 2, function(x) min(x)),
                       upper_bd = apply(post$yhat, MARGIN = 2, function(x) max(x)),
                       mean = apply(post$yhat, MARGIN = 2, function(x) mean(x)),
                       orig = dataset_full_bcf[1:nrow(data_complete_routine_dev),1]) %>%
  as.data.frame() %>%
  mutate(resid = orig - mean,
         resid.low = orig - lower_bd,
         resid.high = orig - upper_bd) 


plot_resid_dev <- resid_dev %>%
  ggplot() +
  theme_bw() +
  geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
  geom_point(aes(x = orig, y = mean)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
  xlim(min(resid_dev$orig), max(resid_dev$orig)) +
  ylim(min(resid_dev$orig), max(resid_dev$orig)) +
  xlab("Observed HbA1c (mmol/mol)") +
  ylab("Predicted HbA1c (mmol/mol)")
  



#########



pdf(file = "Plots/3.2.bcf_effects.pdf")
prop.score$fitted.values %>%
  as.data.frame() %>%
  set_names(c("value")) %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  ggtitle("Propensity scores")

hist_plot(effects.dev, "Dev BCF: treatment effect", -15, 20)

plot_effects_validation

plot_resid_dev

plot_ATE_dev_prop_score

dev.off()









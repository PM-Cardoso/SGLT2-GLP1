####################
## Description:
##  - Functions used throughout the analysis
####################


### Plot of treatment effect histogram, 
###   separate colours for each favoured therapy

hist_plot <- function(data, title, xmin, xmax, xtitle = "HbA1c difference (mmol/mol)", ytitle = "Number of people") {
  ### Input variables
  # data: dataset with column 'mean' corresponding to treatment effect
  # title: title for the plot
  # xmin: lower limit of x axis
  # xmax: upper limit of x axis
  # xtitle: title of x axis
  # ytitle: title of y axis
  
  # define data
  dat <- data %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))
  
  # plot
  plot <- ggplot(data = dat, aes(x = mean, fill = above)) +
    geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(xmin, xmax, by = 1)) +
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    labs(title = title, x = xtitle, y = ytitle) +
    scale_fill_manual(values = c("red", "#f1a340"))+
    theme_classic() +
    theme(legend.position = c(0.80, 0.97)) +
    theme(legend.title = element_blank())
  
  return(plot)
}

### Plot of in-sample prediction vs in-sample/out-sample prediction

hte_plot <- function(data,pred,obs,obslowerci,obsupperci, dataset.type) {
  ### Input variables
  # data: dataset used in fitting,
  # pred: column with predicted values
  # obs: observed values
  # obslowerci: lower bound of CI for prediction
  # obsupperci: upper bound of CI for prediction
  # dataset.type: type of dataset to choose axis: "Development" or "Validation"
  
  # change x-axis title depending on the type of 
  if (dataset.type == "Development") {
    x.axis.title = "In-sample prediction from full model"
  } else if (dataset.type == "Validation") {
    x.axis.title = "Out-of-sample prediction from full model"
  } else {
    stop("'dataset.type' must be 'Development' or 'Validation'")
  }
  
  # plot
  plot <- ggplot(data = data, aes_string(x = pred, y = obs)) +
    geom_point(alpha = 1) + 
    theme_bw() +
    geom_errorbar(aes_string(ymax = obslowerci, ymin = obsupperci), colour = "black", width = .1) +
    ylab("In-sample prediction from sub model") + 
    xlab(x.axis.title) + 
    scale_x_continuous(limits = c(-20, 20), breaks = c(seq(-20, 20, by = 2))) +
    scale_y_continuous(limits = c(-20, 20), breaks = c(seq(-20, 20, by = 2))) +
    geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + 
    geom_hline(yintercept = 0, linetype="dashed", color = "grey60") 
  
  return(plot)
}

### Stability of Treatment effect prediction, sub models fitted to deciles of prediction.

effects_calibration <- function(data, bart_model) {
  ### Input variables
  # data: dataset used in fitting, with columns patid/pateddrug/hba1c_diff.pred
  # bart_model: BART model to be used in the validation
  
  # check whether dataset has "patid" and "pateddrug"
  if ("patid" %in% colnames(data)) {} else {stop("'patid' needs to be included in the dataset")}
  if ("pateddrug" %in% colnames(data)) {} else {stop("'pateddrug' needs to be included in the dataset")}
  if (class(bart_model) != "bartMachine") {stop("'bart_model' needs to be a bartMachine object")}
  
  # split predicted treatment effects into deciles and summarise
  predicted_treatment_effect <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
          N = length(hba1c_diff),
          hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- max(data$hba1c_diff.q)
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  hba1c_diff <- vector()
  lower <- vector()
  upper <- vector();
  
  # iterate through deciles
  for (i in mnumber) {
    # fit decile sub-model
    models[[i]] <- bartMachine::bartMachine(X = data %>%
                                              filter(hba1c_diff.q == i) %>%
                                              select(colnames(bart_model$X)),
                                            y = data %>%
                                              filter(hba1c_diff.q == i) %>%
                                              select("posthba1c_final") %>%
                                              unlist(),
                                            use_missing_data = bart_model$use_missing_data,
                                            num_trees = 50,
                                            num_burn_in = 2000,
                                            num_iterations_after_burn_in = 1000)
    
    # get posteriors for sub-model SGLT2
    effect_dev_SGLT2 <- bartMachine::bart_machine_get_posterior(models[[i]], data %>%
                                                                  filter(hba1c_diff.q == i) %>%
                                                                  select(colnames(bart_model$X)) %>%
                                                                  mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
    
    # get posteriors for sub-model GLP1
    effect_dev_GLP1 <- bartMachine::bart_machine_get_posterior(models[[i]], data %>%
                                                                 filter(hba1c_diff.q == i) %>%
                                                                 select(colnames(bart_model$X)) %>%
                                                                 mutate(drugclass = factor("GLP1", levels = levels(data$drugclass))))
    
    # calculate treatment effect SGLT2-GLP1
    effect_dev <- effect_dev_SGLT2$y_hat_posterior_samples - effect_dev_GLP1$y_hat_posterior_samples %>%
      as.data.frame()
    
    # summarise treatment effects into a single variable + interval
    effects_summary_dev <- cbind(
      mean = apply(effect_dev, MARGIN = 1, function(x) mean(c(x)))
    ) %>%
      as.data.frame()
    
    # pull all values together
    hba1c_diff <- append(hba1c_diff, mean(effects_summary_dev$mean))
    lower <- append(lower, quantile(effects_summary_dev$mean, probs = c(0.05)))
    upper <- append(upper, quantile(effects_summary_dev$mean, probs = c(0.95)))
    
  }
  
  # returned data.frame with results
  final.dataset <- data.frame(predicted_treatment_effect,
                              cbind(hba1c_diff, lower, upper)) %>% 
    dplyr::mutate(obs = hba1c_diff, lci = lower, uci = upper)
  
  return(final.dataset)
}



#  Stability investigation, fitting sub-models to deciles of treatment effect

plot_full_effects_validation <- function(data.dev, data.val, bart_model) {
  ##
  ## This function fits sub-models to deciles of treatment effect
  ##
  ##### Input variables
  # data.dev - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # data.val - Validation dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # bart_model - Model to be used for validation
  
  # calculate effects calibration of Development dataset
  t.dev <- effects_calibration(data = data.dev,
                               bart_model = bart_model)
  
  # plot sub-model effects
  plot_submodel_dev <- hte_plot(t.dev, "hba1c_diff.pred", "obs", "lci", "uci", "Development")
  
  # calculate effects calibration of Validation dataset
  t.val <- effects_calibration(data = data.val,
                               bart_model = bart_model)
  
  # plot sub-model effects
  plot_submodel_val <- hte_plot(t.val, "hba1c_diff.pred", "obs", "lci", "uci", "Validation")
  
  
  # Plot
  plot <- patchwork::wrap_plots(
    # Plot 1
    plot_submodel_dev,
    # Plot 2
    plot_submodel_val
  ) + patchwork::plot_annotation(tag_levels = "A", # labels A = development, B = validation
                                 title = "Effect submodels", # title of full plot
                                 theme = theme(plot.title = element_text(hjust = 0.5))) # center title of full plot
  
  return(plot)
}

## Calculate residuals

calc_resid <- function(data, posteriors, outcome_variable) {
  ##### Input variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  # outcome_variable - variable with outcome values
  
  # calculate standard deviation of residuals
  resid.SD <- apply(posteriors$y_hat_posterior_samples, MARGIN = 2, function(x) (data[,outcome_variable] - x)^2) %>%
    colSums() %>%
    as.data.frame() %>%
    set_names(c("SD")) %>%
    mutate(SD = sqrt(SD/(nrow(data)-2)))
  
  # calculate standardised residuals
  resid <- posteriors$y_hat_posterior_samples
  for (i in 1:nrow(data)) {
    resid[i,] <- (data[i, outcome_variable] - resid[i,])/resid.SD[,1]
  }
  
  # return data.frame with residuals information for each data entry
  cred_pred <- cbind(lower_bd = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) min(x)),
                     upper_bd = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) max(x)),
                     mean = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) mean(x)),
                     orig = data[,outcome_variable]) %>%
    as.data.frame() %>%
    mutate(resid = orig - mean,
           resid.low = orig - lower_bd,
           resid.high = orig - upper_bd) %>%
    cbind(std.resid = apply(resid, MARGIN = 1, function(x) mean(x)),
          std.resid.low = apply(resid, MARGIN = 1, function(x) min(x)),
          std.resid.high = apply(resid, MARGIN = 1, function(x) max(x)))
  
  return(cred_pred)
}


## Calculate assessments of prediction

rsq <- function (x, y) cor(x, y) ^ 2

calc_assessment <- function(data, posteriors, outcome_variable) {
  ##### Input variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  # outcome_variable - variable with y values
  
  # Calculate R2
  r2 <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) rsq(data[,outcome_variable], x)) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  # Calculate RSS: residual sum of squares
  RSS <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) sum((data[,outcome_variable] - x)^2)) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  # Calculate RMSE: root mean square error
  RMSE <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) sqrt(sum((data[,outcome_variable] - x)^2)/nrow(data))) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  # return data.frame with all assessments
  assessment_values <- list(r2 = r2, RSS = RSS, RMSE = RMSE)
  
  return(assessment_values)
}


## Plot predicted vs observed and standardised residuals

resid_plot <- function(pred_dev, pred_val, title) {
  ##### Imput variables
  # pred_dev - predicted/observed values for development dataset
  # pred_val - predicted/observed values for validation dataset
  # title - plot title
  
  # Plot of predicted vs observed for development dataset
  plot_dev_pred <- pred_dev %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
    geom_point(aes(x = orig, y = mean)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
    ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
    xlab("Observed HbA1c (mmol/mol)") +
    ylab("Predicted HbA1c (mmol/mol)")
  
  # Plot of predicted vs observed for validation dataset
  plot_val_pred <- pred_val %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
    geom_point(aes(x = orig, y = mean)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
    ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
    xlab("Observed HbA1c (mmol/mol)") +
    ylab("Predicted HbA1c (mmol/mol)")
  
  # Plot of standardised residuals for development dataset
  plot_dev_std <- pred_dev %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
    geom_point(aes(x = mean, y = std.resid)) +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    stat_smooth(aes(x = mean, y = std.resid)) +
    xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
    ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
    xlab("Average Predicted HbA1c (mmol/mol)") +
    ylab("Standardised Residuals")
  
  # Plot of standardised residuals for validation dataset
  plot_val_std <- pred_val %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
    geom_point(aes(x = mean, y = std.resid)) +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    stat_smooth(aes(x = mean, y = std.resid)) +
    xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
    ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
    xlab("Average Predicted HbA1c (mmol/mol)") +
    ylab("Standardised Residuals")
  
  plot_list <- list(plot_dev_pred, plot_val_pred, plot_dev_std, plot_val_std)
  
  plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = title
    )
  
  return(plot)
  
}


## Calculate treatment effect

calc_effect <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  # get posteriors for SGLT2
  effect_SGLT2 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                            select(
                                                              colnames(bart_model$X)
                                                            ) %>%
                                                            mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
  
  # get posteriors for GLP1
  effect_GLP1 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                           select(
                                                             colnames(bart_model$X)
                                                           ) %>%
                                                           mutate(drugclass = factor("GLP1", levels = levels(data$drugclass))))
  
  # calculate treatment effect for entry
  effect <- effect_SGLT2$y_hat_posterior_samples - effect_GLP1$y_hat_posterior_samples %>%
    as.data.frame()
  
  return(effect)
  
}


## Summarise calculated treatment effect

calc_effect_summary <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  # Calculate treatment effects for entries
  effect <- calc_effect(bart_model, data)
  
  # # Summarise treatment effect for each entry
  effects_summary <- cbind(
    `5%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(effect, MARGIN = 1, function(x) mean(c(x)))
  ) %>%
    as.data.frame()
  
  return(effects_summary)
  
}

## Summarise calculated treatment effect

calc_effect_summary_diff_treat <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  # Calculate treatment effects for entries
  effect <- calc_effect(bart_model, data)
  
  return(effect)
  
}


### Differential treatment effect

# wrapper function for the variables of a model
diff_treatment_effect <- function(bart_model, dataset, rby) {
  ##### Input variables
  # bart_model: bart_model used in model fit
  # dataset: dataset used for differential 
  # rby: number of ntiles
  
  # variables being investigated
  variables <- colnames(bart_model$X)
  
  # number of variables being investigated
  nvars <- length(variables)
  
  # list of differential treatment effects for all variables
  effects.list <- vector("list", length = nvars)
  
  # make names of list elements the same as variables being investigated
  names(effects.list) <- variables
  
  # iterate through the variables
  for (i in 1:nvars) {
    
    # Calculate differential treatment effects
    effects.list[[i]] <- calc_diff_treatment_effect(bart_model, dataset, variables[i], rby)
    
    # print out update on stage of calculation
    print(paste0("Calculation of differential treatment effects for ", variables[i], ": DONE"))
  }
  
  return(effects.list)
  
}

# Calculate differential treatment effect
calc_diff_treatment_effect <- function(bart_model, dataset, variable, rby) {
  ##### Input variables
  # bart_model: bart_model used in model fit
  # dataset: dataset used for differential 
  # variable: variable being investigated
  # rby: number of ntiles
  
  # load all data for range of variable values; name: final.all.extra.vars
  load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  
  
  # different approaches whether the variable is continuous or categorical
  if (is.numeric(dataset[, variable])) {
    # if variable is continuous
    
    # ntile values of variable
    range <- quantile(final.all.extra.vars[,variable], probs = c(seq(0, 1, length.out = rby)), na.rm = TRUE)
    
  } else {
    # if variable is categorical
    
    # possible values of the variable
    range <- levels(final.all.extra.vars[,variable])
    
  }
  
  # new dataset with all values from deciles in dataset
  new.dataset <- dataset
  
  # create a long dataset with all possible combinations of range and dataset
  for (i in 1:length(range)) {
    # first decile only change values
    if (i == 1) {
      new.dataset[,variable] <- range[i]
    # other deciles append new values to full dataset
    } else {
      interim.dataset <- dataset
      interim.dataset[,variable] <- range[i]
      new.dataset <- rbind(new.dataset, interim.dataset)
    }
  }
    
  # if variable is categorical, turn new.dataset variable into factor
  if (!is.numeric(dataset[, variable])) {
    
    # turn new.dataset variable into factor
    new.dataset[,variable] <- factor(new.dataset[,variable], levels = range)
  }
  
  # summary of differential treatment effect
  effects_summary <- calc_effect_summary_diff_treat(bart_model, new.dataset) %>%
    cbind(ntile = rep(1:length(range)),
          ntile.value = rep(range)) %>%
    gather(key, mean, -ntile, -ntile.value)
    
  return(effects_summary)
}


# Plot differential treatment effect

plot_diff_treatment_effect <- function(effects, variable, xtitle, k = 1, thinning = NULL, ymin = NULL, ymax = NULL) {
  ##### Input variables
  # effects: effects summary calculated from calc_diff_treatment_effect function
  # variable: variable being investigated
  # xtitle: title of x axis
  # thinning: thin number of samples plotted
  # ymin, ymax: limits of y axis in plot
  
  if (!is.null(thinning)) {
    
    effects <- effects %>%
      filter(key %in% sample(unique(effects$key), thinning))
    
  }
  
  # load all data for range of variable values; name: final.all.extra.vars
  load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  
  # different approaches whether the variable is continuous or categorical
  if (is.numeric(final.all.extra.vars[, variable])) {
    # if variable is continuous
    
    # plot histogram of all values in variable
    plot_hist <- ggplot() +
      theme_void() +
      geom_histogram(aes(x = final.all.extra.vars[, variable])) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(colour = "white"),
            axis.text.y = element_text(colour = "white"),
            axis.ticks.y = element_line(colour = "white"))
    
    # plot differential treatment effects for the range of values
    plot_diff <- effects %>%
      ggplot() +
      geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") +
      stat_smooth(aes(x = ntile.value, y = mean), col = "red", method = "gam", formula = y ~ s(x, bs = "cs", k=k), level = 0.95) +
      xlab(xtitle) + ylab("Treatment Effect")
      
    
    
    # some variables require logging the x-axis due to extreme values
    if (variable == "preast" | variable == "prebil" | variable == "prealt") {
      
      # log scale of histogram plot
      plot_hist <- plot_hist +
        scale_x_log10()
      
      # log scale of differential effects plot
      plot_diff  <- plot_diff +
        scale_x_log10() +
        xlab(paste0(xtitle, " (log)"))
      
    }
    
  } else {
    # if variable is categorical
    
    # plot histogram of all values in variable
    plot_hist <- ggplot() +
      theme_void() +
      geom_bar(aes(x = final.all.extra.vars[, variable])) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(colour = "white"),
            axis.text.y = element_text(colour = "white"),
            axis.ticks.y = element_line(colour = "white"))
    
    # plot differential treatment effects for the range of values
    plot_diff <- effects %>%
      group_by(ntile) %>%
      mutate(mean.value = mean(mean),
             lower.value = quantile(mean, probs = c(0.05)),
             upper.value = quantile(mean, probs = c(0.95))) %>%
      ungroup() %>%
      mutate(ntile = as.double(ntile)) %>%
      ggplot() +
      geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") +
      geom_point(aes(x = ntile, y = mean.value), col = "red") +
      geom_errorbar(aes(ymin = lower.value, ymax = upper.value, x = ntile), alpha = 0.1) +
      xlab(xtitle) + ylab("Treatment Effect") +
      scale_x_continuous(labels = levels(final.all.extra.vars[, variable]), breaks = 1:length(levels(final.all.extra.vars[,variable])))
    
    
  }
  
  # set limits for treatment effect axis
  if (!is.null(ymin) & !is.null(ymax)) {
    plot_diff <- plot_diff +
      ylim(ymin, ymax)
  }
  
  # plot of combined histogram + differential effects
  plot.diff.marg <- patchwork::wrap_plots(list(plot_diff, plot_hist), ncol = 1, heights = c(0.85, 0.15))

  return(plot.diff.marg)
  
}


# Differential treatment response

# wrapper function for the variables of a model
diff_treatment_response <- function(bart_model, dataset, rby) {
  ##### Input variables
  # bart_model: bart_model used in model fit
  # dataset: dataset used for differential 
  # rby: number of ntiles
  
  # variables being investigated
  variables <- colnames(bart_model$X)
  
  # number of variables being investigated
  nvars <- length(variables)
  
  # list of differential treatment effects for all variables
  effects.list <- vector("list", length = nvars)
  
  # make names of list elements the same as variables being investigated
  names(effects.list) <- variables
  
  # iterate through the variables
  for (i in 1:nvars) {
    
    # Calculate differential treatment effects
    effects.list[[i]] <- calc_diff_treatment_response(bart_model, dataset, variables[i], rby)
    
    # print out update on stage of calculation
    print(paste0("Calculation of differential treatment effects for ", variables[i], ": DONE"))
  }
  
  return(effects.list)
  
}

# Calculate differential treatment response
calc_diff_treatment_response <- function(bart_model, dataset, variable, rby) {
  ##### Input variables
  # bart_model: bart_model used in model fit
  # dataset: dataset used for differential 
  # variable: variable being investigated
  # rby: number of ntiles
  
  # load all data for range of variable values; name: final.all.extra.vars
  load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  
  
  # different approaches whether the variable is continuous or categorical
  if (is.numeric(dataset[, variable])) {
    # if variable is continuous
    
    # ntile values of variable
    range <- quantile(final.all.extra.vars[,variable], probs = c(seq(0, 1, length.out = rby)), na.rm = TRUE)
    
  } else {
    # if variable is categorical
    
    # possible values of the variable
    range <- levels(final.all.extra.vars[,variable])
    
  }
  
  # new dataset with all values from deciles in dataset
  new.dataset <- dataset
  
  # create a long dataset with all possible combinations of range and dataset
  for (i in 1:length(range)) {
    # first decile only change values
    if (i == 1) {
      new.dataset[,variable] <- range[i]
      # other deciles append new values to full dataset
    } else {
      interim.dataset <- dataset
      interim.dataset[,variable] <- range[i]
      new.dataset <- rbind(new.dataset, interim.dataset)
    }
  }
  
  # if variable is categorical, turn new.dataset variable into factor
  if (!is.numeric(dataset[, variable])) {
    
    # turn new.dataset variable into factor
    new.dataset[,variable] <- factor(new.dataset[,variable], levels = range)
  }
  
  # calculate treatment response
  response <- calc_response(bart_model, new.dataset)
  
  # summarise treatment response
  response_summary <- calc_response_summary(response)
  
  # combine into a single dataset
  response_summary_dataset <- rbind(
    response_summary[["SGLT2"]] %>%
      mutate(ntile = rep(1:length(range), each = nrow(dataset)),
             ntile.value = rep(range, each = nrow(dataset)),
             key = "SGLT2"),
    response_summary[["GLP1"]] %>%
      mutate(ntile = rep(1:length(range), each = nrow(dataset)),
             ntile.value = rep(range, each = nrow(dataset)),
             key = "GLP1")
  ) %>%
    as.data.frame()
  
  return(response_summary_dataset)
}


## Calculate treatment response

calc_response <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  # get posteriors for SGLT2
  effect_SGLT2 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                            select(
                                                              colnames(bart_model$X)
                                                            ) %>%
                                                            mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
  
  # get posteriors for GLP1
  effect_GLP1 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                           select(
                                                             colnames(bart_model$X)
                                                           ) %>%
                                                           mutate(drugclass = factor("GLP1", levels = levels(data$drugclass))))
  
  # calculate treatment effect for entry
  response <- list(SGLT2 = effect_SGLT2,
                   GLP1 = effect_GLP1)
  
  return(response)
  
}

## Summarise Treatment response
calc_response_summary <- function(response) {
  ##### Input variables
  # response posteriors
  
  # summarise SGLT2
  summary_SGLT2 <- cbind(
    `5%` = apply(response[["SGLT2"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(response[["SGLT2"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(response[["SGLT2"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(response[["SGLT2"]]$y_hat_posterior_samples, MARGIN = 1, function(x) mean(c(x)))
  ) %>%
    as.data.frame()
  
  # summarise GLP1
  summary_GLP1 <- cbind(
    `5%` = apply(response[["GLP1"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(response[["GLP1"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(response[["GLP1"]]$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(response[["GLP1"]]$y_hat_posterior_samples, MARGIN = 1, function(x) mean(c(x)))
  ) %>%
    as.data.frame()
  
  
  response_summary <- list(SGLT2 = summary_SGLT2,
                           GLP1 = summary_GLP1)
  
  return(response_summary)
}


# Plot differential treatment response

plot_diff_treatment_response <- function(response, pre_hba1c, variable, xtitle, k = 1, ymin = NULL, ymax = NULL) {
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
    
    
    
    plot_diff <- response %>%
      ggplot() +
      stat_smooth(aes(x = ntile.value, y = mean, colour = key)) +
      scale_colour_manual(values=c("red","#f1a340")) +
      xlab(xtitle) + ylab("Predicted HbA1c change (mmol/mol)") + theme(legend.position = "none")
    
    
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
      xlab(xtitle) + ylab("Predicted HbA1c change (mmol/mol)") + 
      theme(legend.position = "none") +
      scale_x_continuous(labels = levels(final.all.extra.vars[, variable]), breaks = 1:length(levels(final.all.extra.vars[,variable])))
    
    
  }
  
  if (!is.null(ymin) & !is.null(ymax)) {
    plot_diff <- plot_diff +
      ylim(ymin, ymax)
  }
  
  # plot of combined histogram + differential effects
  plot.diff.marg <- patchwork::wrap_plots(list(plot_diff, plot_hist), ncol = 1, heights = c(0.85, 0.15))
  
  return(plot.diff.marg)
  
}


# Evaluating ATE from model, unadjusted

ATE_validation <- function(data) {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # bart_model: bart model used for full model in order to take variables used

  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- max(data$hba1c_diff.q)
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector();
  
  # iterate through deciles
  for (i in mnumber) {
    
    # calculate mean of outcome HbA1c in SGLT2
    sglt2.mean <- data %>%
      filter(hba1c_diff.q == i) %>%
      filter(drugclass == "SGLT2") %>%
      select(posthba1c_final) %>%
      colMeans()
    
    # calculate mean of outcome HbA1c in GLP1
    glp1.mean <- data %>%
      filter(hba1c_diff.q == i) %>%
      filter(drugclass == "GLP1") %>%
      select(posthba1c_final) %>%
      colMeans()
    
    # calculate treatment effect SGLT2 - GLP1
    mean.value <- sglt2.mean - glp1.mean
    
    # pull all values together
    hba1c_diff.obs <- append(hba1c_diff.obs, mean.value)
    lower.obs <- append(lower.obs, mean.value)
    upper.obs <- append(upper.obs, mean.value)
  }
  
  # returned data.frame with results 
  t <- data.frame(predicted_treatment_effect,
                  cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>% 
    dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)
  
  return(t)
  
}

# 
# calc_ATE_prop_score <- function(dataset, seed = NULL) {
#   ##### Input variables
#   # bart_model: bart model used for full model in order to take variables used
#   # dataset: dataset for which we calculate propensity scores
#   
#   # extracting selected variables for individuals in dataset
#   data.new <- dataset %>%
#     select(patid, pateddrug) %>%
#     left_join(final.all.extra.vars %>%
#                 select(patid, 
#                        pateddrug,
#                        drugclass,
#                        yrdrugstart,
#                        prebmi,
#                        t2dmduration,
#                        drugline,
#                        prehba1cmmol,
#                        egfr_ckdepi,
#                        ncurrtx,
#                        Category), by = c("patid", "pateddrug"))
#   
#   # fit propensity model with the variables that influence therapy indication
#   prop_model <- bartMachine::bartMachine(X = data.new %>%
#                                            select(yrdrugstart,
#                                                   prebmi,
#                                                   t2dmduration,
#                                                   drugline,
#                                                   prehba1cmmol,
#                                                   egfr_ckdepi,
#                                                   ncurrtx,
#                                                   Category),
#                                          y = data.new[,"drugclass"],
#                                          use_missing_data = TRUE,
#                                          impute_missingness_with_rf_impute = FALSE,
#                                          impute_missingness_with_x_j_bar_for_lm = TRUE,
#                                          num_trees = 200,
#                                          num_burn_in = 1000,
#                                          num_iterations_after_burn_in = 200,
#                                          seed = seed)
# 
#   # keep propensity model
#   prop_model
# 
#   return(prop_model)
# }



calc_ATE_validation <- function(data, variable, prop_model) {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # variable - variable with y values
  # prop_model - propensity score model
  
  # # load all data for range of variable values; name: final.all.extra.vars
  # load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  # 
  # # calculate propensity score
  # prop_model <- calc_ATE_prop_score(data, seed)
  # 
  # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2,
  #   should be the way around)
  prop_score <- 1 - prop_model$p_hat_train
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- max(data$hba1c_diff.q)
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector();
  
  # join dataset and propensity score
  data.new <- data %>%
    cbind(prop_score)
  
  # formula
  formula <- paste0(variable, " ~ factor(drugclass) + prop_score")
  
  # iterate through deciles
  for (i in mnumber) {
    # fit linear regression for decile
    models[[i]] <- lm(as.formula(formula),data=data.new,subset=hba1c_diff.q==i)
    
    # collect treatment effect from regression
    hba1c_diff.obs <- append(hba1c_diff.obs,models[[i]]$coefficients[2])
    
    # calculate confidence intervals
    confint_all <- confint(models[[i]], levels=0.95)
    
    # collect lower bound CI
    lower.obs <- append(lower.obs,confint_all[2,1])
    
    # collect upper bound CI
    upper.obs <- append(upper.obs,confint_all[2,2])
    
  }
  
  
  # join treatment effects for deciles in a data.frame
  effects <- data.frame(predicted_treatment_effect,cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>%
    dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)
  
  # returned list with fitted propensity model + decile treatment effects  
  t <- list(prop_model = prop_model, effects = effects)
  
  return(t)
}


#Function to output ATE by subgroup
ATE_plot <- function(data,pred,obs,obslowerci,obsupperci, ymin, ymax) {
  ###
  # data: dataset used in fitting,
  # pred: column with predicted values
  # obs: observed values
  # obslowerci: lower bound of CI for prediction
  # obsupperci: upper bound of CI for prediction
  # dataset.type: type of dataset to choose axis: "Development" or "Validation"
  
  if (missing(ymin)) {
    ymin <- plyr::round_any(floor(min(data[obslowerci])), 2, f = floor)
  }
  if (missing(ymax)) {
    ymax <- plyr::round_any(ceiling(max(data[obsupperci])), 2, f = ceiling)
  }
  
  # Plot predicted treatment effects vs observed treatment effects
  plot <- ggplot(data = data,aes_string(x = pred,y = obs)) +
    geom_point(alpha = 1) + 
    theme_bw() +
    geom_errorbar(aes_string(ymin = obslowerci, ymax = obsupperci), colour = "black", width = 0.1) +
    ylab("Decile average treatment effect") + 
    xlab("Predicted conditional average treatment effect") +
    scale_x_continuous(limits = c(ymin, ymax), breaks = c(seq(ymin, ymax, by = 2))) +
    scale_y_continuous(limits = c(ymin, ymax), breaks = c(seq(ymin, ymax, by = 2))) +
    geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") 
  
  return(plot)
  
}


### prop matching

calc_ATE_validation_prop_matching <- function(data, variable, prop_model) {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # variable - variable with y values
  # prop_model - propensity score model
  
  # # load all data for range of variable values; name: final.all.extra.vars
  # load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  # 
  # # calculate propensity score
  # prop_model <- calc_ATE_prop_score(data, seed)
  # 
  # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
  prop_score <- 1 - prop_model$p_hat_train
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- max(data$hba1c_diff.q)
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector();
  
  # join dataset and propensity score
  data.new <- data %>%
    cbind(prop_score)
  
  # formula
  formula <- paste0(variable, " ~ factor(drugclass)")
  # iterate through deciles
  for (i in mnumber) {
    
    # dataset individuals in decile that had GLP1
    rows.glp1 <- which(data.new$hba1c_diff.q == i & data.new$drugclass == "GLP1")
    
    # dataset individuals in decile that had SGLT2
    rows.sglt2 <- which(data.new$hba1c_diff.q == i & data.new$drugclass == "SGLT2")
    
    # list of matched SGLT2 rows to GLP1 
    matched.glp1 <- vector(mode = "numeric", length = length(rows.glp1))
    
    # iterate through rows of GLP1
    for (l in 1:length(rows.glp1)) {
      # closest SGLT2 row to GLP1
      chosen.row <- which.min(abs(prop_score[rows.sglt2] - prop_score[rows.glp1[l]]))
      
      # check if distance is less than 0.05 (caliper distance)
      if (prop_score[rows.sglt2[chosen.row]] - prop_score[rows.glp1[l]] < 0.05) {
        # if chosen row is within caliper distance
        
        # update list of matched rows
        matched.glp1[l] <- rows.sglt2[chosen.row]
        
        # remove row from being matched again
        rows.sglt2 <- rows.sglt2[-chosen.row]
        
      } else {
        # if chosen row is outside caliper distance
        
        # update list of matched rows with NA
        matched.glp1[l] <- NA
        
      }
    }
    
    # rows without NA in list of SGLT2 rows matched
    not.na.rows <- !is.na(matched.glp1)
    
    # only keep rows with matched entries
    matched.glp1 <- matched.glp1[not.na.rows]
    rows.glp1 <- rows.glp1[not.na.rows]
    
    # fit linear regression for decile in the matched dataset
    models[[i]] <- lm(as.formula(formula),data=data.new[c(rows.glp1, matched.glp1),],subset=hba1c_diff.q==i)
    
    # collect treatment effect from regression
    hba1c_diff.obs <- append(hba1c_diff.obs,models[[i]]$coefficients[2])
    
    # calculate confidence intervals
    confint_all <- confint(models[[i]], levels=0.95)
    
    # collect lower bound CI
    lower.obs <- append(lower.obs,confint_all[2,1])
    
    # collect upper bound CI
    upper.obs <- append(upper.obs,confint_all[2,2])
    
  }
  
  # join treatment effects for deciles in a data.frame
  effects <- data.frame(predicted_treatment_effect,cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>%
    dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)
  
  # returned list with fitted propensity model + decile treatment effects  
  t <- list(prop_model = prop_model, effects = effects)
  
  return(t)
}


### inverse propensity score weighting 

calc_ATE_validation_inverse_prop_weighting <- function(data, variable, prop_model) {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # variable - variable with y values
  # prop_model - propensity score model
  
  # # load all data for range of variable values; name: final.all.extra.vars
  # load("Samples/SGLT2-GLP1/datasets/cprd_19_sglt2glp1_allcohort.Rda")
  # 
  # # calculate propensity score
  # prop_model <- calc_ATE_prop_score(data, seed)
  # 
  # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
  prop_score <- 1 - prop_model$p_hat_train
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- max(data$hba1c_diff.q)
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  hba1c_diff.obs <- vector(); lower.obs <- vector(); upper.obs <- vector();
  
  # join dataset and propensity score
  data.new <- data %>%
    cbind(calc_prop = prop_score)
  
  # weights for SGLT2 Z = 1
  sglt2.data <- data.new %>%
    filter(drugclass == "SGLT2") %>%
    mutate(calc_prop = 1/(calc_prop))
  
  # weights for GLP1 Z = 0
  glp1.data <- data.new %>%
    filter(drugclass == "GLP1") %>%
    mutate(calc_prop = 1/(1-calc_prop))
  
  data.new <- rbind(sglt2.data, glp1.data)
  
  # formula
  formula <- paste0(variable, " ~ factor(drugclass)")
  
  # iterate through deciles
  for (i in mnumber) {
    # fit linear regression for decile
    models[[i]] <- lm(as.formula(formula),data=data.new,subset=hba1c_diff.q==i, weights = calc_prop)
    
    # collect treatment effect from regression
    hba1c_diff.obs <- append(hba1c_diff.obs,models[[i]]$coefficients[2])
    
    # calculate confidence intervals
    confint_all <- confint(models[[i]], levels=0.95)
    
    # collect lower bound CI
    lower.obs <- append(lower.obs,confint_all[2,1])
    
    # collect upper bound CI
    upper.obs <- append(upper.obs,confint_all[2,2])
    
  }
  
  # join treatment effects for deciles in a data.frame
  effects <- data.frame(predicted_treatment_effect,cbind(hba1c_diff.obs,lower.obs,upper.obs)) %>%
    dplyr::mutate(obs=hba1c_diff.obs,lci=lower.obs,uci=upper.obs)
  
  # returned list with fitted propensity model + decile treatment effects  
  t <- list(prop_model = prop_model, effects = effects)
  
  return(t)
}



###################################################
###################################################
###################################################
###################################################

### Calculate TOC (forked from github/grf-labs/grf)

boot_grf <- function(data, statistic, R, clusters, half.sample = TRUE, ...) {
  samples.by.cluster <- split(seq_along(clusters), clusters)
  n <- length(samples.by.cluster) # number of clusters
  if (n <= 1 || (half.sample && floor(n / 2) <= 1)) {
    stop("Cannot bootstrap sample with only one effective unit.")
  }
  if (half.sample) {
    n.bs <- floor(n / 2)
    index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, n.bs, replace = FALSE)], use.names = FALSE), simplify = FALSE)
  } else {
    index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, replace = TRUE)], use.names = FALSE), simplify = FALSE)
  }
  
  t0 <- statistic(data, seq_len(NROW(data)), ...)
  t0 <- matrix(unlist(t0), ncol = length(t0))
  
  res <- lapply(seq_len(R), function(i) statistic(data, index.list[[i]], ...))
  t <- matrix(, R, length(t0))
  for (r in seq_len(R)) {
    t[r, ] <- unlist(res[[r]])
  }
  
  list(t0 = t0, t = t)
}

estimate_rate <- function(data, indices, q, wtd.mean) {
  prio <- data[indices, 3]
  sort.idx <- order(prio, decreasing = TRUE)
  sample.weights <- data[indices, 2][sort.idx]
  
  num.ties <- tabulate(prio) # count by rank in increasing order
  num.ties <- num.ties[num.ties != 0] # ignore potential ranks not present in BS sample
  # if continuous scoring then no need for ~slower aggregation
  if (all(num.ties == 1)) {
    DR.scores.sorted <- (data[indices, 1] / data[indices, 2])[sort.idx]
  } else {
    grp.sum <- rowsum(data[indices, 1:2][sort.idx, ], prio[sort.idx], reorder = FALSE)
    DR.avg <- grp.sum[, 1] / grp.sum[, 2]
    DR.scores.sorted <- rep.int(DR.avg, rev(num.ties))
  }
  sample.weights.cumsum <- cumsum(sample.weights)
  sample.weights.sum <- sample.weights.cumsum[length(sample.weights)]
  ATE <- sum(DR.scores.sorted * sample.weights) / sample.weights.sum
  TOC <- cumsum(DR.scores.sorted * sample.weights) / sample.weights.cumsum - ATE
  RATE <- wtd.mean(TOC, sample.weights)
  
  nw <- q * sample.weights.sum
  idx <- findInterval(nw + 1e-15, sample.weights.cumsum) # epsilon tol. since finite precision may cause not equal
  denominator.adj <- nw - sample.weights.cumsum[pmax(idx, 1)] # \sum weight remainder, zero when no fractional obs. needed
  numerator.adj <- denominator.adj * DR.scores.sorted[pmin(idx + 1, max(idx))]
  idx[idx == 0] <- 1
  uniq.idx <- unique(idx)
  grid.id <- rep.int(seq_along(uniq.idx), c(uniq.idx[1], diff(uniq.idx)))
  DR.scores.grid <- rowsum(cbind(DR.scores.sorted * sample.weights, sample.weights), grid.id, reorder = FALSE)
  TOC.grid <- (cumsum(DR.scores.grid[, 1])[grid.id[idx]] + numerator.adj) /
    (cumsum(DR.scores.grid[, 2])[grid.id[idx]] + denominator.adj) - ATE
  
  c(RATE, TOC.grid, use.names = FALSE)
}

## Function for TOC
toc_function <- function(data, tau, W.hat, Y.hat, R = 200, q = seq(0.1, 1, by = 0.1), target = c("AUTOC", "QINI")) {
  subset.clusters <- 1:nrow(data)
  subset.weights <- rep(1, nrow(data))/ sum( rep(1, nrow(data)))
  priorities <- as.data.frame(tau, fix.empty.names = FALSE)
  empty.names <- colnames(priorities) == ""
  colnames(priorities)[empty.names] <- c("priority1", "priority2")[1:ncol(priorities)][empty.names]
  priorities[,] <- lapply(priorities, function(x) as.integer(as.factor(x)))
  
  W.orig <- model.matrix(~drugclass, data) %>% as.data.frame() %>%select(starts_with("drugclass")) %>% unlist()
  Y.orig <- data[,"posthba1c_final"]
  
  debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
  Y.residual <- Y.orig - (Y.hat + tau * (W.orig - W.hat))
  
  DR.scores <- tau + debiasing.weights * Y.residual
  
  if (target == "AUTOC") {
    wtd.mean <- function(x, w) sum(x * w) / sum(w)
  } else if (target == "QINI") {
    wtd.mean <- function(x, w) sum(cumsum(w) / sum(w) * w * x) / sum(w)
  }
  
  
  boot.output <- boot_grf(
    data = data.frame(DR.scores * subset.weights, subset.weights, priorities),
    # In case of two priorities do a paired bootstrap estimating both prios on same sample.
    statistic = function(data, indices, q, wtd.mean)
      lapply(c(4, 3)[1:ncol(priorities)], function(j) estimate_rate(data[, -j], indices, q, wtd.mean)),
    R = R,
    clusters = subset.clusters,
    half.sample = TRUE,
    q = q,
    wtd.mean = wtd.mean
  )
  
  dim(boot.output[["t"]]) <- c(R, dim(boot.output[["t0"]]))
  point.estimate <- boot.output[["t0"]]
  std.errors <- apply(boot.output[["t"]], c(2, 3), sd)
  if (ncol(priorities) > 1) {
    point.estimate <- cbind(point.estimate, point.estimate[, 1] - point.estimate[, 2])
    std.errors <- cbind(std.errors, apply(boot.output[["t"]][,, 1] - boot.output[["t"]][,, 2], 2, sd))
  }
  point.estimate[abs(point.estimate) < 1e-15] <- 0
  std.errors[abs(std.errors) < 1e-15] <- 0
  
  if (R < 2) {
    std.errors[] <- 0
  }
  priority <- c(colnames(priorities), paste(colnames(priorities), collapse = " - "))[1:length(point.estimate[1, ])]
  
  output <- list()
  class(output) <- "rank_average_treatment_effect"
  output[["estimate"]] <- point.estimate[1, ]
  output[["std.err"]] <- std.errors[1, ]
  output[["target"]] <- paste(priority, "|", target)
  output[["TOC"]] <- data.frame(estimate = c(point.estimate[-1, ]),
                                std.err = c(std.errors[-1, ]),
                                q = q,
                                priority = priority[rep(1:length(priority), each = length(q))],
                                stringsAsFactors = FALSE)
  
  return(output)
  
}


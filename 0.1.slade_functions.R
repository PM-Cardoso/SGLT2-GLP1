####################
## Description:
##  - Functions used throughout the analysis
####################


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


### Plot of treatment effect histogram

# Plots
hist_plot <- function(data, title, xmin, xmax) {
  ###
  # data: dataset with column 'mean' corresponding to treatment effect
  # title: title for the plot
  # xmin: lower limit of x axis
  # xmax: upper limit of x axis
  
  #define data
  dat <- data %>% dplyr::select(mean) %>% mutate(above=ifelse(mean> 0, "Favours GLP1", "Favours SGLT2"))
  c_low <- quantile(dat$mean,.001)
  c_upp <- quantile(dat$mean,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  c_low <- min(dat$mean,.001)
  c_upp <- quantile(dat$mean,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  
  #plot
  ggplot(data=dat, aes(x=mean,fill=above)) +
    geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(xmin,xmax,by=1)) +
    geom_vline(aes(xintercept=0), linetype="dashed")+
    labs(title=title,x="HbA1c difference (mmol/mol)", y = "Number of people") +
    scale_fill_manual(values=c("#998ec3","#f1a340"))+
    theme_classic() +
    theme(legend.position = c(0.80, 0.97)) +
    theme(legend.title = element_blank())
}


### Plots of predicted vs observed treatment effect

#Function to output HTE by subgroup
hte_plot <- function(data,pred,obs,obslowerci,obsupperci) {
  ###
  # data: dataset used in fitting,
  # pred: column with predicted values
  # obs: observed values
  # obslowerci: lower bound of CI for prediction
  # obsupperci: upper bound of CI for prediction
  
  
  #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
  ymin  <- -20;  ymax <- 20
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_point(alpha=1) + theme_bw() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    ylab("Q: Predicted HbA1c difference (mmol/mol)") + xlab("Predicted HbA1c difference (mmol/mol)") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") 
}


### Treatment effect calibration for predicted vs observed treatment effect

effects_calibration <- function(data, bart_model) {
  ###
  # data: dataset used in fitting, with columns patid/pateddrug/hba1c_diff.pred
  # bart_model: BART model to be used in the validation
  
  # check whether dataset has "patid" and "pateddrug"
  if ("patid" %in% colnames(data)) {} else {stop("'patid' needs to be included in the dataset")}
  if ("pateddrug" %in% colnames(data)) {} else {stop("'pateddrug' needs to be included in the dataset")}
  if (class(bart_model) != "bartMachine") {stop("'bart_model' needs to be a bartMachine object")}
  
  # split predicted treatment effects into deciles
  predicted_observed_complete_routine <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
          N = length(hba1c_diff),
          hba1c_diff.pred = mean(hba1c_diff))
  
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.unadj <- vector(); lower.unadj <- vector(); upper.unadj <- vector();
  
  for (i in mnumber) {
    # fit decile model
    models[[i]] <- bartMachine::bartMachine(X = data %>%
                                              filter(hba1c_diff.q == i) %>%
                                              select(colnames(bart_model$X)),
                                            y = data %>%
                                              filter(hba1c_diff.q == i) %>%
                                              select("posthba1c_final") %>%
                                              unlist(),
                                            use_missing_data = bart_model$use_missing_data,
                                            impute_missingness_with_rf_impute = bart_model$impute_missingness_with_rf_impute,
                                            impute_missingness_with_x_j_bar_for_lm = bart_model$impute_missingness_with_x_j_bar_for_lm,
                                            num_trees = 50,
                                            num_burn_in = 2000,
                                            num_iterations_after_burn_in = 1000)
    
    effect_dev_SGLT2 <- bartMachine::bart_machine_get_posterior(models[[i]], data %>%
                                                                  filter(hba1c_diff.q == i) %>%
                                                                  select(colnames(bart_model$X)) %>%
                                                                  mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
    
    effect_dev_GLP1 <- bartMachine::bart_machine_get_posterior(models[[i]], data %>%
                                                                 filter(hba1c_diff.q == i) %>%
                                                                 select(colnames(bart_model$X)) %>%
                                                                 mutate(drugclass = factor("GLP1", levels = levels(data$drugclass))))
    
    effect_dev <- effect_dev_SGLT2$y_hat_posterior_samples - effect_dev_GLP1$y_hat_posterior_samples %>%
      as.data.frame()
    
    effects_summary_dev <- cbind(
      `5%` = apply(effect_dev, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
      `50%` = apply(effect_dev, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
      `95%` = apply(effect_dev, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
      mean = apply(effect_dev, MARGIN = 1, function(x) mean(c(x)))
    ) %>%
      as.data.frame()
    
    hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,mean(effects_summary_dev$mean))
    lower.unadj <- append(lower.unadj,mean(effects_summary_dev$`5%`))
    upper.unadj <- append(upper.unadj,mean(effects_summary_dev$`95%`))
    
    
  }
  
  #Final data.frame  
  t <- data.frame(predicted_observed_complete_routine,
                   cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj)) %>% 
    dplyr::mutate(obs=hba1c_diff.obs.unadj,lci=lower.unadj,uci=upper.unadj)
  
  return(t)
}

## Calculate residuals

calc_resid <- function(data, posteriors) {
  ##### Imput variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  
  resid.SD <- apply(posteriors$y_hat_posterior_samples, MARGIN = 2, function(x) (data$posthba1c_final - x)^2) %>%
    colSums() %>%
    as.data.frame() %>%
    set_names(c("SD")) %>%
    mutate(SD = sqrt(SD/(nrow(data)-2)))
  
  resid <- posteriors$y_hat_posterior_samples
  for (i in 1:nrow(data)) {
    resid[i,] <- (data$posthba1c_final[i] - resid[i,])/resid.SD[,1]
  }
  
  cred_pred <- cbind(lower_bd = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) min(x)),
                     upper_bd = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) max(x)),
                     mean = apply(posteriors$y_hat_posterior_samples, MARGIN = 1, function(x) mean(x)),
                     orig = data[,"posthba1c_final"]) %>%
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

calc_assessment <- function(data, posteriors) {
  ##### Imput variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  
  r2 <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) rsq(data[,"posthba1c_final"], x)) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  RSS <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) sum((data[,"posthba1c_final"] - x)^2)) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  RMSE <- posteriors$y_hat_posterior_samples %>%
    apply(MARGIN = 2, function(x) sqrt(sum((data[,"posthba1c_final"] - x)^2)/nrow(data))) %>%
    quantile(probs = c(0.05, 0.5, 0.95))
  
  assessment_values <- list(r2 = r2, RSS = RSS, RMSE = RMSE)
  
  return(assessment_values)
}

## Plot predicted vs observed and standardised residuals

resid_plot <- function(pred_dev, pred_val, title) {
  ##### Imput variables
  # pred_dev - predicted/observed values for development dataset
  # pred_val - predicted/observed values for validation dataset
  # title - plot title
  
  cowplot::plot_grid(
    
    #title
    cowplot::ggdraw() +
      cowplot::draw_label(title)
    
    ,
    
    cowplot::plot_grid(
      
      pred_dev %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
        geom_point(aes(x = orig, y = mean)) +
        geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        xlab("Observed HbA1c (mmol/mol)") +
        ylab("Predicted HbA1c (mmol/mol)")
      
      ,
      
      pred_val %>%
        ggplot() +
        theme_bw() +
        geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
        geom_point(aes(x = orig, y = mean)) +
        geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
        xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
        xlab("Observed HbA1c (mmol/mol)") +
        ylab("Predicted HbA1c (mmol/mol)")
      
      ,
      
      pred_dev %>%
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
      
      ,
      
      pred_val %>%
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
      
      , ncol = 2, nrow = 2, labels = c("A", "B", "", "")
      
    )
    
    , ncol = 1, nrow = 2, rel_heights = c(0.1,1)
    
  )

}

## Calculate treatment effect

calc_effect <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  
  effect_SGLT2 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                            select(
                                                              colnames(bart_model$X)
                                                            ) %>%
                                                            mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
  
  effect_GLP1 <- bartMachine::bart_machine_get_posterior(bart_model, data %>%
                                                           select(
                                                             colnames(bart_model$X)
                                                           ) %>%
                                                           mutate(drugclass = factor("GLP1", levels = levels(data$drugclass))))
  
  effect <- effect_SGLT2$y_hat_posterior_samples - effect_GLP1$y_hat_posterior_samples %>%
    as.data.frame()
  
  return(effect)
  
}


calc_effect_summary <- function(bart_model, data) {
  ##### Input variables
  # bart_model - bart model used for fitting
  # data - data being investigated
  
  effect <- calc_effect(bart_model, data)
  
  effects_summary <- cbind(
    `5%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
    `50%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
    `95%` = apply(effect, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
    mean = apply(effect, MARGIN = 1, function(x) mean(c(x)))
  ) %>%
    as.data.frame()
  
  
}

plot_full_effects_validation <- function(data.dev, data.val, bart_model) {
  ##### Input variables
  # data.dev - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # data.val - Validation dataset with variables + treatment effect quantiles (hba1c_diff.q)
  # bart_model - Model to be used for validation
  
  # Effects calibration of Development dataset
  t.dev <- effects_calibration(data = data.dev,
                               bart_model = bart_model)
  
  plot_predicted_observed_dev <- hte_plot(t.dev, "hba1c_diff.pred", "obs", "lci", "uci")
  
  # Effects calibration of Validation dataset
  t.val <- effects_calibration(data = data.val,
                               bart_model = bart_model)
  
  plot_predicted_observed_val <- hte_plot(t.val, "hba1c_diff.pred", "obs", "lci", "uci")
  
  
  # Plot
  plot <- cowplot::plot_grid(
    
    #title
    cowplot::ggdraw() +
      cowplot::draw_label("Effects Validation")
    
    ,
    
    #effects plot
    cowplot::plot_grid(
      
      # Development plot
      plot_predicted_observed_dev
      
      ,
      
      # Validation plot
      plot_predicted_observed_val
      
      , ncol = 2, nrow = 1, labels = c("A", "B")
    )
    
    , ncol = 1, nrow = 2, rel_heights = c(0.1, 1)
  )
  
  return(plot)
}

### Instability assessment

instability_modelling <- function(bart_model, m, n) {
  ##### Input variables
  # bart_model: bart_model being investigated for instability
  # m: number of subsets from the dataset
  # n: number of rows included in each subset of the dataset
  
  # set a list of subsets of length m
  subsets_rows <- as.list(1:m)
  
  
  for (i in 1:m) {
    # total number of rows
    n_ind <- nrow(bart_model$X)
    # sample n rows from original dataset
    rows_ind <- sample(1:n_ind, n)
    # save subset rows
    subsets_rows[[i]] <- rows_ind
  }
  
  # set a list of bartMachine models for subsets
  subsets_models <- as.list(1:m)
  
  
  for (i in 1:m) {
    bart_sub_model <- bartMachine::bartMachine(X = bart_model$X[subsets_rows[[i]],],
                                               y = bart_model$y[subsets_rows[[i]]],
                                               use_missing_data = bart_model$use_missing_data,
                                               impute_missingness_with_rf_impute = bart_model$impute_missingness_with_rf_impute,
                                               impute_missingness_with_x_j_bar_for_lm = bart_model$impute_missingness_with_x_j_bar_for_lm,
                                               num_trees = bart_model$use_missing_data,
                                               num_burn_in = bart_model$num_burn_in,
                                               num_iterations_after_burn_in = bart_model$num_iterations_after_burn_in,
                                               serialize = TRUE
    )
    subsets_model[[i]] <- bart_sub_model
  }
  
  subset_modelling <- list(
    original = bart_model,
    m = m,
    n = n,
    rows = subsets_rows,
    models = subsets_models
  )
  
  return(subset_modelling)
}

instability_predictor_vs_observed_plot <- function(bart_model_subset) {
  ##### Input variables
  # bart_model_subset: list of rows and models used in the instability assessment
  
  p <- ggplot() +
    theme_bw() +
    stat_smooth(aes(x = bart_model_subset[["original"]]$y, 
                    y = bart_model_subset[["original"]]$y_hat_train),
                method = "lm", colour = "blue")
  
  for (i in 1:bart_model_subset[["m"]]) {
    p <- p + stat_smooth(aes(x = bart_model_subset[["models"]][[i]]$y, 
                             y = bart_model_subset[["models"]][[i]]$y_hat_train),
                         method = "lm", colour = "blue", alpha = 0.3)
  }
  
  p <- p + geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), lwd=0.75) +
    xlim(min(bart_model_subset[["original"]]$y,
             bart_model_subset[["original"]]$y_hat_train), 
         max(bart_model_subset[["original"]]$y,
             bart_model_subset[["original"]]$y_hat_train)) +
    ylim(min(bart_model_subset[["original"]]$y,
             bart_model_subset[["original"]]$y_hat_train), 
         max(bart_model_subset[["original"]]$y,
             bart_model_subset[["original"]]$y_hat_train)) +
    xlab("Observed HbA1c (mmol/mol)") +
    ylab("Predicted HbA1c (mmol/mol)")
}



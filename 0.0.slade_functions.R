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
hist_plot <- function(data,sx,sy,y, title, xmin, xmax) {
  #label for hist
  annotation <- data.frame(
    x = c(sx,sy),
    y = c(y),
    label = c("Favours SGLT2", "Favours GLP1")
  )
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
  
  #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
  ymin  <- -15;  ymax <- 15
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_point(alpha=1) + theme_bw() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    ylab("Observed HbA1c difference (mmol/mol)") + xlab("Predicted HbA1c difference (mmol/mol)") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") 
}


### Treatment effect calibration for predicted vs observed treatment effect

effects_calibration <- function(data, model = "BART", formula = "formula1", dataset = NULL) {
  ###
  # data: dataset used in fitting, with columns patid/pateddrug/hba1c_diff.pred
  # model: model used to calculate observed treatment effect: BARt or linear regression
  # formula: variables used in the calculation of observed effect
  
  
  # check whether dataset has "patid" and "pateddrug"
  if ("patid" %in% colnames(data)) {} else {stop("'patid' needs to be included in the dataset")}
  if ("pateddrug" %in% colnames(data)) {} else {stop("'pateddrug' needs to be included in the dataset")}
  if (is.null(dataset)) {stop("'dataset' needs to be specified whether dataset used is Dev or Val")}
  
  if (formula == "formula1") {
    vars_selected <- c("patid", "pateddrug", "posthba1c_final", "drugclass")
    if (model == "Linear") {
      linear_formula <- "posthba1c_final~factor(drugclass)"
    }
  } else if (formula == "formula2") {
    vars_selected <- c("patid", "pateddrug", "posthba1c_final", "drugclass", "prehba1cmmol", "ncurrtx", "drugline", "hba1cmonth", "egfr_ckdepi", "prealt")
    if (model == "Linear") {
      linear_formula <- "posthba1c_final~factor(drugclass)+prehba1cmmol+ncurrtx+drugline+rms::rcs(hba1cmonth,3)+egfr_ckdepi+log(prealt)"
    }
  } else if (formula == "formula3") {
    vars_selected <- c("patid", "pateddrug", "posthba1c_final", "drugclass", "prehba1cmmol", "ncurrtx", "drugline", "hba1cmonth", "egfr_ckdepi", "prealt", "agetx", "prebmi")
    if (model == "Linear") {
      linear_formula <- "posthba1c_final~factor(drugclass)+rms::rcs(prehba1cmmol,3)+ncurrtx+drugline+rms::rcs(hba1cmonth,3)+rms::rcs(egfr_ckdepi,3)+rms::rcs(log(prealt),3)+rms::rcs(agetx,3)+rms::rcs(prebmi,3)"
    }
  }
  
  
  if (dataset == "Dev") {
    patient_values <- data %>%
      select(patid, pateddrug, hba1c_diff.q) %>%
      left_join(final.dev %>%
                  select(vars_selected), by = c("patid", "pateddrug"))
  }
  
  
  predicted_observed_complete_routine_dev_t1 <- data %>%
    plyr::ddply("hba1c_diff.q", dplyr::summarise,
          N = length(hba1c_diff),
          hba1c_diff.pred = mean(hba1c_diff))
  
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.unadj <- vector();lower.unadj <- vector();upper.unadj <- vector();
  
  
  if (model == "BART") {
    
    # #Unadj
    for (i in mnumber) {
      models[[i]] <- bartMachine::bartMachine(X = patient_values %>%
                                                filter(hba1c_diff.q == i) %>%
                                                select(-hba1c_diff.q,
                                                       -patid,
                                                       -pateddrug,
                                                       -posthba1c_final),
                                              y = patient_values %>%
                                                filter(hba1c_diff.q == i) %>%
                                                select(posthba1c_final) %>%
                                                unlist(),
                                              use_missing_data = TRUE,
                                              impute_missingness_with_rf_impute = FALSE,
                                              impute_missingness_with_x_j_bar_for_lm = TRUE,
                                              num_trees = 50,
                                              num_burn_in = 300,
                                              num_iterations_after_burn_in = 100
      )
      effect_dev_SGLT2 <- bartMachine::bart_machine_get_posterior(models[[i]], patient_values %>%
                                                       filter(hba1c_diff.q == i) %>%
                                                       select(-hba1c_diff.q,
                                                              -patid,
                                                              -pateddrug,
                                                              -posthba1c_final) %>%
                                                       mutate(drugclass = factor("SGLT2", levels = levels(data$drugclass))))
      
      effect_dev_GLP1 <- bartMachine::bart_machine_get_posterior(models[[i]], patient_values %>%
                                                      filter(hba1c_diff.q == i) %>%
                                                      select(-hba1c_diff.q,
                                                             -patid,
                                                             -pateddrug,
                                                             -posthba1c_final) %>%
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
    
  } else {
    # When model == Linear
    
    for(i in mnumber) {
      models[[i]] <- lm(as.formula(linear_formula),data=patient_values,subset=hba1c_diff.q==i)
      hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,models[[i]]$coefficients[2])
      confint_all <- confint(models[[i]], levels=0.95)
      lower.unadj <- append(lower.unadj,confint_all[2,1])
      upper.unadj <- append(upper.unadj,confint_all[2,2])
    }
    
  }
  
  #Final data.frame  
  t1 <- data.frame(predicted_observed_complete_routine_dev_t1,cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj))
  
  return(t1)
}


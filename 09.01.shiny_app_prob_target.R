####################
## Description:
##  - In this file we make a Shiny App for predicting the probability of 
##      achieving a target HbA1c.
####################

# if probability says 100%, change to >99%
# find patients with several outcomes for examples



## increase memory usage to 3gb of RAM
options(java.parameters = "-Xmx3000m")

library(bartMachine)
require(tidyverse)
require(shiny)

# Define UI for inputting patient values
ui <- fluidPage(
  
  # App title
  titlePanel("HbA1c prediction"),
  
  fluidRow(
    column(3,
           h3("Clinical information:"),
           h4("HbA1c target (mmol/mol)"),
           textInput("target_num",
                     label = "",
                     value = "52",
                     width = "30%"),
           br(),
           actionButton("calculate", "Calculate"),
           
           br(),
           
           radioButtons("output_type",
                        label = h3("Output type:"),
                        choices = list("Absolute probability" = 1,
                                       "Outcome interval" = 2,
                                       "Probability density" = 3),
                        selected = 1)
           ),
    
    column(2,
           selectInput("sex_select",
                       label = h6("Sex"),
                       choices = list("",
                                      "Male",
                                      "Female",
                                      "Missing"),
                       selected = NULL),
           numericInput("egfr_num",
                        label = h6("eGFR (ml/min/1.73 m2)"),
                        value = 96,
                        min = 40,
                        max = 170),
           numericInput("alt_num",
                        label = h6("ALT (U/L)"),
                        value = 23,
                        min = 2,
                        max = 4100),
           numericInput("hba1c_num",
                        label = h6("HbA1c (mmol/mol)"),
                        value = 65,
                        min = 50,
                        max = 150)
          ),
    column(2,
           selectInput("smoker_select",
                       label = h6("Smoking status"),
                       choices = list("",
                                      "Non-smoker",
                                      "Ex-smoker",
                                      "Active smoker",
                                      "Missing"),
                       selected = NULL),
           numericInput("age_num",
                        label = h6("Age (years)"),
                        value = 65,
                        min = 18,
                        max = 120),
           numericInput("hdl_num",
                        label = h6("HDL (mmol/L)"),
                        value = 0.97,
                        min = 0,
                        max = 6),
           numericInput("bmi_num",
                        label = h6("BMI (kg/m2)"),
                        value = 40.3,
                        min = 10,
                        max = 90)
           ),
    column(2,
           numericInput("platelets_num",
                        label = h6("Platelets"),
                        value = 225,
                        min = 0,
                        max = 1500),
           numericInput("t2dmduration_num",
                        label = h6("Time to prescription (years)"),
                        value = 8.37,
                        min = 0,
                        max = 50),
           numericInput("alb_num",
                        label = h6("Albumin (mg/g)"),
                        value = 40,
                        min = 0,
                        max = 70),
           numericInput("sys_num",
                        label = h6("Systolic Blood Pressure (mmHg)"),
                        value = 120,
                        min = 50,
                        max = 250)
           ),
    column(2,
           numericInput("ast_num",
                        label = h6("Aspartate aminotransferase (U/L)"),
                        value = 28.9,
                        min = 0,
                        max = 400),
           numericInput("bil_num",
                        label = h6("BIL"),
                        value = 8,
                        min = 0,
                        max = 135),
           numericInput("score_num",
                        label = h6("CVD score"),
                        value = 0.01,
                        min = 0,
                        max = 25)
           )
  ),
  br(), 
  
  # Conditional panel in case output_type = 1 (absolute probability)
  conditionalPanel(condition = "input.output_type == 1",
                   
                   # title of panel
                   plotOutput("absolute_plot", height = "160px", width = "800px")
  ),
  
  # Conditional panel in case output_type = 2 (outcome interval)
  conditionalPanel(condition = "input.output_type == 2",
                   
                   # title of panel
                   plotOutput("outcome_interval", height = "160px", width = "800px")
  ),
  
  
  # Conditional panel in case output_type = 3 (probability density)
  conditionalPanel(condition = "input.output_type == 3",
                   
                   # title of panel
                   plotOutput("probability_plot", height = "250px", width = "800px")
  )   
            
  
)


# Define server logic to plot various variables against mpg ----
server <- function(input, output, session) {
  
  # target HbA1c mmol/mol
  target <- eventReactive(input$calculate, {
    target <- as.numeric(input$target_num)
  })
  
  # join patient's data
  patient <- eventReactive(input$calculate, {
    # join patient's data
    patient <- NULL
    # therapy
    patient$drugclass <- factor("SGLT2", levels = c("GLP1", "SGLT2"))
    # egfr
    patient$egfr_ckdepi <- as.numeric(input$egfr_num)
    # hba1cmonth
    patient$hba1cmonth <- as.numeric(12)
    # prealt
    patient$prealt <- as.numeric(input$alt_num)
    # prehba1cmmol
    patient$prehba1cmmol <- as.numeric(input$hba1c_num)
    # score.excl.mi
    patient$score.excl.mi <- as.numeric(input$score_num)
    # Category (smoker)
    if (input$smoker_select == "Missing") {
      patient$Category <- factor(NA, levels = c("Active smoker", "Ex-smoker", "Non-smoker"))
    } else {
      patient$Category <- factor(input$smoker_select, levels = c("Active smoker", "Ex-smoker", "Non-smoker"))
    }
    # drugline
    patient$drugline <- factor("2", levels = c("2", "3", "4", "5"))
    # ncurrtx
    patient$ncurrtx <- factor("1", levels = c("0", "1", "2", "3"))
    # yrdrugstart
    patient$yrdrugstart <- as.numeric(2022)
    # agetx
    patient$agetx <- as.numeric(input$age_num)
    # malesex
    if (input$sex_select == "Male") {
      patient$malesex <- factor("1", levels = c("0", "1"))
    } else if (input$sex_select == "Female") {
      patient$malesex <- factor("0", levels = c("0", "1"))
    } else {
      patient$malesex <- factor(NA, levels = c("0", "1"))
    }
    # prehdl
    patient$prehdl <- as.numeric(input$hdl_num)
    # prebmi
    patient$prebmi <- as.numeric(input$bmi_num)
    # prebil
    patient$prebil <- as.numeric(input$bil_num)
    # preplatelets
    patient$preplatelets <- as.numeric(input$platelets_num)
    # t2dmduration
    patient$t2dmduration <- as.numeric(input$t2dmduration_num)
    # prealb
    patient$prealb <- as.numeric(input$alb_num)
    # presys
    patient$presys <- as.numeric(input$sys_num)
    # preast
    patient$preast <- as.numeric(input$ast_num)
  
    # turn into data.frame
    patient <- as.data.frame(patient)
  })

  eventReactive(input$reload, {
    session$reload()
  })
  
  posteriors_hba1c <- eventReactive(input$calculate, {
    
    patient <- patient()
    
    bart_model <- readRDS("Shiny App/bart_model.rds")
    
    patient$drugclass <- factor("SGLT2", levels = c("GLP1", "SGLT2"))
    patient <- rbind(patient,
                     patient %>%
                       mutate(drugclass = factor("GLP1", levels = c("GLP1", "SGLT2"))))
    
    posteriors_hba1c <-  bartMachine::bart_machine_get_posterior(bart_model, patient)
    
    posteriors_hba1c
  })
  
  posteriors_weight <- eventReactive(input$calculate, {
    
    patient <- patient()
    
    bart_model_weights <- readRDS("Shiny App/bart_model_weights.rds")
    
    patient <- patient %>%
      mutate(drugclass = factor("SGLT2", levels = c("GLP1", "SGLT2")))
    patient <- rbind(patient,
                     patient %>%
                       mutate(drugclass = factor("GLP1", levels = c("GLP1", "SGLT2"))))
    
    posteriors_weight <-  bartMachine::bart_machine_get_posterior(bart_model_weights, patient)
    
    posteriors_weight
  })
  
  # make output plot of probability SGLT2
  output$absolute_plot <- renderPlot({
    
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors_hba1c <- posteriors_hba1c()
    
    # SGLT2 probability of achieving target and surpassing it
    probability_SGLT2 <- length(which(posteriors_hba1c$y_hat_posterior_samples[1,] < target))/length(posteriors_hba1c$y_hat_posterior_samples[1,])
    
    # GLP1 probability of achieving target and surpassing it
    probability_GLP1 <- length(which(posteriors_hba1c$y_hat_posterior_samples[2,] < target))/length(posteriors_hba1c$y_hat_posterior_samples[2,])
    
    # reshape data
    df <- data.frame(
      SGLT2 = probability_SGLT2*100,
      GLP1 = probability_GLP1*100
    ) %>%
      gather(key, val) %>%
      mutate(
        key = factor(key, levels = c("GLP1", "SGLT2")),
        Total = 100) %>%
      mutate(label = paste0(key, ": ", val,"%"))
    
    if (probability_SGLT2 > 0.999) {
      df$label[1] <- "SGLT2: > 99.9%"
    } else if (probability_SGLT2 < 0.001) {
      df$label[1] <- "SGLT2: <0.1%"
    }
    if (probability_GLP1 > 0.999) {
      df$label[2] <- "GLP1: > 99.9%"
    } else if (probability_GLP1 < 0.001) {
      df$label[2] <- "GLP1: <0.1%"
    }

    plot_target <- ggplot(df, aes(key, val)) +
      geom_col(fill = c("#f1a340","red")) +
      geom_col(aes(y = Total), colour = "white", alpha = 0.1) +
      geom_text(
        aes(y = 5, label = label),
        hjust = 0,
        # fontface = "bold",
        colour = "black",
        size = 10) +
      ggtitle("12-month probability of achieving target HbA1c") +
      coord_flip() +
      scale_colour_manual(values=c("red","#f1a340")) +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        title = element_text(size = 13))

    posteriors_weight <- posteriors_weight()

    plot_weight <- cbind(val = posteriors_weight$y_hat, key = c("SGLT2", "GLP1")) %>%
      as.data.frame() %>%
      mutate(val = as.numeric(val)) %>%
      ggplot(aes(x = key, y = val, fill = key)) +
      geom_col() +
      geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5) +
      geom_text(
        aes(y = val + 0.1, label = paste0(key, ": ", signif(val, digits = 2), " kg")),
        hjust = 0,
        # fontface = "bold",
        colour = "black",
        size = 10) +
      scale_fill_manual(values=c("red","#f1a340")) +
      coord_flip() +
      theme_classic() +
      ggtitle("Estimated 6-month weight change") +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            title = element_text(size = 13))
    
    cowplot::plot_grid(
      
      plot_target,
      
      plot_weight,
      
      ncol = 2, nrow = 1
    )
    
  })
  
  output$outcome_interval <- renderPlot({
    
    # load predictions
    posteriors_hba1c <- posteriors_hba1c()
    
    plot_outcome_interval <- cbind(
      `5%` = apply(posteriors_hba1c$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
      `50%` = apply(posteriors_hba1c$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
      `95%` = apply(posteriors_hba1c$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
      key = c("SGLT2", "GLP1")
    ) %>%
      as.data.frame() %>%
      mutate(`5%` = as.numeric(`5%`),
             `50%` = as.numeric(`50%`),
             `95%` = as.numeric(`95%`),
             key = as.factor(key)) %>%
      ggplot(aes(x = key, y = `50%`)) +
      geom_crossbar(aes(ymin = `5%`, ymax = `95%`, fill = key, colour = key)) +
      geom_text(
        aes(y = `50%`-1.3, label = key),
        hjust = 0,
        # fontface = "bold",
        colour = "black",
        size = 10) +
      scale_fill_manual(values=c("red","#f1a340")) +
      scale_colour_manual(values=c("red","#f1a340")) +
      coord_flip() +
      theme_classic() +
      ggtitle("Estimated 12-month average HbA1c") +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            title = element_text(size = 13))
    
    posteriors_weight <- posteriors_weight()
    
    plot_weight_interval <- cbind(
      `5%` = apply(posteriors_weight$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.05))),
      `50%` = apply(posteriors_weight$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.50))),
      `95%` = apply(posteriors_weight$y_hat_posterior_samples, MARGIN = 1, function(x) quantile(c(x), probs = c(0.95))),
      key = c("SGLT2", "GLP1")
    ) %>%
      as.data.frame() %>%
      mutate(`5%` = as.numeric(`5%`),
             `50%` = as.numeric(`50%`),
             `95%` = as.numeric(`95%`),
             key = as.factor(key)) %>%
      ggplot(aes(x = key, y = `50%`)) +
      geom_crossbar(aes(ymin = `5%`, ymax = `95%`, fill = key, colour = key)) +
      geom_text(
        aes(y = `50%`-1, label = key),
        hjust = 0,
        # fontface = "bold",
        colour = "black",
        size = 10) +
      scale_fill_manual(values=c("red","#f1a340")) +
      scale_colour_manual(values=c("red","#f1a340")) +
      coord_flip() +
      theme_classic() +
      ggtitle("Estimated 6-month weight change") +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            title = element_text(size = 13))
    
    cowplot::plot_grid(
      
      plot_outcome_interval,
      
      plot_weight_interval,
      
      ncol = 2, nrow = 1
    )
    
  })
  
  output$probability_plot <- renderPlot({
    
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors_hba1c <- posteriors_hba1c()
    
    plot_outcome <- posteriors_hba1c$y_hat_posterior_samples %>%
      t() %>%
      as.data.frame() %>%
      set_names("SGLT2", "GLP1") %>%
      gather(key, value) %>%
      ggplot() +
      geom_density(aes(x = value, fill = key), alpha = 0.5, colour = "black") +
      labs(x = "Average HbA1c (mmol/mol)", y = "Density") +
      scale_fill_manual(values=c("red","#f1a340")) +
      ggtitle("Estimated 12-month average HbA1c") +
      theme_classic() +
      theme(legend.position = c(0.80, 0.87),
            legend.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            title = element_text(size = 13))
      
    posteriors_weight <- posteriors_weight()
    
    plot_weight_hist <- posteriors_weight$y_hat_posterior_samples %>%
      t() %>%
      as.data.frame() %>%
      set_names("SGLT2", "GLP1") %>%
      gather(key, value) %>%
      ggplot() +
      geom_density(aes(x = value, fill = key), alpha = 0.5, colour = "black") +
      labs(x = "Average weight reduction (kg)", y = "Density") +
      scale_fill_manual(values=c("red","#f1a340")) +
      ggtitle("Estimated 6-month weight change") +
      theme_classic() +
      theme(legend.title = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            title = element_text(size = 13),
            legend.position = "none")
    
    cowplot::plot_grid(
      
      plot_outcome,
      
      plot_weight_hist,
      
      ncol = 2, nrow = 1
    )
    
    
  })
  
}

# shinyApp(ui, server)
runGadget(ui, server, viewer = browserViewer(browser = getOption("browser")))


# ## Code for testing purposes
# 
# # Used in slade to ensure the library being used is my personal library
# .libPaths(.libPaths()[c(2,1,3)])
# 
# 
# ## increase memery usage to 50gb of RAM
# options(java.parameters = "-Xmx50g")
# 
# library(tidyverse)
# library(bartMachine)
# 
# bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/cvd_new/bart_model_final.rds")
# 
# 
# patient <- cbind(
#   drugclass = "SGLT2",
#   egfr_ckdepi = 96,
#   hba1cmonth = 12,
#   prealt = 23,
#   prehba1cmmol = 65,
#   score.excl.mi = 0.01,
#   Category = "Non-smoker",
#   drugline = "2",
#   ncurrtx = "1",
#   yrdrugstart = 2017,
#   agetx = 65,
#   malesex = "1",
#   prehdl = 0.97,
#   prebmi = 40.3,
#   prebil = 8,
#   preplatelets = 225,
#   t2dmduration = 8.4,
#   prealb = 40,
#   presys = 120,
#   preast = NA
# ) %>%
#   as.data.frame() %>%
#   mutate(drugclass = factor(drugclass, levels = c("GLP1", "SGLT2")),
#          egfr_ckdepi = as.numeric(egfr_ckdepi),
#          hba1cmonth = as.numeric(hba1cmonth),
#          prealt = as.numeric(prealt),
#          prehba1cmmol = as.numeric(prehba1cmmol),
#          score.excl.mi = as.numeric(score.excl.mi),
#          Category = factor(Category, levels = c("Active smoker", "Ex-smoker", "Non-smoker")),
#          drugline = factor(drugline, levels = c("2", "3", "4", "5")),
#          ncurrtx = factor(ncurrtx, levels = c("0", "1", "2", "3")),
#          yrdrugstart = as.numeric(yrdrugstart),
#          agetx = as.numeric(agetx),
#          malesex = factor(malesex, levels = c("0", "1")),
#          prehdl = as.numeric(prehdl),
#          prebmi = as.numeric(prebmi),
#          prebil = as.numeric(prebil),
#          preplatelets = as.numeric(preplatelets),
#          t2dmduration = as.numeric(t2dmduration),
#          prealb = as.numeric(prealb),
#          presys = as.numeric(presys),
#          preast = as.numeric(preast)
#   )
# 
# test <- bart_model_final
# 
# library(synthpop)
# 
# syn.dataset <- syn(cbind(test$X, y = test$y))
# test$X <- syn.dataset$syn[1:100,-21]
# test$y <- syn.dataset$syn[1:100,"y"]
# 
# # results were still the same because they are deterministic
# posteriors_SGLT2 <- bartMachine::bart_machine_get_posterior(test, patient)
# summary(t(posteriors_SGLT2$y_hat_posterior_samples))
# 
# saveRDS(test, "bart_model.rds")
# 
# # model for weight reduction
# 
# bart_model_weights <- readRDS("Samples/SGLT2-GLP1/Weight_reduction/bart_model_weight.rds")
# 
# test <- bart_model_weights
# 
# library(synthpop)
# 
# syn.dataset <- syn(cbind(test$X, y = test$y))
# test$X <- syn.dataset$syn[1:100,-21]
# test$y <- syn.dataset$syn[1:100,"y"]
# 
# # posteriors
# posteriors_weights <- bartMachine::bart_machine_get_posterior(test, patient)
# summary(t(posteriors_weights$y_hat_posterior_samples))
# 
# saveRDS(test, "bart_model_weights.rds")




####################
## Description:
##  - In this file we make a Shiny App for predicting the probability of 
##      achieving a target HbA1c.
####################

library(shiny)

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
           actionButton("calculate", "Calculate")
           
           ),
    
    column(2,
           selectInput("sex_select",
                       label = h6("Sex"),
                       choices = list("",
                                      "Male",
                                      "Female"),
                       selected = NULL),
           numericInput("egfr_num",
                        label = h6("eGFR"),
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
                                      "Active smoker"),
                       selected = NULL),
           numericInput("age_num",
                        label = h6("Age"),
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
                        label = h6("Time to prescription"),
                        value = 8.37,
                        min = 0,
                        max = 50),
           numericInput("alb_num",
                        label = h6("ALB"),
                        value = 40,
                        min = 0,
                        max = 70),
           numericInput("sys_num",
                        label = h6("Systolic Blood Pressure"),
                        value = 120,
                        min = 50,
                        max = 250)
           ),
    column(2,
           numericInput("ast_num",
                        label = h6("Aspartate aminotransferase"),
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
                        max = 1)
           )
  ),
  br(), 
  
  
  # title of panel
  h2("Chance of achieving HbA1c target at 12 months"),
  
  h3(textOutput("probability_SGLT2_text")),
  plotOutput("probability_SGLT2_plot", height = "100px", width = "600px"),
  h3(textOutput("probability_GLP1_text")),
  plotOutput("probability_GLP1_plot", height = "100px", width = "600px")
  
)


# Define server logic to plot various variables against mpg ----
server <- function(input, output, session) {
  
  
  ## increase memory usage to 3gb of RAM
  options(java.parameters = "-Xmx3g")
  
  require(bartMachine)
  require(tidyverse)
  
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
    # score
    patient$score <- as.numeric(input$score_num)
    # Category (smoker)
    patient$Category <- factor(input$smoker_select, levels = c("Active smoker", "Ex-smoker", "Non-smoker"))
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
    } else {
      patient$malesex <- factor("0", levels = c("0", "1"))
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
  
  posteriors <- eventReactive(input$calculate, {
    
    patient <- patient()
    
    bart_model <- readRDS("Shiny App/bart_model.rds")
    
    patient$drugclass <- factor("SGLT2", levels = c("GLP1", "SGLT2"))
    patient <- rbind(patient,
                     patient %>%
                       mutate(drugclass = factor("GLP1", levels = c("GLP1", "SGLT2"))))
    
    posteriors <-  bartMachine::bart_machine_get_posterior(bart_model, patient)
    
  })
  # make output text of probability SGLT2
  output$probability_SGLT2_text <- renderText({
  
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors <- posteriors()
    
    # SGLT2 probability of achieving target and surpassing it
    probability_SGLT2 <- length(which(posteriors$y_hat_posterior_samples[1,] < target))/length(posteriors$y_hat_posterior_samples[1,])
    
    text <- paste0("SGLT2-inhibitors: ")
    
  })
  
  # make output plot of probability SGLT2
  output$probability_SGLT2_plot <- renderPlot({
    
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors <- posteriors()
    
    # SGLT2 probability of achieving target and surpassing it
    probability_SGLT2 <- length(which(posteriors$y_hat_posterior_samples[1,] < target))/length(posteriors$y_hat_posterior_samples[1,])
    
    # reshape data
    df <- data.frame(
      SGLT2 = probability_SGLT2*100
    ) %>%
      gather(key, val) %>%
      mutate(
        key = factor(key, rev(unique(key))),
        Total = 100)

    plot <- ggplot(df, aes(key, val)) +
      geom_col(fill = "forestgreen") +
      geom_col(aes(y = Total), alpha = 0.5, colour = "black") +
      geom_text(
        aes(y = 5, label = paste0(val,"%")),
        hjust = 0,
        fontface = "bold",
        colour = "white",
        size = 10) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

    plot
  })

  # make output text of probability GLP1
  output$probability_GLP1_text <- renderText({
    
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors <- posteriors()
    
    # SGLT2 probability of achieving target and surpassing it
    probability_GLP1 <- length(which(posteriors$y_hat_posterior_samples[2,] < target))/length(posteriors$y_hat_posterior_samples[2,])
    
    text <- paste0("GLP1-agnostics: ")

  })
  
  # make output plot of probability GLP1
  output$probability_GLP1_plot <- renderPlot({
    
    # load target HbA1c 
    target <- target()
    
    # load predictions
    posteriors <- posteriors()
    
    # GLP1 probability of achieving target and surpassing it
    probability_GLP1 <- length(which(posteriors$y_hat_posterior_samples[2,] < target))/length(posteriors$y_hat_posterior_samples[2,])
    
    # reshape data
    df <- data.frame(
      GLP1 = probability_GLP1*100
    ) %>%
      gather(key, val) %>%
      mutate(
        key = factor(key, rev(unique(key))),
        Total = 100)
    
    plot <- ggplot(df, aes(key, val)) +
      geom_col(fill = "forestgreen") +
      geom_col(aes(y = Total), alpha = 0.5, colour = "black") +
      geom_text(
        aes(y = 5, label = paste0(val, "%")),
        hjust = 0,
        fontface = "bold",
        colour = "white",
        size = 10) +
      coord_flip() +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
    
    plot
    
    
  })
  
}


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
# bart_model_final <- readRDS("Samples/SGLT2-GLP1/Final_model/With_grf_no_prop/bart_model_final.rds")
# 
# 
# patient <- cbind(
#   drugclass = "SGLT2",
#   egfr_ckdepi = 96,
#   hba1cmonth = 12,
#   prealt = 23,
#   prehba1cmmol = 65,
#   score = 0.01,
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
#          score = as.numeric(score),
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

####################
## Description:
##  - In this file we make a Shiny App for predicting the probability of 
##      achieving a target HbA1c.
####################

library(shiny)

#:----------------------------------------------------------------------
### UI


# Define UI for inputting patient values
ui <- fluidPage(
  
  ### colour for the rows
  tags$head(tags$style(HTML(".grey-row {background-color: #f2f2f2;}"))),
  ### colour for titles
  tags$head(tags$style(HTML(".title-blue {color: #0000CC;}"))),
  ### size of title for hba1c plot
  tags$head(tags$style("#hba1c_outcome_error_bar_text{font-size: 20px;}")),
  tags$head(tags$style("#hba1c_outcome_error_bar_title{font-size: 23px;}")),
  ### size of title for additional outcomes
  tags$head(tags$style("#weight_outcome_text{font-size: 20px;}")),
  tags$head(tags$style("#weight_outcome_title{font-size: 23px;}")),
  tags$head(tags$style("#discontinuation_outcome_text{font-size: 20px;}")),
  tags$head(tags$style("#discontinuation_outcome_title{font-size: 23px;}")),
  ### inputs height
  
  fluidRow(
    
    
    column(10, 
           # Title of the page
           h2("SGLT2-inhibitor and GLP-1 receptor agonist Research Tool Decision Aid", 
              style="font-weight:bold;"),
           # PSA for using the app
           h5("Please note this is a beta version provided for academic research and validation purposes and should not be used for clinical decision making.", 
              style="font-weight:bold;"),
           # Info about what the app does
           h5("The model uses an individual's clinical information to provide individualised estimates of likely achieved blood glucose control (HbA1c) benefit on SGLT2-inhibitor or GLP-1 receptor agonist therapy."),
           
    ),
    
    column(2,
           br(),
           img(src='uni_logo.JPG', align = "right", width = "200"),
    )
    
  ),
  
  
  column(4,
         
         fluidRow(
           class = "grey-row",
           column(12,
                  h4(class = "title-blue", "Clinical information:")
           )
         ),
         
         fluidRow(
           class = "grey-row",
           column(6,
                  h5("Biomarkers:", style="font-weight:bold;"),
                  selectInput("sex_select",
                              label = h6("Sex"),
                              choices = list("Male",
                                             "Female"),
                              selected = "Female"),
                  numericInput("age_num",
                               label = h6("Age ( years )"),
                               value = 65,
                               min = 18,
                               max = 120),
                  numericInput("t2dmduration_num",
                               label = h6("T2D duration ( years )"),
                               value = 8,
                               min = 0,
                               max = 102),
                  numericInput("bmi_num",
                               label = h6("BMI ( kg / m", tags$sup("2"), ")"),
                               value = 40,
                               min = 15,
                               max = 100)
           ),
           column(6,
                  br(),
                  br(),
                  numericInput("hba1c_num",
                               label = h6("Baseline HbA1c ( mmol / mol )"),
                               value = 65,
                               min = 53,
                               max = 120),
                  numericInput("alt_num",
                               label = h6("ALT ( U / L )"),
                               value = 35,
                               min = 0,
                               max = 200),
                  numericInput("egfr_num",
                               label = h6("eGFR ( ml / min / 1.73 m", tags$sup("2"), ") *"),
                               value = 90,
                               min = 45,
                               max = 300),
                  numericInput("creatinine_num",
                               label = h6("Serum creatinine ( \u03BCmol / L ) [optional] *"),
                               value = NA,
                               min = 0,
                               max = 400)
                  
           )
           
         ),
         
         fluidRow(
           class = "grey-row",
           column(6,
                  h5("Cardiovascular conditions:", style="font-weight:bold;"),
                  selectInput("prepad_select",
                              label = h6("Peripheral arterial disease"),
                              choices = list("No",
                                             "Yes"),
                              selected = "No"),
                  selectInput("preheartfailure_select",
                              label = h6("Heart failure"),
                              choices = list("No",
                                             "Yes"),
                              selected = "No"),
                  selectInput("preihd_select",
                              label = h6("Ischaemic heart disease"),
                              choices = list("No",
                                             "Yes"),
                              selected = "No")
           ),
           column(6,
                  h5("Microvascular complications:", style="font-weight:bold;"),
                  selectInput("preneuropathy_select",
                              label = h6("Neuropathy"),
                              choices = list("No",
                                             "Yes"),
                              selected = "No"),
                  selectInput("preretinopathy_select",
                              label = h6("Retinopathy"),
                              choices = list("No",
                                             "Yes"),
                              selected = "No")
                  
                  
           )
           
           
         ),
         
         fluidRow(
           class = "grey-row",
           column(4,
                  actionButton("calculate", "Calculate"),
                  br(),
                  br()
           )
         ),
         
         fluidRow(
           class = "grey-row",
           column(12,
                  h5("Notes:"),
                  "* Either eGFR or serum creatinine can be provided. eGFR will be calculated from serum creatinine (sex and age is required) if a serum creatinine is entered. If the derived eGFR is less than 45, model outputs will be generated fixing eGFR at 45.",
                  tags$a(href="https://www.kidney.org/professionals/kdoqi/gfr_calculator/formula", 
                         "Find out more here!", target="_blank"),
                  br(),
                  br(),
                  "For any enquiries, please contact pml204@exeter.ac.uk",
                  br(),
                  br()
                  
                  
           )
         )
         
         
  ),
  
  column(8,
         
         # Conditional panel in case of error (value outside range)
         conditionalPanel(
           condition = "input.age_num < 18 || input.age_num > 120 ||
                      input.bmi_num < 15 || input.bmi_num > 100 ||
                      input.hba1c_num < 53 || input.hba1c_num > 120 ||
                      input.egfr_num < 45 || input.egfr_num > 300 ||
                      input.creatinine_num < 0 || input.creatinine_num > 400 ||
                      input.t2dmduration_num < 0 || input.t2dmduration_num > 102 ||
                      input.alt_num < 0 || input.alt_num > 200",
           h4("Error:", style="font-weight:bold;")
         ),
         conditionalPanel(
           condition = "input.t2dmduration_num < 0 || input.t2dmduration_num > 102",
           h4("T2D duration value outside range. Please enter a value between 0 and 102.")
         ),
         conditionalPanel(
           condition = "input.age_num < 18 || input.age_num > 120",
           h4("Age value outside range. Please enter a value between 18 and 120.")
         ),
         conditionalPanel(
           condition = "input.bmi_num < 15 || input.bmi_num > 100",
           h4("BMI value outside range. Please enter a value between 15 and 100.")
         ),
         conditionalPanel(
           condition = "input.hba1c_num < 53 || input.hba1c_num > 120",
           h4("Baseline HbA1c value outside range. Please enter a value between 53 and 120.")
         ),
         conditionalPanel(
           condition = "input.alt_num < 0 || input.alt_num > 200",
           h4("ALT value outside range. Please enter a value between 0 and 200.")
         ),
         conditionalPanel(
           condition = "input.egfr_num < 45 || input.egfr_num > 300",
           h4("eGFR value outside range. Please enter a value between 45 and 300.")
         ),
         conditionalPanel(
           condition = "input.creatinine_num < 0 || input.creatinine_num > 400",
           h4("Serum creatinine value outside range. Please enter a value between 0 and 400.")
         ),
         
         fluidRow(
           conditionalPanel(
             condition = "input.calculate > 0 &
          !(input.age_num < 18 || input.age_num > 120 ||
          input.bmi_num < 15 || input.bmi_num > 100 ||
          input.hba1c_num < 53 || input.hba1c_num > 120 ||
          input.egfr_num < 45 || input.egfr_num > 300 ||
          input.creatinine_num < 0 || input.creatinine_num > 400 ||
          input.t2dmduration_num < 0 || input.t2dmduration_num > 102 ||
          input.alt_num < 0 || input.alt_num > 200)",
             
             
             # HbA1c panel
             column(12,
                    ### HbA1c outcome error bars
                    strong(textOutput("hba1c_outcome_error_bar_title")),
                    # title
                    strong(textOutput("hba1c_outcome_error_bar_text")),
                    # plots
                    plotOutput("hba1c_outcome_error_bar", height = "200px"),
                    ### Weight outcome
                    strong(textOutput("weight_outcome_title")),
                    # title
                    strong(textOutput("weight_outcome_text")),
                    # plots
                    plotOutput("weight_outcome", height = "200px"),
                    ### Discontinuation outcome
                    strong(textOutput("discontinuation_outcome_title")),
                    # title
                    strong(textOutput("discontinuation_outcome_text")),
                    # plots
                    plotOutput("discontinuation_outcome", height = "200px")
             )
             
             
           )
         )
         
  )
  
)


#:--------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  require(plyr)
  require(ggtext)
  require(tidyverse)
  require(bcf)
  
  # join patient's data
  patient <- eventReactive(input$calculate, {
    # join patient's data
    patient <- NULL
    # therapy
    patient$drugclass <- factor("SGLT2", levels = c("GLP1", "SGLT2"))
    # sex
    patient$sex <- factor(input$sex_select, levels = c("Female", "Male"))
    # agetx
    if (is.na(input$age_num)) {patient$agetx <- as.numeric(65)} else {
      patient$agetx <- as.numeric(input$age_num)
    }
    # prebmi
    if (is.na(input$bmi_num)) {patient$prebmi <- as.numeric(40)} else {
      patient$prebmi <- as.numeric(input$bmi_num)
    }
    # prehba1cmmol
    if (is.na(input$hba1c_num)) {patient$prehba1c <- as.numeric(65)} else {
      patient$prehba1c <- as.numeric(input$hba1c_num)
    }
    # egfr (calculated from creatinine or use given egfr)
    if (!is.na(input$creatinine_num)) {
      creatinine_mg_dL <- input$creatinine_num*0.0113
      egfr_value <- as.numeric(ifelse(creatinine_mg_dL<=0.7 & !!patient$sex=="Female",(142 * ((creatinine_mg_dL/0.7)^-0.241) * (0.9938^!!patient$agetx) * 1.012),
                                      ifelse(creatinine_mg_dL>0.7 & !!patient$sex=="Female",(142 * ((creatinine_mg_dL/0.7)^-1.2) * (0.9938^!!patient$agetx) * 1.012),
                                             ifelse(creatinine_mg_dL<=0.9 & !!patient$sex=="Male",(142 * ((creatinine_mg_dL/0.9)^-0.302) * (0.9938^!!patient$agetx)),
                                                    ifelse(creatinine_mg_dL>0.9 & !!patient$sex=="Male",(142 * ((creatinine_mg_dL/0.9)^-1.2) * (0.9938^!!patient$agetx)),NA)))))
      patient$preegfr <- as.numeric(ifelse(egfr_value < 45, 45, ifelse(egfr_value > 300, 300, egfr_value)))
    } else {
      if (is.na(input$egfr_num)) {patient$preegfr <- as.numeric(90)} else {
        patient$preegfr <- as.numeric(input$egfr_num)
      }
    }
    # t2dmduration
    if (is.na(input$t2dmduration_num)) {patient$t2dmduration <- as.numeric(8)} else {
      patient$t2dmduration <- as.numeric(input$t2dmduration_num)
    }
    # alt
    if (is.na(input$alt_num)) {patient$prealt <- as.numeric(35)} else {
      patient$prealt <- as.numeric(input$alt_num)
    }
    # Peripheral arterial disease
    patient$prepad <- factor(input$prepad_select, levels = c("No", "Yes"))
    # Heart failure
    patient$preheartfailure <- factor(input$preheartfailure_select, levels = c("No", "Yes"))
    # Ischaemic heart disease
    patient$preihd <- factor(input$preihd_select, levels = c("No", "Yes"))
    # Neuropathy
    patient$preneuropathy <- factor(input$preneuropathy_select, levels = c("No", "Yes"))
    # Retinopathy
    patient$preretinopathy <- factor(input$preretinopathy_select, levels = c("No", "Yes"))
    # ncurrtx
    patient$ncurrtx <- factor("1", levels = c("1", "2", "3", "4", "5+"))
    # drugline
    patient$drugline <- factor("2", levels = c("2", "3", "4", "5+"))
    # hba1cmonth
    patient$hba1cmonth <- as.numeric(12)
    # turn into data.frame
    patient <- as.data.frame(patient)
  })
  
  ## if creatinine is provided, eGFR is recalculated from creatinine and input is updated.
  observeEvent(input$calculate, {
    if(!is.na(input$creatinine_num)) {
      
      patient <- patient()
      
      updateNumericInput(
        session = session,
        inputId = "egfr_num",
        value = round(patient$preegfr, digits = 1)
      ) 
    }
  })
  
  
  # Posteriors treatment effect
  posterior_tau_mu <- eventReactive(input$calculate, {
    
    # Show a message to the user that their data is being processed
    showModal(modalDialog("Your data is being processed. Please wait...",
                          footer = NULL,
                          tags$head(
                            tags$style(HTML("
                              .modal {
                                position: fixed;
                                top: 50%;
                                left: 50%;
                                transform: translate(-50%, -50%);
                              }
                            "))
                          )))
    
    # load patient
    patient <- patient() %>%
      rbind(patient() %>%
              mutate(drugclass = factor("GLP1", levels = c("GLP1", "SGLT2"))))
    
    bcf_model <- readRDS("bcf_model.rds")
    
    predictions <-  predict(object = bcf_model,
                            x_predict_control = patient %>%
                              select(
                                c("agetx", "t2dmduration", "drugline", "ncurrtx", "hba1cmonth", "prehba1c", "preegfr", "prealt", "prepad")
                              ) %>%
                              mutate_all(list(~as.numeric(.))) %>%
                              as.matrix(),
                            x_predict_moderate = patient %>%
                              select(
                                c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preheartfailure" ,"preihd", "preneuropathy", "prepad" ,"preretinopathy")
                              ) %>%
                              mutate_all(list(~as.numeric(.))) %>%
                              as.matrix(),
                            pi_pred = c(0.5, 0.5),
                            z_pred = patient %>%
                              select(drugclass) %>%
                              mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                              unlist(),
                            save_tree_directory = ".", 
                            log_file = NULL)
    
    
    return(predictions)
    
  })
  
  
  # hide the message when the calculation is done
  observeEvent(posterior_tau_mu(), {
    removeModal()
  })
  
  
  
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  #:------------------------------  HbA1c outcome  ----------------------------------
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  
  
  
  output$hba1c_outcome_error_bar <- renderPlot({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame() %>%
      mutate(
        low_qt = as.numeric(
          c(
            quantile(predictions$yhat[,1], 0.025),
            quantile(predictions$yhat[,2], 0.025)
          )),
        high_qt = as.numeric(
          c(
            quantile(predictions$yhat[,1], 0.975),
            quantile(predictions$yhat[,2], 0.975)
          ))
      )
    
    
    
    predictions_hba1c %>%
      ggplot(aes(x = drugclass, y = value, colour = drugclass, ymin = low_qt, ymax = high_qt)) +
      geom_hline(yintercept = patient$prehba1c, colour = "black", linewidth = 2, alpha = 1) +
      geom_point(size = 10) +
      geom_errorbar(linewidth = 2, width = 0.7) +
      geom_richtext(aes(x = 1.5, y = patient$prehba1c, label = "Baseline HbA1c"), size = 5, colour = "black", fill = "white", angle = 0, hjust = "inward") +
      scale_colour_manual(values = c("#f1a340", "dodgerblue2")) +scale_x_discrete(labels = c("SGLT2i" = paste0("SGLT2i:\n", round(predictions_hba1c$value[1], digits = 1), " mmol/mol"),
                                                                                             "GLP1-RA" = paste0("GLP1-RA:\n", round(predictions_hba1c$value[2], digits = 1), " mmol/mol"))) +
      theme_classic() +
      labs(
        y = "Outcome HbA1c (mmol/mol)"
      ) +
      theme(legend.position = "none",
            title = element_text(size = 30),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.ticks.y = element_blank()) +
      coord_flip()
    
    
  })
  
  
  output$hba1c_outcome_error_bar_text <- renderText({
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
  
    paste0("Additional average HbA1c reduction from baseline on ", as.character(predictions_hba1c[which.min(predictions_hba1c$value), "drugclass"]), ": ", round(abs(round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)), digits = 1), " mmol/mol")
    
  })
  
  output$hba1c_outcome_error_bar_title <- renderText({
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    paste0("12-month HbA1c outcome:")
    
  })
  
  
  
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  #:------------------------------  Weight outcome  ---------------------------------
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  
  
  
  output$weight_outcome <- renderPlot({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    weight_outcomes <- readRDS("weight_outcome.rds") %>%
      filter(sex == patient$sex)
    
    if (treatment_effect <= -5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-20,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(5,31]") %>% select(-sex, -intervals)
    }
    
    weight_outcomes %>%
      ggplot(aes(x = drugclass, y = mean, colour = drugclass, ymin = lci, ymax = uci)) +
      geom_point(size = 10) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 2, alpha = 1) +
      geom_richtext(aes(x = 1.5, y = 0, label = "No weight change"), size = 5, colour = "black", fill = "white", angle = 0, hjust = "inward") +
      geom_errorbar(linewidth = 2, width = 0.7) +
      scale_colour_manual(values = c("#f1a340", "dodgerblue2")) +
      scale_x_discrete(labels = c("SGLT2i" = paste0("SGLT2i:\n", round(weight_outcomes$mean[1], digits = 1), " kg"),
                                  "GLP1-RA" = paste0("GLP1-RA:\n", round(weight_outcomes$mean[2], digits = 1), " kg"))) +
      theme_classic() +
      labs(
        y = "Average weight change (kg)"
      ) +
      theme(legend.position = "none",
            title = element_text(size = 30),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.ticks.y = element_blank()) +
      coord_flip()
    
  })
  
  
  output$weight_outcome_text <- renderText({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    weight_outcomes <- readRDS("weight_outcome.rds") %>%
      filter(sex == patient$sex)
    
    if (treatment_effect <= -5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-20,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(5,31]") %>% select(-sex, -intervals)
    }
    
    paste0("Additional average weight change from ", as.character(weight_outcomes[which.min(weight_outcomes$mean), "drugclass"]), ": ", round(abs(round(weight_outcomes$mean[1], digits = 1) - round(weight_outcomes$mean[2], digits = 1)), digits = 1), " kg")
    
  })
  
  output$weight_outcome_title <- renderText({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    weight_outcomes <- readRDS("weight_outcome.rds") %>%
      filter(sex == patient$sex)
    
    if (treatment_effect <= -5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-20,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      weight_outcomes <- weight_outcomes %>% filter(intervals == "(5,31]") %>% select(-sex, -intervals)
    }
    
    paste0("12-month weight outcome:")
    
  })
  
  
  
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  #:-------------------------  Discontinuation outcome  -----------------------------
  #:---------------------------------------------------------------------------------
  #:---------------------------------------------------------------------------------
  
  output$discontinuation_outcome <- renderPlot({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    discontinuation_outcomes <- readRDS("discontinuation_outcome.rds") %>%
      filter(sex == patient$sex) %>%
      mutate(mean = mean*100,
             lci = lci*100,
             uci = uci*100)
    
    if (treatment_effect <= -5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-21,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(5,40]") %>% select(-sex, -intervals)
    }
    
    discontinuation_outcomes %>%
      ggplot(aes(x = drugclass, y = mean, colour = drugclass, ymin = lci, ymax = uci)) +
      geom_point(size = 10) +
      geom_errorbar(linewidth = 2, width = 0.7) +
      scale_colour_manual(values = c("#f1a340", "dodgerblue2")) +
      scale_x_discrete(labels = c("SGLT2i" = paste0("SGLT2i:\n", round(discontinuation_outcomes$mean[1], digits = 1), "%"),
                                  "GLP1-RA" = paste0("GLP1-RA:\n", round(discontinuation_outcomes$mean[2], digits = 1), "%"))) +
      theme_classic() +
      ylim(0, max(discontinuation_outcomes[,-4])) +
      labs(
        y = "Discontinuation (%)"
      ) +
      theme(legend.position = "none",
            title = element_text(size = 30),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.ticks.y = element_blank()) +
      coord_flip()
    
    
  })
  
  
  output$discontinuation_outcome_text <- renderText({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    discontinuation_outcomes <- readRDS("discontinuation_outcome.rds") %>%
      filter(sex == patient$sex) %>%
      mutate(mean = mean*100,
             lci = lci*100,
             uci = uci*100)
    
    if (treatment_effect <= -5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-21,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(5,40]") %>% select(-sex, -intervals)
    }
    
    paste0("Reduced average discontinuation from ", as.character(discontinuation_outcomes[which.min(discontinuation_outcomes$mean), "drugclass"]), ": ", round(abs(round(discontinuation_outcomes$mean[1], digits = 1) - round(discontinuation_outcomes$mean[2], digits = 1)), digits = 1), "%")
    
  })

  
  output$discontinuation_outcome_title <- renderText({
    
    patient <- patient()
    
    predictions <- posterior_tau_mu()
    
    predictions_hba1c <- predictions$yhat %>%
      colMeans() %>%
      as.data.frame() %>%
      mutate(drugclass = factor(c("SGLT2i", "GLP1-RA"), levels = c("SGLT2i", "GLP1-RA"))) %>%
      set_names(c("value", "drugclass")) %>%
      as.data.frame()
    
    treatment_effect <- round(predictions_hba1c$value[1], digits = 1) - round(predictions_hba1c$value[2], digits = 1)
    
    discontinuation_outcomes <- readRDS("discontinuation_outcome.rds") %>%
      filter(sex == patient$sex) %>%
      mutate(mean = mean*100,
             lci = lci*100,
             uci = uci*100)
    
    if (treatment_effect <= -5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-21,-5]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= -3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-5,-3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 0 ) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(-3,0]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 3) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(0,3]") %>% select(-sex, -intervals)
    } else if (treatment_effect <= 5) {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(3,5]") %>% select(-sex, -intervals)
    } else {
      discontinuation_outcomes <- discontinuation_outcomes %>% filter(intervals == "(5,40]") %>% select(-sex, -intervals)
    }
    
    paste0("6-month discontinuation outcome:")
    
  })
  
  
}



shinyApp(ui, server)
# runGadget(ui, server, viewer = browserViewer(browser = getOption("browser")))

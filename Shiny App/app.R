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
  tags$head(tags$style(HTML(".title-blue { color: #0000CC; }"))),
  
  # Title of the page
  h1("SGLT2-inhibitor and GLP-1 receptor agonist Treatment Selection Decision Aid", 
     style="font-weight:bold;"),
  
  # PSA for using the app
  h4("Prease note this is a beta version provided for academic research and vlaidation purposes and should not be used for clinical decision making.", 
     style="font-weight:bold;"),
  
  br(),
  
  # Info about what the app does
  h4("The model uses an individual's clinical information to provide individualised estimates of likely achieved blood glucose control (HbA1c) benefit on SGLT2-inhibitor or GLP-1 receptor agonist therapy."),
  
  # App title
  titlePanel("HbA1c prediction"),
  
  fluidRow(
    class = "grey-row",
    column(2,
           h3(class = "title-blue", "Clinical information:"),
           
           br(),
           
           actionButton("calculate", "Calculate"),
           
           br(),
           
           radioButtons("output_type",
                        label = h3("Output type:"),
                        choices = list("Absolute probability" = 1,
                                       "Histogram" = 2),
                        selected = 1)
    ),
    
    column(2,
           h4("Biomarkers:"),
           selectInput("sex_select",
                       label = h6("Sex"),
                       choices = list("Male",
                                      "Female"),
                       selected = "Male"),
           numericInput("age_num",
                        label = h6("Age (years)"),
                        value = 65,
                        min = 18,
                        max = 120),
           numericInput("t2dmduration_num",
                        label = h6("Type 2 diabetes duration (years)"),
                        value = 8.37,
                        min = 0,
                        max = 102),
           numericInput("bmi_num",
                        label = h6("BMI (kg/m", tags$sup("2"), ")"),
                        value = 40.3,
                        min = 15,
                        max = 100)
    ),
    column(2,
           br(),
           br(),
           numericInput("hba1c_num",
                        label = h6("Baseline HbA1c (mmol/mol)"),
                        value = 65,
                        min = 25,
                        max = 120),
           numericInput("egfr_num",
                        label = h6("eGFR (ml/min/1.73m", tags$sup("2"), ")"),
                        value = 96,
                        min = 45,
                        max = 300),
           numericInput("alt_num",
                        label = h6("ALT (U/L)"),
                        value = 23,
                        min = 0,
                        max = 200)
    ),
    column(2,
           h4("Cardiovascular conditions:"),
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
    column(2,
           h4("Microvascular complications:"),
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
  
  # Conditional panel in case of error (value outside range)
  conditionalPanel(
    condition = "input.age_num < 18 || input.age_num > 120 || 
                   input.t2dmduration_num < 0 || input.t2dmduration_num > 102 ||
                   input.bmi_num < 15 || input.bmi_num > 100 ||
                   input.hba1c_num < 25 || input.hba1c_num > 120 ||
                   input.egfr_num < 45 || input.egfr_num > 300 ||
                   input.alt_num < 0 || input.alt_num > 200",
    h4("Error:", style="font-weight:bold;")
  ),
  conditionalPanel(
    condition = "input.age_num < 18 || input.age_num > 120",
    h4("Age value outside range. Please enter a value between 18 and 120.")
  ),
  conditionalPanel(
    condition = "input.t2dmduration_num < 0 || input.t2dmduration_num > 102",
    h4("Type 2 diabetes duration value outside range. Please enter a value between 0 and 102.")
  ),
  conditionalPanel(
    condition = "input.bmi_num < 15 || input.bmi_num > 100",
    h4("BMI value outside range. Please enter a value between 15 and 100.")
  ),
  conditionalPanel(
    condition = "input.hba1c_num < 25 || input.hba1c_num > 120",
    h4("Baseline HbA1c value outside range. Please enter a value between 25 and 120.")
  ),
  conditionalPanel(
    condition = "input.egfr_num < 45 || input.egfr_num > 300",
    h4("eGFR value outside range. Please enter a value between 45 and 300.")
  ),
  conditionalPanel(
    condition = "input.alt_num < 0 || input.alt_num > 200",
    h4("ALT value outside range. Please enter a value between 0 and 200.")
  ),
  
  # Conditional panel in case output_type = 1 (absolute probability)
  conditionalPanel(condition = "input.output_type == 1 & input.calculate > 0 & 
                   !(input.age_num < 18 || input.age_num > 120 || 
                   input.t2dmduration_num < 0 || input.t2dmduration_num > 102 ||
                   input.bmi_num < 15 || input.bmi_num > 100 ||
                   input.hba1c_num < 25 || input.hba1c_num > 120 ||
                   input.egfr_num < 45 || input.egfr_num > 300 ||
                   input.alt_num < 0 || input.alt_num > 200)",
                   
                   # title of panel
                   plotOutput("effect_probability", height = "250px", width = "800px")
  ),
  
  # Conditional panel in case output_type = 2 (outcome interval)
  conditionalPanel(condition = "input.output_type == 2 & input.calculate > 0 & 
                   !(input.age_num < 18 || input.age_num > 120 || 
                   input.t2dmduration_num < 0 || input.t2dmduration_num > 102 ||
                   input.bmi_num < 15 || input.bmi_num > 100 ||
                   input.hba1c_num < 25 || input.hba1c_num > 120 ||
                   input.egfr_num < 45 || input.egfr_num > 300 ||
                   input.alt_num < 0 || input.alt_num > 200)",
                   
                   # title of panel
                   plotOutput("effect_histogram", height = "400px", width = "800px")
  ),
  
  
)


#:--------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  require(tidyverse)
  require(rms)
  require(patchwork)
  require(ggtext)
  
  # join patient's data
  patient <- reactive({
    # join patient's data
    error = 0
    patient <- NULL
    # therapy
    patient$drugclass <- factor("SGLT2", levels = c("GLP1", "SGLT2"))
    # sex
    patient$sex <- factor(input$sex_select, levels = c("Female", "Male"))
    # agetx
    if (input$age_num < 18 || input$age_num > 99) {error= error + 1}
    patient$agetx <- as.numeric(input$age_num)
    # t2dmduration
    patient$t2dmduration <- as.numeric(input$t2dmduration_num)
    # prebmi
    patient$prebmi <- as.numeric(input$bmi_num)
    # prehba1cmmol
    patient$prehba1c <- as.numeric(input$hba1c_num)
    # egfr
    patient$preegfr <- as.numeric(input$egfr_num)
    # prealt
    patient$prealt <- as.numeric(input$alt_num)
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
    # hba1cmonth
    patient$hba1cmonth <- as.numeric(12)
    # drugline
    patient$drugline <- factor("2", levels = c("2", "3", "4", "5+"))
    # ncurrtx
    patient$ncurrtx <- factor("1", levels = c("1", "2", "3", "4", "5+"))
    
    # turn into data.frame
    patient <- as.data.frame(patient)
  })
  
  
  # Posteriors treatment effect
  posterior_effect <- reactive({
    
    # # Show a message to the user that their data is being processed (NOT NEEDED YET SINCE THE MODEL IS FAST)
    # showModal(modalDialog("Your data is being processed. Please wait...", 
    #                       footer = NULL,
    #                       tags$head(
    #                         tags$style(HTML("
    #                           .modal {
    #                             position: fixed;
    #                             top: 50%;
    #                             left: 50%;
    #                             transform: translate(-50%, -50%);
    #                           }
    #                         "))
    #                       )))
    
    # load patient
    patient <- patient()
    
    patient.matrix <- model.matrix(~ rcs(patient$agetx, c(45.06724, 58.54890, 71.66600)) + sex + ncurrtx + rms::rcs(prehba1c, c(60, 74, 99)) + rms::rcs(prebmi, c(26.6, 33.5, 43.4)) + rms::rcs(preegfr, c(72.43617, 97.22303, 111.99286)) + preheartfailure + preihd + preneuropathy + prepad + preretinopathy, data = patient) %>%
      as.matrix()
    
    # calculate effects
    effects <- patient.matrix %*% t(readRDS("best_linear_projection.rds"))
    
    effects
    
  })
  
  # # hide the message when the calculation is done
  # observeEvent(posterior_effect(), {
  #   removeModal()
  # })
  
  # Probability of benefit on therapy
  output$effect_probability <- renderPlot({

    posterior_effect <- posterior_effect()

    # SGLT2 probability of achieving target and surpassing it
    probability_SGLT2 <- length(which(posterior_effect < 0))/length(posterior_effect)

    # GLP1 probability of achieving target and surpassing it
    probability_GLP1 <- length(which(posterior_effect > 0))/length(posterior_effect)

    # reshape data
    df <- data.frame(
      SGLT2 = round(probability_SGLT2*100),
      GLP1 = round(probability_GLP1*100)
    ) %>%
      gather(key, val) %>%
      mutate(
        key = factor(key, levels = c("GLP1", "SGLT2"), labels = c("GLP-1RA", "SGLT2i")),
        Total = 100) %>%
      mutate(label = paste0(key, ": ", val,"%"))

    if (probability_SGLT2 > 0.999) {
      df$label[1] <- "SGLT2i: > 99.9%"
    } else if (probability_SGLT2 < 0.001) {
      df$label[1] <- "SGLT2i: <0.1%"
    }
    if (probability_GLP1 > 0.999) {
      df$label[2] <- "GLP-1RA: > 99.9%"
    } else if (probability_GLP1 < 0.001) {
      df$label[2] <- "GLP-1RA: <0.1%"
    }


    plot_effects <- df %>%
      filter(key == "SGLT2i") %>%
      ggplot(aes(key, val)) +
      geom_col(aes(y = Total), fill = "dodgerblue2", alpha = 0.6) +
      geom_col(fill = c("#f1a340"), alpha = 1) +
      coord_flip() +
      labs(
        title = "Predicted optimal therapy",
        subtitle = paste0("<span style='color:#f1a340;'>", df$label[1], "</span>  <span style='color:dodgerblue2;'>", df$label[2], "</span>")
      ) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        title = element_text(size = 20),
        plot.subtitle = element_markdown(lineheight = 1.1))



    plot_effects
    
  })


  # Probability of benefit on therapy
  output$effect_histogram <- renderPlot({

    posterior_effect <- posterior_effect()
    
    # SGLT2 probability of achieving target and surpassing it
    probability_SGLT2 <- length(which(posterior_effect < 0))/length(posterior_effect)
    
    # GLP1 probability of achieving target and surpassing it
    probability_GLP1 <- length(which(posterior_effect > 0))/length(posterior_effect)
    
    # reshape data
    df <- data.frame(
      SGLT2 = round(probability_SGLT2*100),
      GLP1 = round(probability_GLP1*100)
    ) %>%
      gather(key, val) %>%
      mutate(
        key = factor(key, levels = c("GLP1", "SGLT2"), labels = c("GLP-1RA", "SGLT2i")),
        Total = 100) %>%
      mutate(label = paste0(key, ": ", val,"%"))
    
    if (probability_SGLT2 > 0.999) {
      df$label[1] <- "SGLT2i: > 99.9%"
    } else if (probability_SGLT2 < 0.001) {
      df$label[1] <- "SGLT2i: <0.1%"
    }
    if (probability_GLP1 > 0.999) {
      df$label[2] <- "GLP-1RA: > 99.9%"
    } else if (probability_GLP1 < 0.001) {
      df$label[2] <- "GLP-1RA: <0.1%"
    }
    
    dat <- t(posterior_effect) %>%
      as.data.frame() %>%
      set_names("effects") %>%
      mutate(above= ifelse(effects < 0, "Favours SGLT2i", "Favours GLP1-RA")) %>%
      mutate(above = factor(above, levels = c("Favours SGLT2i", "Favours GLP1-RA")))
    
    # plot
    plot <- ggplot(data = dat, aes(x = effects, fill = above)) +
      geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(floor(min(posterior_effect)), ceiling(max(posterior_effect)), by = 0.5)) +
      geom_vline(aes(xintercept = 0), linetype = "dashed")+
      labs(
        title = "Predicted treatment benefit",
        subtitle = paste0("<span style='color:#f1a340;'>Negative value corresponds to a benefit on ", df$label[1],"</span><br><span style='color:dodgerblue2;'>Positive value corresponds to a benefit on ", df$label[2],"</span>"),
        x = "HbA1c difference (mmol/mol)"
      ) +
      scale_fill_manual(values = c("#f1a340", "dodgerblue2"))+
      theme_classic() +
      theme(legend.title = element_blank(),
            legend.direction = "horizontal",
            legend.position = "bottom",
            legend.box = "horizontal",
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            title = element_text(size = 20),
            plot.subtitle = element_markdown(lineheight = 1.1))

    plot


  })


  
}
  
  

# shinyApp(ui, server)
runGadget(ui, server, viewer = browserViewer(browser = getOption("browser")))











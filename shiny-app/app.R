################################################################################
################################################################################
## Shiny App for the Any-HPV Model ###
################################################################################
################################################################################

# Author: Babatunde Y. Alli
# Date: May 19, 2022


library(recipes)
library(rms)
library(shiny)

# Define UI for the application 
ui <- fluidPage(

        # Application title
        titlePanel("Oral HPV Trial Enrichment Model"),
    
        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(
            h3("Input Values"),
            sliderInput("age",
                        "Age",
                        min = 10,
                        max = 90,
                        value = 35),
            selectInput("gender",
                        "Sex",
                        choices = list("Male","Female")),
            selectInput("race",
                        "Race",
                        choices = list("Mexican American", "Non-Hispanic Black", 
                                       "Non-Hispanic White", "Other Hispanic", 
                                       "Other Race")),
            selectInput("marital_status",
                        "Marital status",
                        choices = list("Married/Living with partner", 
                                       "Never married", 
                                       "Widowed/Divorced/Separated")),
            sliderInput("smk_pk_yrs",
                        "Smoking pack-years",
                        min = 0,
                        max = 300,
                        value = 10),
            sliderInput("age_first_sex",
                        "Age at sexual debut",
                        min = 10,
                        max = 90,
                        value = 13),
            sliderInput("lifetime_num_oral_sex_partner",
                        "Lifetime number of people the individual 
                        has performed oral sex on",
                        min = 0,
                        max = 200,
                        value = 5),
            sliderInput("lifetime_num_sex_partner",
                        "Lifetime number of sexual partners",
                        min = 0,
                        max = 200,
                        value = 5),
            sliderInput("p_thresh",
                        "Probability threshold based on the decision curve",
                        min = 0,
                        max = 0.4,
                        value = 0.17)
          ),
          
    
          # Show predicted probability
          mainPanel(style = "flex-direction: column",
            div(style = "display: inline-block; width: 100%;",
                fluidRow(
                  column(12, imageOutput("dcurve", height = "50%"))
                )
            ),
            div(style = "display: inline-block; width: 100%;",
                fluidRow(
                  column(12, wellPanel(p(h3(textOutput("text")))))
                )
            )
          )
      )
)

# Define server logic
server <- function(input, output) {
            model <- readRDS("penalized_model.RDS")
            prepare <- readRDS("prepare.RDS")
      
            set.seed(12345)
            
            input_df <- reactive({
              df <- data.frame(age_first_sex = input$age_first_sex, 
                               lifetime_num_oral_sex_partner = input$lifetime_num_oral_sex_partner, 
                               lifetime_num_sex_partner = input$lifetime_num_sex_partner, 
                               gender = ifelse(input$gender == "Male", 1, 0), 
                               age = input$age, 
                               race = factor(c(input$race), levels = c("Mexican American", 
                                                                       "Non-Hispanic Black", 
                                                                       "Non-Hispanic White", 
                                                                       "Other Hispanic", "Other Race")
                               ), 
                               marital_status=factor(c(input$marital_status), 
                                                     levels = c("Married/Living with partner", 
                                                                "Never married", 
                                                                "Widowed/Divorced/Separated")
                               ), 
                               smk_pk_yrs=input$smk_pk_yrs)
              baked_df <- bake(prepare, new_data = df) #pre-process data
              as.data.frame(baked_df)
              
            })
      
            pred <- reactive({
              lp <- predict(model, input_df(), type="lp")
              plogis(lp) # from linear predictor to probabilities
            })
            
            output$text <- renderText({
              pred = pred()
              if (pred >= input$p_thresh) {
                decision_txt <- paste("This probability is greater or equal to the selected 
                                    probability threshold of:", 
                                    paste(round(input$p_thresh, 2), ".", sep = ""), 
                                    "\nTherefore, the recommended decision is to recruit.")
              } else {
                decision_txt <- paste("This probability is less than the selected 
                                    probability threshold of:", 
                                    paste(round(input$p_thresh, 2), ".", sep = ""), 
                                    "\nTherefore, the recommended decision is NOT to recruit.")
              }
              
              text =  paste("The predicted probability of oral HPV is: ", 
                                         paste(round(pred, 2), ".", sep = ""), 
                                         decision_txt, sep = "\n\n\n\n")
              text
              
            })
            output$dcurve <- renderImage({
              
              list(src = "dcurve.svg")
              
            }, deleteFile = F)
            
  
}

# Run the application 
shinyApp(ui = ui, server = server)

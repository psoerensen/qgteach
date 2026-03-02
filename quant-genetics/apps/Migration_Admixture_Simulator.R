library(shiny)
library(ggplot2)

ui <- fluidPage(
  
  titlePanel("Migration and Admixture Simulator"),
  withMathJax(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput("p1_0", "Initial allele frequency in Population 1 (p₁₀):",
                  min = 0.01, max = 0.99, value = 0.2, step = 0.01),
      
      sliderInput("p2_0", "Initial allele frequency in Population 2 (p₂₀):",
                  min = 0.01, max = 0.99, value = 0.8, step = 0.01),
      
      sliderInput("m", "Migration rate (m):",
                  min = 0, max = 0.5, value = 0.05, step = 0.01),
      
      sliderInput("generations", "Number of generations:",
                  min = 5, max = 200, value = 50),
      
      actionButton("run", "Run Simulation"),
      
      br(), hr(),
      actionButton("toggleDetails", "Show explanation")
    ),
    
    mainPanel(
      
      plotOutput("migrationPlot", height = "450px"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) Model Used"),
        
        helpText("Allele frequency update equations:"),
        helpText("$$p_1' = (1 - m)p_1 + m p_2$$"),
        helpText("$$p_2' = (1 - m)p_2 + m p_1$$"),
        
        p("A fraction m of individuals in each population are migrants from the other population."),
        
        hr(),
        
        h3("2) What to Explore"),
        
        p(strong("Small migration rate:")),
        p("Allele frequencies converge slowly."),
        
        p(strong("Large migration rate:")),
        p("Allele frequencies converge rapidly."),
        
        p(strong("Very different starting frequencies:")),
        p("Notice how gene flow reduces divergence."),
        
        hr(),
        
        h3("3) Key Principle"),
        
        p("Migration increases genetic similarity between populations."),
        p("It counteracts divergence caused by drift or selection."),
        p("Admixture is the genomic signature of gene flow.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    p1 <- numeric(input$generations + 1)
    p2 <- numeric(input$generations + 1)
    
    p1[1] <- input$p1_0
    p2[1] <- input$p2_0
    
    for (t in 1:input$generations) {
      p1[t + 1] <- (1 - input$m) * p1[t] + input$m * p2[t]
      p2[t + 1] <- (1 - input$m) * p2[t] + input$m * p1[t]
    }
    
    data.frame(
      Generation = 0:input$generations,
      Population1 = p1,
      Population2 = p2
    )
  })
  
  output$migrationPlot <- renderPlot({
    
    req(results())
    
    df_long <- reshape2::melt(results(), id.vars = "Generation")
    
    ggplot(df_long, aes(x = Generation, y = value, color = variable)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      ylim(0, 1) +
      labs(title = "Allele Frequency Change Under Migration",
           y = "Allele frequency (p)",
           color = "Population")
  })
  
  observeEvent(input$toggleDetails, {
    if (input$toggleDetails %% 2 == 1) {
      updateActionButton(session, "toggleDetails",
                         label = "Hide explanation")
    } else {
      updateActionButton(session, "toggleDetails",
                         label = "Show explanation")
    }
  })
}

shinyApp(ui, server)
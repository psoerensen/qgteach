library(shiny)
library(ggplot2)

ui <- fluidPage(
  
  titlePanel("Natural Selection Simulator (Viability Selection)"),
  withMathJax(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      sliderInput("p0", "Initial allele frequency (p₀):",
                  min = 0.01, max = 0.99, value = 0.2, step = 0.01),
      
      sliderInput("wAA", "Fitness w_AA:",
                  min = 0, max = 1.5, value = 1, step = 0.01),
      
      sliderInput("wAa", "Fitness w_Aa:",
                  min = 0, max = 1.5, value = 1, step = 0.01),
      
      sliderInput("waa", "Fitness w_aa:",
                  min = 0, max = 1.5, value = 0.5, step = 0.01),
      
      sliderInput("generations", "Number of generations:",
                  min = 10, max = 200, value = 50),
      
      actionButton("run", "Run Simulation"),
      
      br(), hr(),
      actionButton("toggleDetails", "Show explanation")
    ),
    
    mainPanel(
      
      plotOutput("allelePlot", height = "450px"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) Model Used"),
        
        helpText("Mean fitness:"),
        helpText("$$\\bar{w} = p^2 w_{AA} + 2pq w_{Aa} + q^2 w_{aa}$$"),
        
        helpText("Allele frequency update:"),
        helpText("$$p' = \\frac{p^2 w_{AA} + p q w_{Aa}}{\\bar{w}}$$"),
        
        hr(),
        
        h3("2) What to Explore"),
        
        p(strong("Directional selection:")),
        p("Make one genotype have highest fitness."),
        
        p(strong("Recessive beneficial allele:")),
        p("Set w_AA = w_Aa > w_aa."),
        
        p(strong("Overdominance (balancing selection):")),
        p("Set w_Aa highest."),
        
        p(strong("Underdominance:")),
        p("Set w_Aa lowest."),
        
        hr(),
        
        h3("3) Connection to Other Apps"),
        
        p("Unlike drift, this change is deterministic."),
        p("Selection changes allele frequencies through differential fitness."),
        p("Hardy–Weinberg proportions are reshaped each generation before selection acts.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    p <- numeric(input$generations + 1)
    wbar <- numeric(input$generations + 1)
    
    p[1] <- input$p0
    
    for (t in 1:input$generations) {
      
      q <- 1 - p[t]
      
      # Mean fitness
      wbar[t] <- p[t]^2 * input$wAA +
        2 * p[t] * q * input$wAa +
        q^2 * input$waa
      
      # Update allele frequency
      p[t + 1] <- (p[t]^2 * input$wAA +
                     p[t] * q * input$wAa) / wbar[t]
    }
    
    data.frame(
      Generation = 0:input$generations,
      p = p
    )
  })
  
  output$allelePlot <- renderPlot({
    
    req(results())
    
    ggplot(results(), aes(x = Generation, y = p)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      ylim(0, 1) +
      labs(title = "Allele Frequency Change Under Selection",
           y = "Allele frequency (p)",
           x = "Generation")
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
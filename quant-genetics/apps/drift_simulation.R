library(shiny)
library(ggplot2)

# Wright–Fisher simulation
simulate_drift <- function(N, p0, generations) {
  p <- numeric(generations + 1)
  p[1] <- p0
  
  for (t in 1:generations) {
    X <- rbinom(1, size = 2 * N, prob = p[t])
    p[t + 1] <- X / (2 * N)
  }
  p
}

ui <- fluidPage(
  titlePanel("Genetic Drift Simulator (Wright–Fisher Model)"),
  withMathJax(),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("N", "Population Size (N):",
                  min = 5, max = 1000, value = 50, step = 5),
      
      sliderInput("p0", "Initial Allele Frequency (p₀):",
                  min = 0.01, max = 0.99, value = 0.5, step = 0.01),
      
      sliderInput("generations", "Number of Generations:",
                  min = 10, max = 200, value = 50),
      
      sliderInput("replicates", "Number of Replicate Populations:",
                  min = 1, max = 50, value = 10),
      
      actionButton("run", "Run Simulation"),
      
      br(), hr(),
      
      actionButton("toggleDetails", "Show detailed explanation")
    ),
    
    mainPanel(
      plotOutput("driftPlot", height = "500px"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) How the Data Are Simulated"),
        p("The simulation follows the Wright–Fisher model of genetic drift."),
        helpText("$$X_{t+1} \\sim \\text{Binomial}(2N, p_t)$$"),
        helpText("$$p_{t+1} = \\frac{X_{t+1}}{2N}$$"),
        p("In each generation, 2N gametes are sampled at random from the previous generation."),
        p("The number of copies of the allele in the next generation follows a binomial distribution."),
        p("Each colored line in the plot represents one independently simulated population."),
        
        br(),
        
        h4("Expected Heterozygosity Under Drift"),
        helpText("$$H_t = H_0 \\left(1 - \\frac{1}{2N}\\right)^t$$"),
        p("Heterozygosity declines over time because random sampling gradually removes variation."),
        plotOutput("heteroPlot", height = "300px"),
        
        br(),
        verbatimTextOutput("summaryText"),
        
        hr(),
        
        h3("2) How to Use the App to Explore Drift"),
        p(strong("Change population size (N):")),
        p("Compare small and large populations to see how drift strength changes."),
        p(strong("Change the initial allele frequency (p₀):")),
        p("Observe how rare alleles are often lost quickly."),
        p(strong("Increase the number of replicates:")),
        p("Notice how each population evolves differently due to chance."),
        p("Use the simulation to connect visual patterns to the theoretical expectations shown in the slides."),
        
        hr(),
        
        h3("3) General Patterns to Expect"),
        p(strong("Drift is stronger in small populations.")),
        p("Allele frequencies fluctuate more and fixation happens faster when N is small."),
        p(strong("Rare alleles are likely to be lost.")),
        p("When p₀ is small, most populations lose the allele by chance."),
        p(strong("Heterozygosity declines over time.")),
        p("Genetic variation decreases because alleles are randomly fixed or lost."),
        p(strong("No selection is acting.")),
        p("All changes occur purely due to random sampling."),
        p("Hardy–Weinberg equilibrium assumes an infinitely large population. This simulation illustrates what happens when that assumption is violated.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Run simulations when requested
  results <- eventReactive(input$run, {
    replicate(input$replicates,
              simulate_drift(input$N, input$p0, input$generations))
  })
  
  output$driftPlot <- renderPlot({
    req(results())
    
    df <- data.frame(
      Generation = rep(0:input$generations, input$replicates),
      Frequency  = as.vector(results()),
      Replicate  = factor(rep(seq_len(input$replicates),
                              each = input$generations + 1))
    )
    
    ggplot(df, aes(x = Generation, y = Frequency, group = Replicate, color = Replicate)) +
      geom_line(alpha = 0.7) +
      theme_minimal() +
      ylim(0, 1) +
      labs(title = "Allele Frequency Trajectories",
           y = "Allele Frequency (p)",
           x = "Generation") +
      theme(legend.position = "none")
  })
  
  output$heteroPlot <- renderPlot({
    req(results())
    
    generations <- 0:input$generations
    
    p_matrix <- results()
    H_sim <- rowMeans(2 * p_matrix * (1 - p_matrix))
    
    H0 <- 2 * input$p0 * (1 - input$p0)
    H_theory <- H0 * (1 - 1/(2 * input$N))^generations
    
    df <- data.frame(
      Generation = generations,
      Simulated  = H_sim,
      Theory     = H_theory
    )
    
    ggplot(df, aes(x = Generation)) +
      geom_line(aes(y = Simulated), linewidth = 1) +
      geom_line(aes(y = Theory), linetype = "dashed", linewidth = 1) +
      theme_minimal() +
      labs(title = "Decay of Expected Heterozygosity",
           y = "Heterozygosity (H)",
           x = "Generation")
  })
  
  output$summaryText <- renderPrint({
    req(results())
    
    final_freq <- results()[input$generations + 1, ]
    fixation_prob <- mean(final_freq == 1)
    loss_prob     <- mean(final_freq == 0)
    
    cat("Fixation probability (p = 1):", round(fixation_prob, 3), "\n")
    cat("Loss probability (p = 0):", round(loss_prob, 3), "\n")
    cat("\nExpected fixation probability under neutrality:", round(input$p0, 3))
  })
  
  # Toggle button label (Show/Hide)
  observeEvent(input$toggleDetails, {
    if (input$toggleDetails %% 2 == 1) {
      updateActionButton(session, "toggleDetails", label = "Hide detailed explanation")
    } else {
      updateActionButton(session, "toggleDetails", label = "Show detailed explanation")
    }
  })
}

shinyApp(ui, server)

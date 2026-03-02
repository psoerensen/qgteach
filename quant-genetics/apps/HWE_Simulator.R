library(shiny)
library(ggplot2)

# Simulate one generation of random mating
simulate_random_mating <- function(N, p) {
  # Sample 2N gametes
  gametes <- rbinom(2 * N, 1, p)
  
  # Pair into individuals
  genotypes <- gametes[seq(1, 2*N, 2)] + gametes[seq(2, 2*N, 2)]
  
  AA <- sum(genotypes == 2)
  Aa <- sum(genotypes == 1)
  aa <- sum(genotypes == 0)
  
  c(AA = AA, Aa = Aa, aa = aa)
}

ui <- fluidPage(
  
  titlePanel("Hardy–Weinberg Simulation: Random Mating"),
  withMathJax(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      numericInput("AA0", "Initial AA:", 5, min = 0),
      numericInput("Aa0", "Initial Aa:", 8, min = 0),
      numericInput("aa0", "Initial aa:", 3, min = 0),
      
      numericInput("generations", "Number of generations:", 10, min = 1),
      
      br(),
      actionButton("run", "Run Simulation"),
      
      br(), hr(),
      actionButton("toggleDetails", "Show explanation")
    ),
    
    mainPanel(
      
      plotOutput("genoPlot", height = "450px"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) How the Simulation Works"),
        
        helpText("Allele frequency is estimated from the initial population:"),
        helpText("$$p = \\frac{2N_{AA} + N_{Aa}}{2N}$$"),
        
        helpText("Gametes are sampled randomly:"),
        helpText("$$X \\sim \\text{Binomial}(2N, p)$$"),
        
        helpText("Gametes are paired at random to form zygotes."),
        
        hr(),
        
        h3("2) What You Should Observe"),
        
        p(strong("After one generation of random mating:")),
        p("Genotype frequencies approach:"),
        helpText("$$p^2, \\quad 2pq, \\quad q^2$$"),
        
        p(strong("Subsequent generations:")),
        p("Genotype frequencies remain stable if no other forces act."),
        
        hr(),
        
        h3("3) Connection to Other Apps"),
        
        p("This illustrates the Hardy–Weinberg null model."),
        p("The HWE test compares observed data to these expectations."),
        p("The drift app shows what happens when the infinite population assumption is violated.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    N <- input$AA0 + input$Aa0 + input$aa0
    if (N == 0) return(NULL)
    
    # Initial allele frequency
    p <- (2 * input$AA0 + input$Aa0) / (2 * N)
    
    geno_matrix <- matrix(NA, nrow = input$generations + 1, ncol = 3)
    colnames(geno_matrix) <- c("AA", "Aa", "aa")
    
    # Generation 0
    geno_matrix[1, ] <- c(input$AA0, input$Aa0, input$aa0)
    
    for (g in 1:input$generations) {
      geno_matrix[g + 1, ] <- simulate_random_mating(N, p)
    }
    
    geno_matrix
  })
  
  output$genoPlot <- renderPlot({
    
    req(results())
    
    df <- as.data.frame(results())
    df$Generation <- 0:input$generations
    
    df_long <- reshape2::melt(df, id.vars = "Generation")
    
    ggplot(df_long, aes(x = Generation, y = value, color = variable)) +
      geom_line(linewidth = 1) +
      theme_minimal() +
      labs(title = "Genotype Counts Over Generations",
           y = "Count",
           color = "Genotype")
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


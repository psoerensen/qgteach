library(shiny)
library(ggplot2)

ui <- fluidPage(
  
  titlePanel("Hardy–Weinberg Equilibrium Test Explorer"),
  withMathJax(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      numericInput("AA", "Observed count (A/A):", value = 17, min = 0),
      numericInput("AG", "Observed count (A/G):", value = 55, min = 0),
      numericInput("GG", "Observed count (G/G):", value = 12, min = 0),
      
      br(),
      actionButton("run", "Run HWE Test"),
      
      br(), hr(),
      actionButton("toggleDetails", "Show detailed explanation")
    ),
    
    mainPanel(
      
      h3("Observed vs Expected Genotype Counts"),
      tableOutput("resultsTable"),
      
      br(),
      h3("Chi-Square Test Result"),
      verbatimTextOutput("chiResult"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) How the Test Is Computed"),
        
        helpText("Allele frequencies are estimated from the sample:"),
        helpText("$$p = \\frac{2N_{AA} + N_{AG}}{2N}$$"),
        helpText("$$q = 1 - p$$"),
        
        br(),
        
        helpText("Expected genotype frequencies under HWE:"),
        helpText("$$p^2, \\quad 2pq, \\quad q^2$$"),
        
        br(),
        
        helpText("Expected counts:"),
        helpText("$$E_{AA} = Np^2$$"),
        helpText("$$E_{AG} = N(2pq)$$"),
        helpText("$$E_{GG} = Nq^2$$"),
        
        br(),
        
        helpText("Chi-square statistic:"),
        helpText("$$\\chi^2 = \\sum \\frac{(O - E)^2}{E}$$"),
        
        hr(),
        
        h3("2) Interpretation"),
        
        p("The test has 1 degree of freedom."),
        p("If the p-value is small (e.g. < 0.05), we reject Hardy–Weinberg equilibrium."),
        
        p("Deviation from HWE may indicate:"),
        tags$ul(
          tags$li("Non-random mating"),
          tags$li("Selection"),
          tags$li("Population structure"),
          tags$li("Genetic drift"),
          tags$li("Genotyping error")
        ),
        
        hr(),
        
        h3("3) Key Concept"),
        
        p("The Hardy–Weinberg test compares observed genotype counts"),
        p("to expected counts under random mating."),
        p("It transforms a theoretical equilibrium into a formal statistical test.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  results <- eventReactive(input$run, {
    
    N <- input$AA + input$AG + input$GG
    
    if (N == 0) return(NULL)
    
    # Allele frequencies
    p <- (2 * input$AA + input$AG) / (2 * N)
    q <- 1 - p
    
    # Expected frequencies
    exp_AA <- p^2
    exp_AG <- 2 * p * q
    exp_GG <- q^2
    
    # Expected counts
    E_AA <- N * exp_AA
    E_AG <- N * exp_AG
    E_GG <- N * exp_GG
    
    # Chi-square contributions
    chi_AA <- (input$AA - E_AA)^2 / E_AA
    chi_AG <- (input$AG - E_AG)^2 / E_AG
    chi_GG <- (input$GG - E_GG)^2 / E_GG
    
    chi_total <- chi_AA + chi_AG + chi_GG
    
    p_value <- pchisq(chi_total, df = 1, lower.tail = FALSE)
    
    list(
      N = N,
      p = p,
      q = q,
      E_AA = E_AA,
      E_AG = E_AG,
      E_GG = E_GG,
      chi_AA = chi_AA,
      chi_AG = chi_AG,
      chi_GG = chi_GG,
      chi_total = chi_total,
      p_value = p_value
    )
  })
  
  output$resultsTable <- renderTable({
    
    req(results())
    
    data.frame(
      Genotype = c("A/A", "A/G", "G/G"),
      Observed = c(input$AA, input$AG, input$GG),
      Expected = round(c(results()$E_AA,
                         results()$E_AG,
                         results()$E_GG), 3),
      Chi_Square_Contribution = round(c(results()$chi_AA,
                                        results()$chi_AG,
                                        results()$chi_GG), 3)
    )
  })
  
  output$chiResult <- renderPrint({
    
    req(results())
    
    cat("Total sample size (N):", results()$N, "\n\n")
    
    cat("Estimated allele frequencies:\n")
    cat("p (A) =", round(results()$p, 4), "\n")
    cat("q (G) =", round(results()$q, 4), "\n\n")
    
    cat("Chi-square statistic:", round(results()$chi_total, 3), "\n")
    cat("Degrees of freedom: 1\n")
    cat("p-value:", round(results()$p_value, 5), "\n\n")
    
    if (results()$p_value < 0.05) {
      cat("Conclusion: Reject Hardy–Weinberg equilibrium (α = 0.05)")
    } else {
      cat("Conclusion: Do not reject Hardy–Weinberg equilibrium (α = 0.05)")
    }
  })
  
  observeEvent(input$toggleDetails, {
    if (input$toggleDetails %% 2 == 1) {
      updateActionButton(session, "toggleDetails",
                         label = "Hide detailed explanation")
    } else {
      updateActionButton(session, "toggleDetails",
                         label = "Show detailed explanation")
    }
  })
}

shinyApp(ui, server)

library(shiny)
library(ggplot2)

ui <- fluidPage(
  
  titlePanel("Gene Pool and Hardy–Weinberg Explorer"),
  withMathJax(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      numericInput("AA", "Number of AA individuals:", value = 5, min = 0),
      numericInput("Aa", "Number of Aa individuals:", value = 8, min = 0),
      numericInput("aa", "Number of aa individuals:", value = 3, min = 0),
      
      br(),
      actionButton("update", "Update Gene Pool"),
      
      br(), hr(),
      actionButton("toggleDetails", "Show detailed explanation")
    ),
    
    mainPanel(
      
      h3("Gene Pool Visualization"),
      plotOutput("genePoolPlot", height = "350px"),
      
      hr(),
      
      h3("Genotype Frequencies"),
      verbatimTextOutput("genoFreq"),
      
      h3("Allele Frequencies"),
      verbatimTextOutput("alleleFreq"),
      
      conditionalPanel(
        condition = "input.toggleDetails % 2 === 1",
        
        hr(),
        
        h3("1) How Frequencies Are Calculated"),
        
        helpText("Total population size:"),
        helpText("$$N = \\#AA + \\#Aa + \\#aa$$"),
        
        helpText("Total alleles:"),
        helpText("$$2N$$"),
        
        helpText("Allele frequency of A:"),
        helpText("$$p = \\frac{2\\#AA + \\#Aa}{2N}$$"),
        
        helpText("Allele frequency of a:"),
        helpText("$$q = \\frac{2\\#aa + \\#Aa}{2N}$$"),
        
        hr(),
        
        h3("2) How to Use This App"),
        
        p("Change genotype counts and observe how allele frequencies change."),
        p("Notice that different genotype combinations can produce the same allele frequencies."),
        p("This illustrates how the gene pool summarizes population composition."),
        
        hr(),
        
        h3("3) General Principles"),
        
        p(strong("Genotype frequencies sum to 1.")),
        p(strong("Allele frequencies sum to 1.")),
        p("The gene pool is described more compactly by allele frequencies."),
        p("These allele frequencies are the foundation of Hardy–Weinberg equilibrium.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  totals <- eventReactive(input$update, {
    
    N <- input$AA + input$Aa + input$aa
    total_alleles <- 2 * N
    
    A_count <- 2 * input$AA + input$Aa
    a_count <- 2 * input$aa + input$Aa
    
    p <- ifelse(total_alleles > 0, A_count / total_alleles, 0)
    q <- ifelse(total_alleles > 0, a_count / total_alleles, 0)
    
    list(
      N = N,
      total_alleles = total_alleles,
      A_count = A_count,
      a_count = a_count,
      p = p,
      q = q
    )
  })
  
  output$genePoolPlot <- renderPlot({
    
    req(totals())
    
    df <- data.frame(
      Genotype = c("AA", "Aa", "aa"),
      Count = c(input$AA, input$Aa, input$aa)
    )
    
    ggplot(df, aes(x = Genotype, y = Count)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Genotype Counts",
           y = "Number of Individuals",
           x = "")
  })
  
  output$genoFreq <- renderPrint({
    req(totals())
    
    if (totals()$N == 0) {
      cat("Population size is zero.")
      return()
    }
    
    f_AA <- input$AA / totals()$N
    f_Aa <- input$Aa / totals()$N
    f_aa <- input$aa / totals()$N
    
    cat("Population size (N):", totals()$N, "\n\n")
    cat("f(AA) =", round(f_AA, 3), "\n")
    cat("f(Aa) =", round(f_Aa, 3), "\n")
    cat("f(aa) =", round(f_aa, 3), "\n\n")
    cat("Check: sum =", round(f_AA + f_Aa + f_aa, 3))
  })
  
  output$alleleFreq <- renderPrint({
    req(totals())
    
    if (totals()$total_alleles == 0) {
      cat("No alleles present.")
      return()
    }
    
    cat("Total alleles (2N):", totals()$total_alleles, "\n\n")
    cat("Number of A alleles:", totals()$A_count, "\n")
    cat("Number of a alleles:", totals()$a_count, "\n\n")
    cat("p (A) =", round(totals()$p, 3), "\n")
    cat("q (a) =", round(totals()$q, 3), "\n\n")
    cat("Check: p + q =", round(totals()$p + totals()$q, 3))
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

library(shiny)
library(ggplot2)

# =========================================================
# TEXTBOOK HARD SWEEP (clean, dramatic)
# =========================================================
generate_textbook_sweep <- function(
    n_hap = 100,
    L = 201,
    sweep_freq = 1.0,
    core_width = 30,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(L %% 2 == 1)
  
  center <- (L + 1) / 2
  pos <- seq_len(L)
  dist <- abs(pos - center)
  
  H <- matrix(0, nrow = n_hap, ncol = L)
  
  core <- dist <= core_width
  is_sweep <- rbinom(n_hap, 1, sweep_freq) == 1
  
  for (i in seq_len(n_hap)) {
    H[i, ] <- rbinom(L, 1, 0.5)
    if (is_sweep[i]) {
      H[i, core] <- 1
      H[i, center] <- 1
    } else {
      H[i, center] <- 0
    }
  }
  
  list(H = H, center = center, pos = pos, dist = dist)
}

# =========================================================
# REALISTIC HARD SWEEP (noisy)
# =========================================================
generate_realistic_sweep <- function(
    n_hap = 100,
    L = 201,
    sweep_freq = 0.9,
    rec_scale = 30,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(L %% 2 == 1)
  
  center <- (L + 1) / 2
  pos <- seq_len(L)
  dist <- abs(pos - center)
  
  base_p <- runif(L, 0.2, 0.8)
  copy_prob <- exp(-dist / rec_scale)
  copy_prob[center] <- 1
  
  H <- matrix(0, nrow = n_hap, ncol = L)
  is_sweep <- rbinom(n_hap, 1, sweep_freq) == 1
  
  for (i in seq_len(n_hap)) {
    H[i, ] <- rbinom(L, 1, base_p)
    if (is_sweep[i]) {
      match <- rbinom(L, 1, copy_prob) == 1
      H[i, match] <- 1
      H[i, center] <- 1
    } else {
      H[i, center] <- 0
    }
  }
  
  list(H = H, center = center, pos = pos, dist = dist)
}

# =========================================================
# STATISTICS
# =========================================================
gene_diversity <- function(H) {
  p <- colMeans(H)
  2 * p * (1 - p)
}

ld_r2_with_center <- function(H, center) {
  x <- H[, center]
  px <- mean(x)
  vx <- px * (1 - px)
  r2 <- numeric(ncol(H))
  
  for (j in seq_len(ncol(H))) {
    y <- H[, j]
    py <- mean(y)
    vy <- py * (1 - py)
    if (vx == 0 || vy == 0) {
      r2[j] <- NA_real_
    } else {
      cov_xy <- mean((x - px) * (y - py))
      r <- cov_xy / sqrt(vx * vy)
      r2[j] <- r^2
    }
  }
  r2
}

# =========================================================
# UI
# =========================================================
ui <- fluidPage(
  titlePanel("Selective Sweep Visualizer"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode", "Simulation mode:",
                   choices = c("Textbook Hard Sweep", "Realistic Hard Sweep")),
      
      sliderInput("n_hap", "Number of haplotypes:", 40, 200, 100, 10),
      sliderInput("L", "Number of loci:", 101, 401, 201, 20),
      sliderInput("sweep_freq", "Beneficial allele frequency:", 0.6, 1, 1, 0.01),
      sliderInput("width", "Sweep width / recombination scale:", 5, 80, 30, 1),
      
      checkboxInput("sort_haps", "Cluster sweep haplotypes together", TRUE),
      
      actionButton("run", "Generate")
    ),
    
    mainPanel(
      h3("Haplotypes"),
      plotOutput("hapPlot", height = "400px"),
      hr(),
      h3("Gene Diversity"),
      plotOutput("gdPlot", height = "250px"),
      hr(),
      h3("LD with Selected Site (r²)"),
      plotOutput("ldPlot", height = "250px")
    )
  )
)

# =========================================================
# SERVER
# =========================================================
server <- function(input, output) {
  
  sweep_data <- eventReactive(input$run, {
    L <- as.integer(input$L)
    if (L %% 2 == 0) L <- L + 1
    
    if (input$mode == "Textbook Hard Sweep") {
      generate_textbook_sweep(
        n_hap = input$n_hap,
        L = L,
        sweep_freq = input$sweep_freq,
        core_width = input$width
      )
    } else {
      generate_realistic_sweep(
        n_hap = input$n_hap,
        L = L,
        sweep_freq = input$sweep_freq,
        rec_scale = input$width
      )
    }
  })
  
  output$hapPlot <- renderPlot({
    req(sweep_data())
    H <- sweep_data()$H
    center <- sweep_data()$center
    pos <- sweep_data()$pos
    
    # Optional sorting to make sweep block visually clean
    if (input$sort_haps) {
      H <- H[order(H[, center], decreasing = TRUE), ]
    }
    
    df <- data.frame(
      Hap = rep(seq_len(nrow(H)), each = ncol(H)),
      Pos = rep(seq_len(ncol(H)), times = nrow(H)),
      Allele = as.vector(H)
    )
    df$Allele <- as.vector(t(H))
    
    ggplot(df, aes(x = Pos, y = Hap, fill = factor(Allele))) +
      geom_tile() +
      geom_vline(xintercept = center, linetype = "dashed") +
      theme_minimal() +
      labs(fill = "Allele",
           title = "Haplotype structure around selected site")
  })
  
  output$gdPlot <- renderPlot({
    req(sweep_data())
    H <- sweep_data()$H
    center <- sweep_data()$center
    dist <- abs(seq_len(ncol(H)) - center)
    
    gd <- gene_diversity(H)
    df <- data.frame(Distance = dist, GD = gd)
    
    ggplot(df, aes(x = Distance, y = GD)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      ylim(0, 0.5) +
      labs(title = "Diversity trough around selected site",
           y = "2p(1-p)")
  })
  
  output$ldPlot <- renderPlot({
    req(sweep_data())
    H <- sweep_data()$H
    center <- sweep_data()$center
    dist <- abs(seq_len(ncol(H)) - center)
    
    r2 <- ld_r2_with_center(H, center)
    df <- data.frame(Distance = dist, r2 = r2)
    
    ggplot(df, aes(x = Distance, y = r2)) +
      geom_line(linewidth = 1.2) +
      theme_minimal() +
      ylim(0, 1) +
      labs(title = "LD peak near selected site",
           y = expression(r^2))
  })
}

shinyApp(ui, server)


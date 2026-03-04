plot_fitness_landscape <- function(
    wAA = 1,
    wAa = 0.9,
    waa = 0.8,
    base_size = 18
){
  
  p <- seq(0,1,length.out = 500)
  q <- 1 - p
  
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  
  df <- data.frame(
    p = p,
    wbar = wbar
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(p, wbar)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.6,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(min(wbar)*0.98, max(wbar)*1.02)
    ) +
    
    ggplot2::labs(
      title = "Fitness landscape",
      subtitle = paste0(
        "wAA = ", wAA,
        "  |  wAa = ", wAa,
        "  |  waa = ", waa
      ),
      x = "Allele frequency (p)",
      y = "Mean population fitness (w̄)"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_fitness_landscape(
  wAA = 1,
  wAa = 0.9,
  waa = 0.8
)

plot_fitness_landscape(
  wAA = 0.8,
  wAa = 1,
  waa = 0.8
)

plot_fitness_landscape(
  wAA = 1,
  wAa = 0.7,
  waa = 1
)


plot_adaptive_landscape <- function(
    peaks = list(
      c(0.25,0.7),
      c(0.7,0.9)
    ),
    widths = c(0.05,0.08),
    heights = c(1,0.8),
    base_size = 18
){
  
  x <- seq(0,1,length.out = 1000)
  
  fitness <- rep(0,length(x))
  
  for(i in seq_along(peaks)){
    
    mu <- peaks[[i]][1]
    h  <- peaks[[i]][2]
    w  <- widths[i]
    
    fitness <- fitness +
      h*exp(-(x-mu)^2/(2*w^2))
  }
  
  df <- data.frame(
    genotype_space = x,
    fitness = fitness
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(genotype_space,fitness)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.6,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::labs(
      title = "Adaptive landscape",
      subtitle = "Fitness peaks represent adaptive optima",
      x = "Genotype / phenotype space",
      y = "Fitness"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}


plot_adaptive_landscape()

plot_adaptive_landscape(
  peaks = list(
    c(0.2,0.8),
    c(0.5,0.6),
    c(0.8,0.9)
  ),
  widths = c(0.06,0.05,0.04),
  heights = c(0.9,0.7,1)
)

plot_wright_fisher_sampling <- function(
    N = 10,
    p = 0.6,
    seed = 1,
    base_size = 18
){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Parent generation
  # ------------------------------------------------------------
  
  alleles_parent <- c(
    rep("A", round(2*N*p)),
    rep("a", 2*N - round(2*N*p))
  )
  
  alleles_parent <- sample(alleles_parent)
  
  # ------------------------------------------------------------
  # Sample offspring generation
  # ------------------------------------------------------------
  
  alleles_offspring <- sample(
    alleles_parent,
    size = 2*N,
    replace = TRUE
  )
  
  # ------------------------------------------------------------
  # Helper to make grid positions
  # ------------------------------------------------------------
  
  make_grid <- function(alleles, y_offset){
    
    n <- length(alleles)
    
    data.frame(
      x = rep(seq_len(N), each = 2),
      y = y_offset,
      allele = alleles
    )
  }
  
  parent_df <- make_grid(alleles_parent, 1)
  off_df    <- make_grid(alleles_offspring, 0)
  
  df <- rbind(parent_df, off_df)
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x, y, fill = allele)
  ) +
    
    ggplot2::geom_point(
      shape = 21,
      size = 6,
      color = "black",
      stroke = 0.7
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$allele
    ) +
    
    ggplot2::annotate(
      "text",
      x = N/2,
      y = 1.35,
      label = "Parent generation",
      size = base_size*0.35
    ) +
    
    ggplot2::annotate(
      "text",
      x = N/2,
      y = -0.35,
      label = "Offspring generation (random sampling)",
      size = base_size*0.35
    ) +
    
    ggplot2::coord_cartesian(
      xlim = c(0.5, N+0.5),
      ylim = c(-0.5,1.5)
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

plot_wright_fisher_sampling(
  N = 12,
  p = 0.7
)



plot_fixation_probability <- function(
    N = 100,
    s_range = seq(-0.05, 0.1, length.out = 200),
    base_size = 18
){
  
  # Avoid numerical issues near s = 0
  pfix <- ifelse(
    abs(s_range) < 1e-8,
    1/(2*N),
    (1 - exp(-2*s_range)) / (1 - exp(-4*N*s_range))
  )
  
  df <- tibble::tibble(
    s = s_range,
    pfix = pfix
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(s, pfix)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.6,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = qg_cols$grey["dark"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = "Fixation probability of a new allele",
      subtitle = paste0("Population size N = ", N),
      x = "Selection coefficient (s)",
      y = "Fixation probability"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_fixation_probability(
  N = 100
)


plot_site_frequency_spectrum <- function(
    n = 20,                # sample size (chromosomes)
    n_sites = 2000,        # number of segregating sites
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Neutral expectation ~ 1/k
  # ------------------------------------------------------------
  
  k <- 1:(n-1)
  probs <- 1/k
  probs <- probs/sum(probs)
  
  counts <- rmultinom(1, n_sites, probs)
  
  df <- tibble::tibble(
    freq = k,
    sites = as.numeric(counts)
  )
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(freq, sites)
  ) +
    
    ggplot2::geom_col(
      fill = qg_cols$pop["A"],
      color = "black",
      width = 0.8
    ) +
    
    ggplot2::labs(
      title = "Site Frequency Spectrum",
      subtitle = paste0("Sample size n = ", n),
      x = "Derived allele count",
      y = "Number of sites"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_site_frequency_spectrum()

plot_site_frequency_spectrum(
  n = 50,
  n_sites = 5000
)

plot_coalescent_tree <- function(
    n = 6,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Simulate coalescent waiting times
  # ------------------------------------------------------------
  
  k <- n:2
  times <- stats::rexp(length(k), rate = choose(k,2))
  
  t_coal <- cumsum(times)
  
  # ------------------------------------------------------------
  # Build simple tree structure
  # ------------------------------------------------------------
  
  lineages <- list()
  
  for(i in seq_len(n)){
    lineages[[i]] <- list(x = i, y = 0)
  }
  
  edges <- data.frame()
  
  current_lineages <- 1:n
  heights <- rep(0,n)
  
  next_node <- n + 1
  
  for(i in seq_len(n-1)){
    
    pair <- sample(current_lineages,2)
    
    h <- t_coal[i]
    
    for(p in pair){
      
      edges <- rbind(edges,
                     data.frame(
                       x = p,
                       xend = p,
                       y = heights[p],
                       yend = h
                     )
      )
    }
    
    x_new <- mean(pair)
    
    edges <- rbind(edges,
                   data.frame(
                     x = pair[1],
                     xend = pair[2],
                     y = h,
                     yend = h
                   )
    )
    
    heights[pair] <- h
    heights[next_node] <- h
    
    current_lineages <- setdiff(current_lineages, pair)
    current_lineages <- c(current_lineages, next_node)
    
    next_node <- next_node + 1
  }
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  ggplot2::ggplot(edges) +
    
    ggplot2::geom_segment(
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 1.3,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::labs(
      title = "Coalescent genealogy",
      subtitle = paste0("Sample size n = ", n),
      x = "Sampled chromosomes",
      y = "Time (coalescent units)"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_coalescent_tree(n = 6)
plot_coalescent_tree(n = 12)


plot_coalescent_tree <- function(
    n = 6,
    mu = 0.5,                 # mutation rate per branch unit
    show_mutations = FALSE,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Simulate coalescent waiting times
  # ------------------------------------------------------------
  
  k <- n:2
  wait_times <- stats::rexp(length(k), rate = choose(k,2))
  t_coal <- cumsum(wait_times)
  
  # ------------------------------------------------------------
  # Build edges
  # ------------------------------------------------------------
  
  edges <- data.frame()
  
  lineages <- data.frame(
    id = 1:n,
    x = 1:n,
    y = 0
  )
  
  next_node <- n + 1
  
  mutation_df <- data.frame()
  
  for(i in seq_len(n-1)){
    
    pair <- sample(lineages$id,2)
    
    node_height <- t_coal[i]
    
    parents <- lineages[lineages$id %in% pair,]
    
    # vertical branches
    for(j in 1:2){
      
      branch_length <- node_height - parents$y[j]
      
      edges <- rbind(edges,
                     data.frame(
                       x = parents$x[j],
                       xend = parents$x[j],
                       y = parents$y[j],
                       yend = node_height
                     )
      )
      
      # --------------------------------------------------------
      # Mutations on branch
      # --------------------------------------------------------
      
      if(show_mutations){
        
        n_mut <- stats::rpois(1, mu * branch_length)
        
        if(n_mut > 0){
          
          mut_pos <- runif(n_mut,
                           parents$y[j],
                           node_height)
          
          mutation_df <- rbind(
            mutation_df,
            data.frame(
              x = parents$x[j],
              y = mut_pos
            )
          )
        }
      }
    }
    
    # horizontal join
    edges <- rbind(edges,
                   data.frame(
                     x = parents$x[1],
                     xend = parents$x[2],
                     y = node_height,
                     yend = node_height
                   )
    )
    
    # update lineage list
    new_x <- mean(parents$x)
    
    lineages <- lineages[!lineages$id %in% pair,]
    
    lineages <- rbind(
      lineages,
      data.frame(
        id = next_node,
        x = new_x,
        y = node_height
      )
    )
    
    next_node <- next_node + 1
  }
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p <- ggplot2::ggplot(edges) +
    
    ggplot2::geom_segment(
      ggplot2::aes(x = x, xend = xend,
                   y = y, yend = yend),
      linewidth = 1.4,
      color = qg_cols$pop["A"]
    )
  
  if(show_mutations && nrow(mutation_df) > 0){
    
    p <- p +
      ggplot2::geom_point(
        data = mutation_df,
        ggplot2::aes(x, y),
        size = 2.5,
        color = qg_cols$pop["B"]
      )
  }
  
  p +
    ggplot2::labs(
      title = if(show_mutations)
        "Coalescent genealogy with mutations"
      else
        "Coalescent genealogy",
      subtitle = paste0("Sample size n = ", n),
      x = "Sampled chromosomes",
      y = "Time (coalescent units)"
    ) +
    theme_qg(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_coalescent_tree(
  n = 12,
  show_mutations = TRUE
)


plot_coalescent_recombination <- function(
    n = 6,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  simulate_tree <- function(n){
    
    k <- n:2
    wait <- stats::rexp(length(k), rate = choose(k,2))
    t_coal <- cumsum(wait)
    
    edges <- data.frame()
    
    lineages <- data.frame(
      id = 1:n,
      x = 1:n,
      y = 0
    )
    
    next_node <- n + 1
    
    for(i in seq_len(n-1)){
      
      pair <- sample(lineages$id,2)
      
      parents <- lineages[lineages$id %in% pair,]
      h <- t_coal[i]
      
      for(j in 1:2){
        
        edges <- rbind(edges,
                       data.frame(
                         x = parents$x[j],
                         xend = parents$x[j],
                         y = parents$y[j],
                         yend = h
                       )
        )
      }
      
      edges <- rbind(edges,
                     data.frame(
                       x = parents$x[1],
                       xend = parents$x[2],
                       y = h,
                       yend = h
                     )
      )
      
      new_x <- mean(parents$x)
      
      lineages <- lineages[!lineages$id %in% pair,]
      
      lineages <- rbind(
        lineages,
        data.frame(
          id = next_node,
          x = new_x,
          y = h
        )
      )
      
      next_node <- next_node + 1
    }
    
    edges
  }
  
  tree1 <- simulate_tree(n)
  tree1$region <- "Region 1"
  
  tree2 <- simulate_tree(n)
  tree2$region <- "Region 2"
  
  edges <- rbind(tree1, tree2)
  
  edges$x <- ifelse(
    edges$region == "Region 1",
    edges$x,
    edges$x + n + 3
  )
  
  edges$xend <- ifelse(
    edges$region == "Region 1",
    edges$xend,
    edges$xend + n + 3
  )
  
  ggplot2::ggplot(edges) +
    
    ggplot2::geom_segment(
      ggplot2::aes(x = x, xend = xend,
                   y = y, yend = yend),
      linewidth = 1.3,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::annotate(
      "text",
      x = mean(1:n),
      y = max(edges$y) * 1.05,
      label = "Region 1"
    ) +
    
    ggplot2::annotate(
      "text",
      x = mean(1:n) + n + 3,
      y = max(edges$y) * 1.05,
      label = "Region 2"
    ) +
    
    ggplot2::annotate(
      "segment",
      x = n + 1,
      xend = n + 2,
      y = max(edges$y) * 0.8,
      yend = max(edges$y) * 0.8,
      linetype = "dashed"
    ) +
    
    ggplot2::annotate(
      "text",
      x = n + 1.5,
      y = max(edges$y) * 0.9,
      label = "Recombination breakpoint"
    ) +
    
    ggplot2::labs(
      title = "Recombination creates different genealogies along the genome",
      x = "Sampled chromosomes",
      y = "Coalescent time"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

plot_coalescent_recombination(n = 6)


plot_genealogy_vs_ld <- function(
    L = 50,
    r = 0.05,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Simulate LD decay along genome
  # ------------------------------------------------------------
  
  pos <- 1:L
  
  D <- exp(-r * pos)
  r2 <- D^2
  
  ld_df <- data.frame(
    position = pos,
    r2 = r2
  )
  
  # ------------------------------------------------------------
  # Simulate recombination breakpoints
  # ------------------------------------------------------------
  
  n_breaks <- sample(2:4,1)
  
  breaks <- sort(sample(5:(L-5), n_breaks))
  
  genealogy_df <- data.frame(
    start = c(1, breaks + 1),
    end = c(breaks, L),
    tree = paste0("Tree ", seq_len(n_breaks + 1))
  )
  
  genealogy_df$y <- 1
  
  # ------------------------------------------------------------
  # LD plot
  # ------------------------------------------------------------
  
  p_ld <- ggplot2::ggplot(ld_df,
                          ggplot2::aes(position, r2)) +
    
    ggplot2::geom_line(
      linewidth = 1.5,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::labs(
      title = "Linkage disequilibrium along the genome",
      x = "Genomic position",
      y = expression(r^2)
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------
  # Genealogy blocks
  # ------------------------------------------------------------
  
  p_gene <- ggplot2::ggplot(genealogy_df) +
    
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = 0,
        ymax = 1,
        fill = tree
      ),
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$pop[
        seq_len(nrow(genealogy_df))
      ]
    ) +
    
    ggplot2::labs(
      title = "Local genealogies along the genome",
      x = "Genomic position",
      y = NULL
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  p_ld / p_gene
}

plot_genealogy_vs_ld()


plot_selection_vs_drift()

plot_mutation_selection_balance()

plot_island_model()


plot_mutation_dynamics(
  p0 = 0.9,
  mu = 0.01,     # A -> a
  nu = 0.001,    # a -> A
  generations = 150
)

plot_mutation_dynamics(
  p0 = 1,
  mu = 0.01,
  nu = 0.01,
  generations = 100
)

plot_migration_dynamics(
  p1 = 0.1,
  p2 = 0.9,
  m = 0.05,
  generations = 100
)

plot_migration_dynamics(
  p1 = 0.1,
  p2 = 0.9,
  m = 0.2,
  generations = 40
)

plot_effective_population_size(
  Nm = 50,
  Nf = 50
)

plot_effective_population_size(
  Nm = 5,
  Nf = 95
)



plot_molecular_clock(
  mu = 1e-8,
  L = 1e7,
  generations = 300
)

plot_molecular_clock(
  mu = 1e-8,
  L = 1e7,
  generations = 300,
  stochastic = TRUE
)

plot_domestication_bottleneck(
  N_ancestral = 1000,
  N_bottleneck = 10,
  N_post = 200,
  bottleneck_start = 30,
  bottleneck_duration = 20,
  generations = 120
)



plot_mutation_dynamics <- function(
    p0 = 0.9,
    mu = 0.001,
    nu = 0.0001,
    generations = 200,
    base_size = 18
){
  
  p <- numeric(generations + 1)
  p[1] <- p0
  
  for(t in seq_len(generations)){
    p[t+1] <- p[t]*(1-mu) + (1-p[t])*nu
  }
  
  df <- tibble::tibble(
    generation = 0:generations,
    p = p
  )
  
  p_eq <- nu/(mu+nu)
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation, p)) +
    
    ggplot2::geom_line(
      linewidth = 1.5,
      color = qg_cols$pop["B"]
    ) +
    
    ggplot2::geom_hline(
      yintercept = p_eq,
      linetype = "dashed",
      color = qg_cols$grey["dark"]
    ) +
    
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    
    ggplot2::labs(
      title = "Mutation drives allele frequencies toward equilibrium",
      subtitle = paste0("Equilibrium p* = ", round(p_eq,4)),
      x = "Generation",
      y = "Allele frequency (p)"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_migration_dynamics <- function(
    p1 = 0.2,
    p2 = 0.8,
    m = 0.05,
    generations = 100,
    base_size = 18
){
  
  p1_traj <- numeric(generations + 1)
  p2_traj <- numeric(generations + 1)
  
  p1_traj[1] <- p1
  p2_traj[1] <- p2
  
  for(t in seq_len(generations)){
    p1_traj[t+1] <- (1-m)*p1_traj[t] + m*p2_traj[t]
    p2_traj[t+1] <- (1-m)*p2_traj[t] + m*p1_traj[t]
  }
  
  df <- tibble::tibble(
    generation = rep(0:generations,2),
    p = c(p1_traj,p2_traj),
    pop = factor(
      rep(c("Population 1","Population 2"),
          each = generations+1),
      levels = c("Population 1","Population 2")
    )
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation,p,color=pop)) +
    
    ggplot2::geom_line(linewidth=1.5) +
    
    ggplot2::scale_color_manual(
      values=c(
        "Population 1" = qg_cols$pop["A"],
        "Population 2" = qg_cols$pop["B"]
      )
    ) +
    
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    
    ggplot2::labs(
      title="Migration causes allele frequencies to converge",
      subtitle=paste0("Migration rate m = ",m),
      x="Generation",
      y="Allele frequency (p)",
      color=NULL
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_effective_population_size <- function(
    Nm = 10,
    Nf = 90,
    base_size = 18
){
  
  N <- Nm + Nf
  Ne <- (4*Nm*Nf)/(Nm+Nf)
  
  df <- tibble::tibble(
    type = factor(
      c("Census size (N)","Effective size (Ne)"),
      levels = c("Census size (N)","Effective size (Ne)")
    ),
    size = c(N,Ne)
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(type,size,fill=type)) +
    
    ggplot2::geom_col(width=0.6,color="black") +
    
    ggplot2::scale_fill_manual(
      values=c(
        "Census size (N)" = qg_cols$grey["light"],
        "Effective size (Ne)" = qg_cols$pop["A"]
      )
    ) +
    
    ggplot2::labs(
      title="Effective population size can be smaller than N",
      subtitle=paste0("Nm = ",Nm,", Nf = ",Nf),
      x=NULL,
      y="Population size"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position="none"
    )
}

plot_molecular_clock <- function(
    mu = 1e-8,
    L = 1e6,
    generations = 300,
    stochastic = FALSE,
    base_size = 18,
    seed = 1
){
  
  if(stochastic) set.seed(seed)
  
  t <- 0:generations
  
  if(!stochastic){
    
    divergence <- 2*mu*L*t
    
  } else {
    
    new_mut <- stats::rpois(generations,lambda=2*mu*L)
    divergence <- c(0,cumsum(new_mut))
    
  }
  
  df <- tibble::tibble(
    generation=t,
    divergence=divergence
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation,divergence)) +
    
    ggplot2::geom_line(
      linewidth=1.5,
      color=qg_cols$pop["B"]
    ) +
    
    ggplot2::labs(
      title="Molecular clock: divergence increases over time",
      subtitle=paste0("μ = ",mu,"  |  L = ",format(L,scientific=TRUE)),
      x="Time (generations)",
      y="Number of substitutions"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_domestication_bottleneck <- function(
    N_ancestral = 1000,
    N_bottleneck = 20,
    N_post = 200,
    bottleneck_start = 40,
    bottleneck_duration = 20,
    generations = 120,
    H0 = 0.5,
    base_size = 18
){
  
  H <- numeric(generations+1)
  H[1] <- H0
  
  for(t in seq_len(generations)){
    
    if(t < bottleneck_start){
      N_current <- N_ancestral
    } else if(t < bottleneck_start + bottleneck_duration){
      N_current <- N_bottleneck
    } else {
      N_current <- N_post
    }
    
    H[t+1] <- H[t]*(1-1/(2*N_current))
  }
  
  df <- tibble::tibble(
    generation=0:generations,
    H=H
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation,H)) +
    
    ggplot2::geom_line(
      linewidth=1.5,
      color=qg_cols$pop["A"]
    ) +
    
    ggplot2::geom_vline(
      xintercept=c(
        bottleneck_start,
        bottleneck_start+bottleneck_duration
      ),
      linetype="dashed",
      color=qg_cols$grey["dark"]
    ) +
    
    ggplot2::labs(
      title="Domestication bottleneck reduces genetic diversity",
      subtitle=paste0(
        "Ancestral N=",N_ancestral,
        " → Bottleneck N=",N_bottleneck
      ),
      x="Generation",
      y="Expected heterozygosity"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}


plot_selection_vs_drift <- function(
    p0 = 0.5,
    s = 0.1,
    N = 50,
    generations = 100,
    replicates = 20,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  # deterministic selection
  p_sel <- numeric(generations + 1)
  p_sel[1] <- p0
  
  for(t in seq_len(generations)){
    p <- p_sel[t]
    wbar <- p*(1+s) + (1-p)
    p_sel[t+1] <- (p*(1+s))/wbar
  }
  
  df_sel <- tibble::tibble(
    generation = 0:generations,
    p = p_sel
  )
  
  # drift simulations
  sim <- vector("list", replicates)
  
  for(r in seq_len(replicates)){
    
    p <- numeric(generations + 1)
    p[1] <- p0
    
    for(t in seq_len(generations)){
      p[t+1] <- stats::rbinom(1, 2*N, p[t])/(2*N)
    }
    
    sim[[r]] <- tibble::tibble(
      generation = 0:generations,
      p = p,
      replicate = r
    )
  }
  
  df_drift <- dplyr::bind_rows(sim)
  
  ggplot2::ggplot() +
    
    ggplot2::geom_line(
      data = df_drift,
      ggplot2::aes(generation, p, group = replicate),
      color = qg_cols$grey["light"],
      alpha = 0.5
    ) +
    
    ggplot2::geom_line(
      data = df_sel,
      ggplot2::aes(generation, p),
      linewidth = 1.6,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = "Selection vs genetic drift",
      subtitle = paste0("Selection coefficient s = ", s, "  |  Population size N = ", N),
      x = "Generation",
      y = "Allele frequency (p)"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_mutation_selection_balance <- function(
    q0 = 0.1,
    mu = 1e-5,
    s = 0.1,
    generations = 200,
    base_size = 18
){
  
  q <- numeric(generations + 1)
  q[1] <- q0
  
  for(t in seq_len(generations)){
    q[t+1] <- (q[t] - s*q[t]^2*(1-q[t])) + mu*(1-q[t])
  }
  
  df <- tibble::tibble(
    generation = 0:generations,
    q = q
  )
  
  q_eq <- sqrt(mu/s)
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation,q)) +
    
    ggplot2::geom_line(
      linewidth = 1.5,
      color = qg_cols$pop["B"]
    ) +
    
    ggplot2::geom_hline(
      yintercept = q_eq,
      linetype = "dashed",
      color = qg_cols$grey["dark"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0, max(df$q)*1.1)) +
    
    ggplot2::labs(
      title = "Mutation–selection balance",
      subtitle = paste0("Equilibrium q* ≈ √(μ/s) = ", signif(q_eq,3)),
      x = "Generation",
      y = "Frequency of deleterious allele (q)"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_island_model <- function(
    n_pop = 5,
    m = 0.05,
    generations = 100,
    base_size = 18,
    seed = 1
){
  
  set.seed(seed)
  
  p <- runif(n_pop)
  
  traj <- matrix(NA, nrow = generations + 1, ncol = n_pop)
  traj[1,] <- p
  
  for(t in seq_len(generations)){
    
    p_bar <- mean(traj[t,])
    
    traj[t+1,] <- (1-m)*traj[t,] + m*p_bar
  }
  
  df <- tibble::tibble(
    generation = rep(0:generations, n_pop),
    p = as.vector(traj),
    pop = factor(rep(1:n_pop, each = generations + 1))
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(generation,p,color=pop)) +
    
    ggplot2::geom_line(linewidth = 1.3) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = "Island model migration",
      subtitle = paste0("Migration rate m = ", m),
      x = "Generation",
      y = "Allele frequency (p)",
      color = "Population"
    ) +
    
    theme_qg(base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}


# ============================================================
# Plot Heterozygosity Decay
# ============================================================


plot_heterozygosity_decay(
  N = 20,
  generations = 100,
  H0 = 0.5
)

p_small <- plot_heterozygosity_decay(N = 20, generations = 100)
p_large <- plot_heterozygosity_decay(N = 200, generations = 100)

p_small | p_large


plot_heterozygosity_decay <- function(
    N = 50,
    generations = 100,
    H0 = 0.5,
    base_size = 18
) {
  
  H <- numeric(generations + 1)
  H[1] <- H0
  
  for (t in seq_len(generations)) {
    H[t + 1] <- H[t] * (1 - 1/(2 * N))
  }
  
  df <- tibble::tibble(
    generation = 0:generations,
    H = H
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(generation, H)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.5,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    
    ggplot2::labs(
      title = "Decay of heterozygosity under genetic drift",
      subtitle = paste0("Population size N = ", N),
      x = "Generation",
      y = "Expected heterozygosity (H)"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

# ============================================================
# Plot Selective Sweep
# ============================================================

psweep <- make_sweep_plots(
  mode = "Realistic Hard Sweep",
  n_hap = 100,
  L = 201,
  sweep_freq = 0.9,
  width = 30,
  seed = 1
)

psweep$hapPlot
psweep$gdPlot
psweep$ldPlot

make_sweep_plots <- function(
    mode = c("Textbook Hard Sweep", "Realistic Hard Sweep"),
    n_hap = 100,
    L = 201,
    sweep_freq = 1.0,
    width = 30,
    sort_haps = TRUE,
    seed = NULL,
    base_size = 18
) {
  
  mode <- match.arg(mode)
  
  L <- as.integer(L)
  if (L %% 2 == 0) L <- L + 1
  
  sweep_data <- if (mode == "Textbook Hard Sweep") {
    generate_textbook_sweep(
      n_hap = n_hap,
      L = L,
      sweep_freq = sweep_freq,
      core_width = width,
      seed = seed
    )
  } else {
    generate_realistic_sweep(
      n_hap = n_hap,
      L = L,
      sweep_freq = sweep_freq,
      rec_scale = width,
      seed = seed
    )
  }
  
  H <- sweep_data$H
  center <- sweep_data$center
  
  if (sort_haps) {
    H <- H[order(H[, center], decreasing = TRUE), , drop = FALSE]
  }
  
  # ------------------------------------------------------------
  # Haplotype matrix plot
  # ------------------------------------------------------------
  
  df_hap <- data.frame(
    Hap = rep(seq_len(nrow(H)), each = ncol(H)),
    Pos = rep(seq_len(ncol(H)), times = nrow(H)),
    Allele = as.vector(t(H))
  )
  
  # Convert numeric allele coding to A/a
  df_hap$Allele <- ifelse(df_hap$Allele == 1, "A", "a")
  
  df_hap$Allele <- factor(df_hap$Allele, levels = c("A","a"))
  
  hapPlot <- ggplot2::ggplot(
    df_hap,
    ggplot2::aes(Pos, Hap, fill = Allele)
  ) +
    
    ggplot2::geom_tile() +
    
    ggplot2::geom_vline(
      xintercept = center,
      linetype = "dashed",
      color = qg_cols$grey["dark"]
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$allele
    ) +
    
    ggplot2::labs(
      title = "Haplotype structure around selected site",
      x = "Genomic position",
      y = "Haplotype",
      fill = "Allele"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------
  # Gene diversity plot
  # ------------------------------------------------------------
  
  dist <- abs(seq_len(ncol(H)) - center)
  gd <- gene_diversity(H)
  
  df_gd <- data.frame(
    Distance = dist,
    GD = gd
  )
  
  gdPlot <- ggplot2::ggplot(
    df_gd,
    ggplot2::aes(Distance, GD)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.4,
      color = qg_cols$pop["A"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0, 0.5)) +
    
    ggplot2::labs(
      title = "Diversity trough around selected site",
      x = "Distance from selected site",
      y = "Gene diversity (2p(1-p))"
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------
  # LD decay around sweep
  # ------------------------------------------------------------
  
  r2 <- ld_r2_with_center(H, center)
  
  df_ld <- data.frame(
    Distance = dist,
    r2 = r2
  )
  
  ldPlot <- ggplot2::ggplot(
    df_ld,
    ggplot2::aes(Distance, r2)
  ) +
    
    ggplot2::geom_line(
      linewidth = 1.4,
      color = qg_cols$pop["B"]
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    
    ggplot2::labs(
      title = "Linkage disequilibrium near selected site",
      x = "Distance from selected site",
      y = expression(r^2)
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  list(
    hapPlot = hapPlot,
    gdPlot = gdPlot,
    ldPlot = ldPlot
  )
}


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


# ============================================================
# Plot HWE Test
# ============================================================

plot_hwe_comparison(c(AA = 10, Aa = 7, aa = 50))


plot_hwe_comparison <- function(
    geno_counts,
    base_size = 18
) {
  
  # ------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------
  
  if (is.null(names(geno_counts))) {
    stop("geno_counts must be named, e.g. c(AA=10, Aa=7, aa=10)")
  }
  
  required <- c("AA","Aa","aa")
  
  geno_counts[setdiff(required, names(geno_counts))] <- 0
  geno_counts <- geno_counts[required]
  
  N <- sum(geno_counts)
  
  # ------------------------------------------------------------
  # Allele frequency
  # ------------------------------------------------------------
  
  nA <- 2 * geno_counts["AA"] + geno_counts["Aa"]
  
  p_hat <- as.numeric(nA) / (2 * N)
  q_hat <- 1 - p_hat
  
  # ------------------------------------------------------------
  # Expected counts under HWE
  # ------------------------------------------------------------
  
  expected <- c(
    AA = p_hat^2,
    Aa = 2 * p_hat * q_hat,
    aa = q_hat^2
  ) * N
  
  # ------------------------------------------------------------
  # Chi-square test
  # ------------------------------------------------------------
  
  chisq <- sum((geno_counts - expected)^2 / expected)
  
  pval <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  
  # ------------------------------------------------------------
  # Data for plotting
  # ------------------------------------------------------------
  
  df <- data.frame(
    genotype = rep(required, 2),
    type = rep(c("Observed","Expected"), each = 3),
    n = c(geno_counts, expected)
  )
  
  df$genotype <- factor(df$genotype, levels = required)
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(genotype, n, fill = type)
  ) +
    
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$model,
      breaks = c("Observed","Expected"),
      labels = c("Observed","Expected (HWE)")
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(0, max(df$n) * 1.15)
    ) +
    
    ggplot2::labs(
      title = "Hardy–Weinberg equilibrium test",
      subtitle = paste0(
        "N = ", N,
        " | p = ", round(p_hat, 3),
        " | χ² = ", round(chisq, 3),
        " | P = ", signif(pval, 3)
      ),
      x = NULL,
      y = "Genotype count",
      fill = NULL
    ) +
    
    theme_qg(base_size) +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "top"
    )
}

# ============================================================
# Plot Inbreeding
# ============================================================

## Pedigree With Inbreeding

library(kinship2)

sample.ped <- data.frame(
  id     = 1:6,
  father = c(0, 0, 0, 2, 2, 5),
  mother = c(0, 0, 0, 1, 3, 4),
  sex    = c(2, 1, 2, 2, 1, 2)
)

ped <- pedigree(
  id    = sample.ped$id,
  dadid = sample.ped$father,
  momid = sample.ped$mother,
  sex   = sample.ped$sex
)

alleles <- list(
  "1" = c("A","A"),
  "2" = c("A","a"),
  "3" = c("A","A"),
  "4" = c("A","a"),
  "5" = c("A","a"),
  "6" = c("a","a")
)

plot_pedigree_alleles(ped, alleles, mode="dots")

plot_inbreeding_small_pop(
  N_values = c(10,50,100,500),
  generations = 100
)

plot_inbreeding(p=0.5, F=0.4)

plot_inbreeding_vs_drift(
  N = 50,
  generations = 100
)

plot_pedigree_alleles <- function(
    ped,
    alleles = NULL,
    mode = c("dots", "letters", "both"),
    allele_cols = qg_cols$allele,
    dx_mult = 0.25,
    dy_mult = -0.35,
    cex_dots = 1.4,
    cex_letters = 0.9,
    letters_col = "black",
    letters_font = 2,
    dot_border_col = "black",
    dot_border_lwd = 0.8,
    mar = c(4,4,4,4),
    xpd = NA,
    ...
) {
  
  mode <- match.arg(mode)
  
  if (!inherits(ped, "pedigree"))
    stop("`ped` must be a kinship2::pedigree object.")
  
  if (!is.null(alleles)) {
    
    if (!is.list(alleles))
      stop("`alleles` must be a list")
    
    bad <- names(alleles)[vapply(alleles, length, integer(1)) != 2]
    
    if (length(bad) > 0)
      stop("Each element of `alleles` must have length 2.")
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mar = mar, xpd = xpd)
  
  coords <- plot(ped, ...)
  
  if (is.null(alleles))
    return(invisible(coords))
  
  x <- coords$x
  y <- coords$y
  ids <- as.character(ped$id)
  
  dx <- coords$boxw * dx_mult
  dy <- coords$boxh * dy_mult
  
  get_pair <- function(id) {
    pair <- alleles[[id]]
    if (is.null(pair)) return(c(NA, NA))
    pair
  }
  
  for (i in seq_along(ids)) {
    
    pair <- get_pair(ids[i])
    a1 <- pair[1]
    a2 <- pair[2]
    
    if (is.na(a1) || is.na(a2)) next
    
    if (mode %in% c("dots","both")) {
      
      points(x[i]-dx, y[i]+dy,
             pch=21, bg=allele_cols[a1],
             col=dot_border_col,
             lwd=dot_border_lwd,
             cex=cex_dots)
      
      points(x[i]+dx, y[i]+dy,
             pch=21, bg=allele_cols[a2],
             col=dot_border_col,
             lwd=dot_border_lwd,
             cex=cex_dots)
    }
    
    if (mode %in% c("letters","both")) {
      
      text(x[i]-dx, y[i]+dy,
           labels=a1,
           cex=cex_letters,
           col=letters_col,
           font=letters_font)
      
      text(x[i]+dx, y[i]+dy,
           labels=a2,
           cex=cex_letters,
           col=letters_col,
           font=letters_font)
    }
  }
  
  invisible(coords)
}

plot_inbreeding_small_pop <- function(
    N_values = c(20, 50, 200),
    generations = 100,
    base_size = 18
) {
  
  df <- expand.grid(
    generation = 0:generations,
    N = N_values
  )
  
  df$F <- 1 - (1 - 1/(2*df$N))^df$generation
  df$N_label <- paste0("N = ", df$N)
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(generation, F, color = N_label)
  ) +
    
    ggplot2::geom_line(linewidth = 1.4) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = "Inbreeding increases faster in small populations",
      x = "Generation",
      y = "Inbreeding coefficient (F)",
      color = NULL
    ) +
    
    theme_qg()
}

plot_inbreeding <- function(
    p = 0.5,
    F = 0.3,
    base_size = 18
) {
  
  q <- 1 - p
  
  df <- data.frame(
    genotype = c("AA","Aa","aa"),
    freq = c(
      p^2 + p*q*F,
      2*p*q*(1-F),
      q^2 + p*q*F
    )
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(genotype, freq, fill = genotype)
  ) +
    
    ggplot2::geom_col(
      width = 0.6,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$genotype
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = paste0("Genotype frequencies with inbreeding (F = ", F, ")"),
      x = NULL,
      y = "Frequency"
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
}

plot_inbreeding_vs_drift <- function(
    N = 50,
    generations = 100,
    p0 = 0.5,
    base_size = 18
) {
  
  g <- 0:generations
  
  # Inbreeding coefficient
  F <- 1 - (1 - 1/(2*N))^g
  
  # Expected variance of allele frequency under drift
  var_p <- p0 * (1 - p0) * (1 - (1 - 1/(2*N))^g)
  
  df <- data.frame(
    generation = rep(g, 2),
    value = c(F, var_p),
    measure = factor(
      rep(c("Inbreeding (F)", "Var(p) from drift"),
          each = length(g)),
      levels = c("Inbreeding (F)", "Var(p) from drift")
    )
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(generation, value, color = measure)
  ) +
    
    ggplot2::geom_line(linewidth = 1.4) +
    
    ggplot2::scale_color_manual(
      values = c(
        "Inbreeding (F)" = qg_cols$pop["A"],
        "Var(p) from drift" = qg_cols$pop["B"]
      )
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0, max(df$value))) +
    
    ggplot2::labs(
      title = paste0("Drift and inbreeding in a population (N = ", N, ")"),
      x = "Generation",
      y = "Magnitude",
      color = NULL
    ) +
    
    theme_qg()
}

# ============================================================
# Plot Genotype and Allele Frequencies
# ============================================================


pG <- plot_genotype_freq(p = 0.7)
pA <- plot_allele_freq(p = 0.7)

pA | pG

counts <- c(AA = 10, Aa = 7, aa = 10)

pG <- plot_genotype_freq(geno_counts = counts)
pA <- plot_allele_freq(geno_counts = counts)

pA | pG

# ============================================================
# Plot Genotype Frequencies
# ============================================================

plot_genotype_freq <- function(
    p = NULL,
    geno_counts = NULL,
    stage_label = NULL,
    base_size = 18
) {
  
  required <- c("AA","Aa","aa")
  
  if (!is.null(geno_counts)) {
    
    if (is.null(names(geno_counts)))
      stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")
    
    geno_counts[setdiff(required, names(geno_counts))] <- 0
    geno_counts <- geno_counts[required]
    
    N <- sum(geno_counts)
    
    df <- data.frame(
      genotype = required,
      freq = as.numeric(geno_counts) / N
    )
    
  } else if (!is.null(p)) {
    
    q <- 1 - p
    
    df <- data.frame(
      genotype = required,
      freq = c(p^2, 2*p*q, q^2)
    )
    
  } else {
    stop("Provide either p OR geno_counts")
  }
  
  df$genotype <- factor(df$genotype, levels = required)
  
  title_text <- if (is.null(stage_label)) {
    "Genotype frequencies"
  } else {
    paste0("Genotype frequencies (", stage_label, ")")
  }
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(genotype, freq, fill = genotype)
  ) +
    
    ggplot2::geom_col(
      width = 0.6,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$genotype,
      drop = FALSE
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = title_text,
      y = "Frequency",
      x = NULL
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
}



# ============================================================
# Plot Allele Frequencies
# ============================================================

plot_allele_freq <- function(
    p = NULL,
    geno_counts = NULL,
    stage_label = NULL,
    base_size = 18
) {
  
  if (!is.null(geno_counts)) {
    
    if (is.null(names(geno_counts)))
      stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")
    
    required <- c("AA","Aa","aa")
    geno_counts[setdiff(required, names(geno_counts))] <- 0
    geno_counts <- geno_counts[required]
    
    N <- sum(geno_counts)
    
    nA <- 2*geno_counts["AA"] + geno_counts["Aa"]
    p_hat <- as.numeric(nA) / (2*N)
    
    df <- data.frame(
      allele = c("A","a"),
      freq = c(p_hat, 1 - p_hat)
    )
    
  } else if (!is.null(p)) {
    
    df <- data.frame(
      allele = c("A","a"),
      freq = c(p, 1 - p)
    )
    
  } else {
    stop("Provide either p OR geno_counts")
  }
  
  df$allele <- factor(df$allele, levels = c("A","a"))
  
  title_text <- if (is.null(stage_label)) {
    "Allele frequencies"
  } else {
    paste0("Allele frequencies (", stage_label, ")")
  }
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(allele, freq, fill = allele)
  ) +
    
    ggplot2::geom_col(
      width = 0.6,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$allele,
      drop = FALSE
    ) +
    
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    
    ggplot2::labs(
      title = title_text,
      y = "Frequency",
      x = NULL
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )
}


## Plot Allele and Genotype Frequencies from Counts
pG <- plot_genotype_freq(p = 0.7)
pA <- plot_allele_freq(p = 0.7)
pA | pG

counts <- c(AA = 10, Aa = 7, aa = 10)
pG <- plot_genotype_freq(geno_counts = counts)
pA <- plot_allele_freq(geno_counts = counts)
pA | pG


plot_genotype_freq <- function(p = NULL,
                               geno_counts = NULL,
                               stage_label = NULL,
                               base_size = 18) {

  if (!is.null(geno_counts)) {

    if (is.null(names(geno_counts)))
      stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")

    required <- c("AA","Aa","aa")
    geno_counts[setdiff(required, names(geno_counts))] <- 0
    geno_counts <- geno_counts[required]

    N <- sum(geno_counts)

    df <- data.frame(
      genotype = required,
      freq = as.numeric(geno_counts) / N
    )

  } else if (!is.null(p)) {

    q <- 1 - p

    df <- data.frame(
      genotype = c("AA","Aa","aa"),
      freq = c(p^2, 2*p*q, q^2)
    )

  } else {
    stop("Provide either p OR geno_counts")
  }

  # Title logic
  title_text <- if (is.null(stage_label)) {
    "Genotype frequencies"
  } else {
    paste0("Genotype frequencies (", stage_label, ")")
  }

  ggplot(df, aes(genotype, freq, fill = genotype)) +
    geom_col(width = 0.6, color = "black") +
    scale_fill_manual(values = bio_cols$genotype) +
    coord_cartesian(ylim = c(0,1)) +
    labs(
      title = title_text,
      y = "Frequency",
      x = NULL
    ) +
    theme_bio(base_size)
}

plot_allele_freq <- function(p = NULL,
                             geno_counts = NULL,
                             stage_label = NULL,
                             base_size = 18) {

  if (!is.null(geno_counts)) {

    if (is.null(names(geno_counts)))
      stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")

    required <- c("AA","Aa","aa")
    geno_counts[setdiff(required, names(geno_counts))] <- 0
    geno_counts <- geno_counts[required]

    N <- sum(geno_counts)

    nA <- 2*geno_counts["AA"] + geno_counts["Aa"]
    p_hat <- as.numeric(nA) / (2*N)
    q_hat <- 1 - p_hat

    df <- data.frame(
      allele = c("A","a"),
      freq = c(p_hat, q_hat)
    )

  } else if (!is.null(p)) {

    df <- data.frame(
      allele = c("A","a"),
      freq = c(p, 1-p)
    )

  } else {
    stop("Provide either p OR geno_counts")
  }

  title_text <- if (is.null(stage_label)) {
    "Allele frequencies"
  } else {
    paste0("Allele frequencies (", stage_label, ")")
  }

  ggplot(df, aes(allele, freq, fill = allele)) +
    geom_col(width = 0.6, color = "black") +
    scale_fill_manual(values = bio_cols$allele) +
    coord_cartesian(ylim = c(0,1)) +
    labs(
      title = title_text,
      y = "Frequency",
      x = NULL
    ) +
    theme_bio(base_size)
}




# ============================================================
# Gene Pool Concept (qgplots)
# ============================================================

pop <- generate_population(N = 16, p = 0.7, seed = 1)

geno <- layout_population_grid(pop$geno)

p_ind <- plot_population_individuals(geno, pop$geno_labels)
p_cnt <- plot_population_counts(pop$geno_labels)

p_ind | p_cnt

# ------------------------------------------------------------
# Generate population
# ------------------------------------------------------------

generate_population <- function(
    N = NULL,
    p = NULL,
    geno_counts = NULL,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ------------------------------------------------------------
  # MODE 1: explicit genotype counts
  # ------------------------------------------------------------
  
  if (!is.null(geno_counts)) {
    
    stopifnot(all(c("AA","Aa","aa") %in% names(geno_counts)))
    
    genotypes <- rep(names(geno_counts), geno_counts)
    genotypes <- sample(genotypes)
    
    N <- length(genotypes)
    
  } else {
    
    stopifnot(!is.null(N), !is.null(p))
    
    n_alleles <- 2 * N
    
    alleles <- sample(c(
      rep("A", round(p * n_alleles)),
      rep("a", n_alleles - round(p * n_alleles))
    ))
    
    genotypes <- paste0(
      alleles[seq(1, n_alleles, 2)],
      alleles[seq(2, n_alleles, 2)]
    )
    
    genotypes[genotypes %in% c("aA")] <- "Aa"
  }
  
  geno_labels <- data.frame(
    id = seq_len(N),
    genotype = genotypes,
    stringsAsFactors = FALSE
  )
  
  # Expand to allele copies
  allele_vec <- unlist(strsplit(genotypes, ""))
  
  geno <- data.frame(
    id = rep(seq_len(N), each = 2),
    copy = rep(1:2, N),
    allele = allele_vec,
    stringsAsFactors = FALSE
  )
  
  list(
    geno = geno,
    geno_labels = geno_labels
  )
}


# ------------------------------------------------------------
# Arrange individuals on grid
# ------------------------------------------------------------

layout_population_grid <- function(geno) {
  
  N <- length(unique(geno$id))
  grid_size <- ceiling(sqrt(N))
  
  geno$row <- ceiling(geno$id / grid_size)
  geno$col <- (geno$id - 1) %% grid_size + 1
  
  geno$x <- geno$col + ifelse(geno$copy == 1, -0.12, 0.12)
  geno$y <- -geno$row
  
  geno
}


# ------------------------------------------------------------
# Plot individuals (diploid gene pool)
# ------------------------------------------------------------

plot_population_individuals <- function(
    geno,
    geno_labels,
    base_size = 18
) {
  
  N <- nrow(geno_labels)
  grid_size <- ceiling(sqrt(N))
  
  label_df <- data.frame(
    id = geno_labels$id,
    genotype = geno_labels$genotype
  )
  
  label_df$row <- ceiling(label_df$id / grid_size)
  label_df$col <- (label_df$id - 1) %% grid_size + 1
  label_df$x <- label_df$col
  label_df$y <- -label_df$row + 0.35
  
  ggplot2::ggplot() +
    
    ggplot2::geom_point(
      data = geno,
      ggplot2::aes(x, y, fill = allele),
      shape = 21,
      size = 6,
      color = "black",
      stroke = 0.7
    ) +
    
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x, y, label = genotype),
      size = 4,
      color = qg_cols$grey["dark"]
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$allele
    ) +
    
    ggplot2::coord_equal() +
    
    ggplot2::labs(
      title = "Diploid individuals"
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
      
}


# ------------------------------------------------------------
# Plot genotype counts
# ------------------------------------------------------------

plot_population_counts <- function(
    geno_labels,
    base_size = 18
) {
  
  geno_counts <- as.data.frame(
    table(geno_labels$genotype)
  )
  
  colnames(geno_counts) <- c("genotype","n")
  
  geno_counts$genotype <- factor(
    geno_counts$genotype,
    levels = c("AA","Aa","aa")
  )
  
  ggplot2::ggplot(
    geno_counts,
    ggplot2::aes(genotype, n, fill = genotype)
  ) +
    
    ggplot2::geom_col(
      width = 0.6,
      color = "black"
    ) +
    
    ggplot2::scale_fill_manual(
      values = qg_cols$genotype,
      drop = FALSE
    ) +
    
    ggplot2::labs(
      title = "Genotype counts",
      x = NULL,
      y = "Count"
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      legend.position = "none"
    )
}



# ============================================================
# Sequence / Haplotype Visualization (qgplots)
# ============================================================


## Genetic Variation (SNPs, Microsatellites, InDels)
seqs <- c(
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGTATCGCTCC--TACTTAGAGAG---TTAGTC",
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",
  "GGTATCGCTCC--TACTTAGAGAG---CTTAGT"
)

seqs <- c(
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGTATCGCTCC--TACTTAGAGAG---TTAGTC",
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",
  "GGTATCGCTCC--TACTTAGAGAG---CTTAGT"
)

plot_sequence(
  seqs,
  highlight_regions = list(c(12,15), c(25,27))
)


plot_sequence <- function(
    seqs,
    highlight_regions = NULL,
    show_seg_sites = TRUE,
    base_size = 16,
    title = "Sequence variation across homologous chromosomes"
) {
  
  # ------------------------------------------------------------
  # Checks
  # ------------------------------------------------------------
  
  stopifnot(is.character(seqs))
  stopifnot(length(unique(nchar(seqs))) == 1)
  
  n_seq <- length(seqs)
  seq_length <- nchar(seqs[1])
  
  # ------------------------------------------------------------
  # Convert to data frame (base R only)
  # ------------------------------------------------------------
  
  seq_split <- strsplit(seqs, "")
  
  df <- data.frame(
    id  = rep(seq_along(seqs), each = seq_length),
    pos = rep(seq_len(seq_length), times = n_seq),
    base = unlist(seq_split),
    stringsAsFactors = FALSE
  )
  
  # ------------------------------------------------------------
  # Detect segregating sites
  # ------------------------------------------------------------
  
  seg_sites <- tapply(df$base, df$pos, function(x) length(unique(x)) > 1)
  seg_pos <- as.integer(names(seg_sites[seg_sites]))
  
  df$is_seg <- df$pos %in% seg_pos
  
  # Reverse order for plotting (top sequence first)
  df$id <- factor(df$id, levels = rev(unique(df$id)))
  
  # ------------------------------------------------------------
  # Base colors (consistent + safe)
  # ------------------------------------------------------------
  
  base_cols <- c(
    A = "#D55E00",
    T = "#009E73",
    C = "#0072B2",
    G = "#CC79A7",
    "-" = qg_cols$grey["mid"]
  )
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = pos, y = id)
  ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        label = base,
        color = base,
        fontface = if (show_seg_sites) {
          ifelse(is_seg, "bold", "plain")
        } else {
          "plain"
        }
      ),
      family = "mono",
      size = 4.5
    ) +
    
    ggplot2::scale_color_manual(
      values = base_cols
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = seq(0, seq_length, 5),
      limits = c(0.5, seq_length + 0.5),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    
    ggplot2::coord_fixed(ratio = 0.6, clip = "off") +
    
    ggplot2::labs(
      x = "Genomic position",
      y = NULL,
      title = title
    ) +
    
    theme_qg() +
    
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = base_size * 0.7)
    )
  
  # ------------------------------------------------------------
  # Highlight regions
  # ------------------------------------------------------------
  
  if (!is.null(highlight_regions)) {
    
    pad <- 0.2
    
    ymin = 0.5 - pad
    ymax = n_seq + 0.5 + pad
    
    for (region in highlight_regions) {
      
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = region[1] - 0.5,
          xmax = region[2] + 0.5,
          ymin = ymin,
          ymax = ymax,
          fill = NA,
          color = qg_cols$grey["dark"],
          linewidth = 0.7
        )
    }
  }
  
  return(p)
}

# ============================================================
# LD Decay Under Recombination (qgplots)
# ============================================================

plot_ld_decay(
  D0 = 0.1,
  r = 0.01,
  generations = 50,
  title = "Slow LD decay (tight linkage)"
)

plot_ld_decay(
  D0 = 0.1,
  r = 0.1,
  generations = 50,
  title = "Moderate LD decay"
)

plot_ld_decay(
  D0 = 0.1,
  r = 0.5,
  generations = 20,
  title = "Rapid LD decay (unlinked loci)"
)

# ============================================================


plot_ld_decay <- function(
    D0,
    r,
    generations = 20,
    pA = 0.5,
    pB = 0.5,
    show_r2 = TRUE,
    title = "Decay of Linkage Disequilibrium"
) {
  
  # ------------------------------------------------------------
  # Simulate LD decay
  # ------------------------------------------------------------
  D <- numeric(generations + 1)
  D[1] <- D0
  
  for (t in seq_len(generations)) {
    D[t + 1] <- (1 - r) * D[t]
  }
  
  # ------------------------------------------------------------
  # If showing r²
  # ------------------------------------------------------------
  if (show_r2) {
    
    r2 <- D^2 / (pA * (1 - pA) * pB * (1 - pB))
    
    df_long <- data.frame(
      generation = rep(0:generations, 2),
      value = c(D, r2),
      measure = rep(c("D", "r2"), each = generations + 1)
    )
    
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(
        x = generation,
        y = value
      )
    ) +
      
      # Plot D
      ggplot2::geom_line(
        data = subset(df_long, measure == "D"),
        linewidth = 1.5,
        color = qg_cols$pop["A"]
      ) +
      
      # Plot r2
      ggplot2::geom_line(
        data = subset(df_long, measure == "r2"),
        linewidth = 1.5,
        color = qg_cols$pop["B"]
      ) +
      
      ggplot2::labs(
        title = title,
        subtitle = bquote(D[0] == .(D0) ~ "|" ~ r == .(r)),
        x = "Generation",
        y = "LD measure"
      ) +
      
      theme_qg()
  } else {
    
    df <- data.frame(
      generation = 0:generations,
      D = D
    )
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = generation, y = D)
    ) +
      ggplot2::geom_line(
        linewidth = 1.5,
        color = qg_cols$pop["A"]
      ) +
      ggplot2::labs(
        title = title,
        subtitle = bquote(
          D[0] == .(D0) ~ "|" ~ r == .(r)
        ),
        x = "Generation",
        y = "D"
      ) +
      theme_qg()
  }
  
  p
}

# ============================================================
# Linkage Disequilibrium Visualization (qgplots)
# ============================================================


plot_linkage_disequilibrium(
  pAB = 0.40,
  pAb = 0.10,
  paB = 0.10,
  pab = 0.40,
  title = "Strong LD"
)  

plot_linkage_disequilibrium(
  pAB = 0.25,
  pAb = 0.25,
  paB = 0.25,
  pab = 0.25,
  title = "Linkage Equilibrium"
)

plot_linkage_disequilibrium(
  pAB = 0.30,
  pAb = 0.20,
  paB = 0.20,
  pab = 0.30,
  title = "Moderate LD"
)


# ============================================================

plot_linkage_disequilibrium <- function(
    pAB,
    pAb,
    paB,
    pab,
    title = "Linkage Disequilibrium"
) {
  
  # ------------------------------------------------------------
  # Validate
  # ------------------------------------------------------------
  total <- pAB + pAb + paB + pab
  if (abs(total - 1) > 1e-6) {
    stop("Haplotype frequencies must sum to 1.")
  }
  
  # ------------------------------------------------------------
  # Allele frequencies
  # ------------------------------------------------------------
  pA <- pAB + pAb
  pB <- pAB + paB
  
  # ------------------------------------------------------------
  # LD statistics
  # ------------------------------------------------------------
  D  <- pAB - pA * pB
  r2 <- D^2 / (pA * (1 - pA) * pB * (1 - pB))
  
  # D'
  if (D > 0) {
    Dmax <- min(pA * (1 - pB), (1 - pA) * pB)
  } else {
    Dmax <- min(pA * pB, (1 - pA) * (1 - pB))
  }
  
  Dprime <- ifelse(Dmax > 0, D / Dmax, 0)
  
  # ------------------------------------------------------------
  # Observed vs Expected
  # ------------------------------------------------------------
  haplotypes <- c("AB", "Ab", "aB", "ab")
  
  observed <- c(pAB, pAb, paB, pab)
  
  expected <- c(
    pA * pB,
    pA * (1 - pB),
    (1 - pA) * pB,
    (1 - pA) * (1 - pB)
  )
  
  df <- data.frame(
    haplotype = rep(haplotypes, 2),
    freq = c(observed, expected),
    type = rep(c("Observed", "Expected (Independence)"),
               each = 4)
  )
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = haplotype,
                 y = freq,
                 fill = type)
  ) +
    ggplot2::geom_col(position = "dodge",
                      color = "black") +
    ggplot2::scale_fill_manual(values = qg_cols$model) +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "D = ", round(D, 3),
        " | D' = ", round(Dprime, 3),
        " | r² = ", round(r2, 3)
      ),
      x = "Haplotype",
      y = "Frequency"
    ) +
    theme_qg()
}





# ============================================================
# Chromosome & Recombination Visualizations
# ============================================================


# ------------------------------------------------------------
# Example 1: Basic two-locus recombination
# ------------------------------------------------------------

plot_recombination_two_loci()

# ------------------------------------------------------------
# Example 2: Parental haplotypes
# ------------------------------------------------------------

plot_haplotype_inheritance()

# ------------------------------------------------------------
# Example 3: Custom recombinant gametes
# ------------------------------------------------------------

ggplot2::ggplot() +
  
  # Parental haplotypes
  draw_haplotype(3,
                 c("A","B"),
                 qg_cols$pop["A"]) +
  
  draw_haplotype(2.5,
                 c("a","b"),
                 qg_cols$pop["A"]) +
  
  # Recombinant haplotypes
  draw_haplotype(1.2,
                 c("A","b"),
                 qg_cols$pop["A"],
                 qg_cols$pop["B"]) +

# ------------------------------------------------------------
# Example 4: Recombination vs inheritance
# ------------------------------------------------------------

p1 <- plot_recombination_two_loci()
p2 <- plot_haplotype_inheritance()

patchwork::wrap_plots(p1, p2)

  draw_haplotype(0.7,
                 c("a","B"),
                 qg_cols$pop["B"],
                 qg_cols$pop["A"]) +
  
  ggplot2::coord_cartesian(xlim = c(0,10), ylim = c(0,4)) +
  
  theme_qg_schematic()

  
# ============================================================

theme_qg_schematic <- function() {
  theme_qg() +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

draw_haplotype <- function(
    y,
    alleles,
    color_left,
    color_right = NULL,
    xA = 3,
    xB = 7,
    x_cross = 5,
    allele_size = 5,
    locus_height = 0.15
) {
  
  layers <- list()
  
  if (is.null(color_right)) {
    
    layers <- c(layers, list(
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = 10, y = y, yend = y),
        linewidth = 3,
        color = color_left
      )
    ))
    
  } else {
    
    layers <- c(layers, list(
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = x_cross, y = y, yend = y),
        linewidth = 3,
        color = color_left
      ),
      ggplot2::geom_segment(
        ggplot2::aes(x = x_cross, xend = 10, y = y, yend = y),
        linewidth = 3,
        color = color_right
      )
    ))
  }
  
  layers <- c(layers, list(
    ggplot2::geom_segment(
      ggplot2::aes(x = xA, xend = xA,
                   y = y - locus_height,
                   yend = y + locus_height),
      linewidth = 1
    ),
    ggplot2::geom_segment(
      ggplot2::aes(x = xB, xend = xB,
                   y = y - locus_height,
                   yend = y + locus_height),
      linewidth = 1
    ),
    ggplot2::annotate("text",
                      x = xA, y = y + 0.3,
                      label = alleles[1],
                      size = allele_size),
    ggplot2::annotate("text",
                      x = xB, y = y + 0.3,
                      label = alleles[2],
                      size = allele_size)
  ))
  
  layers
}

plot_recombination_two_loci <- function() {
  
  ggplot2::ggplot() +
    
    draw_haplotype(3, c("A","B"), qg_cols$pop["A"]) +
    draw_haplotype(2, c("a","b"), qg_cols$pop["B"]) +
    
    draw_haplotype(0, c("A","b"),
                   qg_cols$pop["A"],
                   qg_cols$pop["B"]) +
    
    draw_haplotype(-1, c("a","B"),
                   qg_cols$pop["B"],
                   qg_cols$pop["A"]) +
    
    ggplot2::geom_segment(
      ggplot2::aes(x = 5, xend = 5, y = 1.6, yend = 3.4),
      linetype = "dashed"
    ) +
    
    ggplot2::annotate("text",
                      x = 5, y = 3.9,
                      label = "Parental haplotypes") +
    
    ggplot2::annotate("text",
                      x = 5, y = 1,
                      label = "Recombinant haplotypes") +
    
    ggplot2::coord_cartesian(
      xlim = c(-1, 11),
      ylim = c(-2, 5)
    ) +
    
    theme_qg_schematic()
}

plot_haplotype_inheritance <- function() {
  
  ggplot2::ggplot() +
    
    ggplot2::annotate("text",
                      x = 5, y = 5.2,
                      label = "Mother",
                      fontface = "bold") +
    
    draw_haplotype(4.5, c("A","B"), qg_cols$pop["A"]) +
    draw_haplotype(4.0, c("a","b"), qg_cols$pop["A"]) +
    
    ggplot2::annotate("text",
                      x = 5, y = 3.3,
                      label = "Father",
                      fontface = "bold") +
    
    draw_haplotype(2.8, c("A","b"), qg_cols$pop["B"]) +
    draw_haplotype(2.3, c("a","B"), qg_cols$pop["B"]) +
    
    ggplot2::coord_cartesian(
      xlim = c(-1, 11),
      ylim = c(0, 6)
    ) +
    
    theme_qg_schematic()
}



# ============================================================
# Selection Dynamics (qgplots)
# ============================================================

p_dir_traj <- plot_selection_trajectory(
  p0 = 0.7, T = 30,
  wAA = 1, wAa = 0.9, waa = 0.6,
  title = "Directional selection"
)

p_het_traj <- plot_selection_trajectory(
  p0 = 0.7, T = 30,
  wAA = 0.8, wAa = 1.0, waa = 0.8,
  title = "Heterozygote advantage"
)

patchwork::wrap_plots(p_dir_traj, p_het_traj)

# ------------------------------------------------------------
# Internal: One selection step
# ------------------------------------------------------------

step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2 * wAA + 2 * p * q * wAa + q^2 * waa
  (p^2 * wAA + p * q * wAa) / wbar
}


# ------------------------------------------------------------
# Internal: Simulate selection trajectory
# ------------------------------------------------------------

simulate_selection <- function(
    p0,
    T,
    wAA,
    wAa,
    waa
) {
  
  p <- numeric(T + 1)
  p[1] <- p0
  
  for (t in seq_len(T)) {
    p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
  }
  
  data.frame(
    generation = 0:T,
    frequency  = p
  )
}


# ------------------------------------------------------------
# Selection trajectory plot
# ------------------------------------------------------------

plot_selection_trajectory <- function(
    p0,
    T,
    wAA,
    wAa,
    waa,
    title = "Selection"
) {
  
  df <- simulate_selection(p0, T, wAA, wAa, waa)
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = generation, y = frequency)
  ) +
    ggplot2::geom_line(
      linewidth = 1.5,
      color = qg_cols$pop["A"]
    ) +
    ggplot2::geom_hline(
      yintercept = p0,
      linetype = "dashed",
      color = qg_cols$grey["mid"]
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = title,
      x = "Generation",
      y = "Allele frequency (p)"
    ) +
    theme_qg()
}


# ============================================================
# Wright–Fisher Drift Simulation (qgplots)
# ============================================================

p_small <- plot_drift(N = 10,  p0 = 0.05, generations = 100)
p_large <- plot_drift(N = 200, p0 = 0.05, generations = 100)

patchwork::wrap_plots(p_small, p_large)

# ------------------------------------------------------------
# Internal simulation function
# ------------------------------------------------------------

simulate_drift <- function(
    N,
    generations = 100,
    replicates = 20,
    p0 = 0.7
) {
  
  results <- vector("list", replicates)
  
  for (r in seq_len(replicates)) {
    
    p <- numeric(generations + 1)
    p[1] <- p0
    
    for (t in seq_len(generations)) {
      p[t + 1] <- stats::rbinom(1, 2 * N, p[t]) / (2 * N)
    }
    
    results[[r]] <- data.frame(
      generation = 0:generations,
      frequency  = p,
      replicate  = r
    )
  }
  
  dplyr::bind_rows(results)
}


# ------------------------------------------------------------
# Drift Plot
# ------------------------------------------------------------

plot_drift <- function(
    N,
    generations = 100,
    replicates = 20,
    p0 = 0.7,
    show_mean = TRUE
) {
  
  df <- simulate_drift(N, generations, replicates, p0)
  
  g <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = generation,
      y = frequency,
      group = replicate
    )
  ) +
    ggplot2::geom_line(
      alpha = 0.4,
      color = qg_cols$grey["dark"]
    ) +
    ggplot2::geom_hline(
      yintercept = p0,
      linetype = "dashed",
      color = qg_cols$grey["mid"]
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0("Genetic Drift (N = ", N, ")"),
      x = "Generation",
      y = "Allele frequency (p)"
    ) +
    theme_qg()
  
  if (show_mean) {
    g <- g +
      ggplot2::stat_summary(
        ggplot2::aes(group = 1),
        fun = mean,
        geom = "line",
        linewidth = 1.3,
        color = qg_cols$pop["A"]
      )
  }
  
  g
}







# ============================================================
# qgplots — Admixture Examples
# Conceptual → Genomic → Statistical Consequence
# ============================================================


# ------------------------------------------------------------
# 1) Conceptual Admixture (Bucket Illustration)
# ------------------------------------------------------------
# Two source populations contribute to an admixed population.
# The number of colored balls represents ancestry proportions.

plot_admixture_buckets(
  n_each = 12,
  mix_red = 9,
  mix_blue = 3,
  seed = 1,
  title = "Conceptual illustration of admixture (75% from Population A)"
)



# ------------------------------------------------------------
# 2) Genome-wide Ancestry Representation (ADMIXTURE-style)
# ------------------------------------------------------------
# Individuals are represented by stacked ancestry proportions.
# The first three groups are nearly pure.
# The fourth group is admixed.

plot_admixture_barplot(
  group_sizes = c(15, 15, 15, 15),
  alpha_list = list(
    c(20,1,1),
    c(1,20,1),
    c(1,1,20),
    c(4,3,3)
  ),
  group_names = c("Population A",
                  "Population B",
                  "Population C",
                  "Admixed"),
  seed = 123,
  title = "ADMIXTURE-style ancestry plot (K = 3)"
)


# ------------------------------------------------------------
# 3) Statistical Consequence: Wahlund Effect
# ------------------------------------------------------------
# Two subpopulations are each in Hardy–Weinberg equilibrium.
# When pooled, the population deviates from HWE.
# This produces heterozygote deficiency (Wahlund effect).

plot_admixture_wahlund(
  p1 = 0.2,   # Allele frequency in Population 1
  p2 = 0.8,   # Allele frequency in Population 2
  m  = 0.5    # Equal mixing
)



# ------------------------------------------------------------
# Alternative: Using observed genotype counts instead
# ------------------------------------------------------------
# Allele frequencies are estimated from data,
# then pooled and tested against HWE expectation.

plot_admixture_wahlund(
  geno1 = c(AA = 40, Aa = 20, aa = 40),
  geno2 = c(AA = 10, Aa = 30, aa = 60)
)

# ============================================================
# qgplots: Global Colors & Theme
# ============================================================

#' Global color palette for qgplots
#'
#' Centralized color definitions used across all plots.
#' Modify here to update the entire visual identity.
#'
#' @export
qg_cols <- list(
  
  # Core population colors
  pop = c(
    A = "#D94848",   # red
    B = "#2B6CB0",   # blue
    C = "#2F9E44",   # green
    D = "#7B2CBF"    # purple
  ),
  
  # Alleles
  allele = c(
    A = "#D94848",
    a = "#2B6CB0"
  ),
  
  # Genotypes
  genotype = c(
    AA = "#D94848",
    Aa = "#9E9AC8",
    aa = "#2B6CB0"
  ),
  
  # Model comparison (Observed vs Expected)
  model = c(
    Observed = "#1B9E77",  # green
    Expected = "#D95F02"   # orange
  ),
  
  # Neutral greys
  grey = c(
    light = "#F5F5F5",
    mid   = "#CCCCCC",
    dark  = "#333333"
  )
  
)

#' qgplots default theme
#'
#' Clean, publication-ready theme for teaching genetics.
#'
#' @param base_size Base font size
#' @export
theme_qg <- function(base_size = 14) {
  
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      
      # Titles
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      
      # Axis titles
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      
      # Remove clutter
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      
      # Legend
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      
      # Clean background
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )
}

#' Default discrete scale for qgplots
#' @export
scale_fill_qg <- function(...) {
  ggplot2::scale_fill_manual(values = qg_cols$pop, ...)
}
# ------------------------------------------------------------
# Internal helper: bucket point generator
# ------------------------------------------------------------

.make_bucket_points <- function(n, bucket, color, x0, x1, y0, y1) {
  data.frame(
    bucket = bucket,
    color  = color,
    x = stats::runif(n, x0 + 0.08, x1 - 0.08),
    y = stats::runif(n, y0 + 0.10, y1 - 0.10)
  )
}


# ============================================================
# 1. Simple bucket illustration of admixture
# ============================================================

plot_admixture_buckets <- function(
    n_each = 10,
    mix_red = 7,
    mix_blue = 3,
    seed = NULL,
    title = "Illustration of admixture"
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  buckets <- data.frame(
    bucket = c("Population A", "Population B", "Admixed"),
    xmin   = c(0, 1.3, 2.6),
    xmax   = c(1.0, 2.3, 3.6),
    ymin   = 0,
    ymax   = 1.2
  )
  
  pts <- rbind(
    .make_bucket_points(n_each, "Population A", "A",
                        buckets$xmin[1], buckets$xmax[1], 0, 1.2),
    .make_bucket_points(n_each, "Population B", "B",
                        buckets$xmin[2], buckets$xmax[2], 0, 1.2),
    .make_bucket_points(mix_red, "Admixed", "A",
                        buckets$xmin[3], buckets$xmax[3], 0, 1.2),
    .make_bucket_points(mix_blue, "Admixed", "B",
                        buckets$xmin[3], buckets$xmax[3], 0, 1.2)
  )
  
  # --- Robust palette extraction (prevents "No shared levels" warning) ---
  cols <- qg_cols$pop[c("A", "B")]
  if (anyNA(cols)) {
    stop("qg_cols$pop must contain named entries 'A' and 'B'.")
  }
  names(cols) <- c("A", "B")  # force correct matching to pts$color values
  
  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = buckets,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = qg_cols$grey["light"],
      color = qg_cols$grey["dark"]
    ) +
    ggplot2::geom_text(
      data = buckets,
      ggplot2::aes(x = (xmin + xmax)/2, y = ymax + 0.12, label = bucket),
      fontface = "bold"
    ) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = x, y = y, color = color),
      size = 4
    ) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::coord_equal(xlim = c(-0.1, 3.7),
                         ylim = c(-0.1, 1.45)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_qg() +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}



# ============================================================
# 2. ADMIXTURE-style ancestry barplot
# ============================================================

plot_admixture_barplot <- function(
    group_sizes,
    alpha_list,
    group_names = NULL,
    colors = NULL,
    seed = NULL,
    title = NULL,
    show_ids = FALSE,
    rotate_ids = FALSE
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  K <- length(alpha_list[[1]])
  
  if (is.null(group_names)) {
    group_names <- paste0("Group ", seq_along(group_sizes))
  }
  
  # ----------------------------------------------------------
  # Colors
  # ----------------------------------------------------------
  if (is.null(colors)) {
    base_cols <- unname(qg_cols$pop)
    colors <- base_cols[seq_len(K)]
  }
  names(colors) <- paste0("Cluster", seq_len(K))
  
  # ----------------------------------------------------------
  # Simulate ancestry proportions
  # ----------------------------------------------------------
  props_list <- lapply(seq_along(group_sizes), function(i) {
    MCMCpack::rdirichlet(group_sizes[i], alpha_list[[i]])
  })
  
  props <- do.call(rbind, props_list)
  
  df <- data.frame(
    individual = 1:nrow(props),
    group = rep(group_names, group_sizes),
    props
  )
  
  colnames(df)[3:(K+2)] <- paste0("Cluster", seq_len(K))
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = starts_with("Cluster"),
    names_to = "ancestry",
    values_to = "proportion"
  )
  
  # ----------------------------------------------------------
  # Axis control
  # ----------------------------------------------------------
  axis_text_x <- if (show_ids) {
    if (rotate_ids) {
      ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    } else {
      ggplot2::element_text()
    }
  } else {
    ggplot2::element_blank()
  }
  
  axis_ticks_x <- if (show_ids) {
    ggplot2::element_line()
  } else {
    ggplot2::element_blank()
  }
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = factor(individual),
                               y = proportion,
                               fill = ancestry)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::labs(x = "Individuals",
                  y = "Ancestry proportion",
                  title = title,
                  fill = "Cluster") +
    theme_qg() +
    ggplot2::theme(
      axis.text.x = axis_text_x,
      axis.ticks.x = axis_ticks_x
    )
}

# ============================================================
# 3. Wahlund effect (HWE deviation due to admixture)
# ============================================================

plot_admixture_wahlund <- function(
    p1 = NULL,
    p2 = NULL,
    geno1 = NULL,
    geno2 = NULL,
    m = 0.5,
    N_theoretical = 1000
) {
  
  required <- c("AA","Aa","aa")
  
  clean_counts <- function(g) {
    g[setdiff(required, names(g))] <- 0
    g[required]
  }
  
  get_p <- function(g) {
    N <- sum(g)
    nA <- 2*g["AA"] + g["Aa"]
    as.numeric(nA)/(2*N)
  }
  
  geno_hwe <- function(p){
    q <- 1 - p
    c(AA = p^2, Aa = 2*p*q, aa = q^2)
  }
  
  if (!is.null(geno1) && !is.null(geno2)) {
    
    geno1 <- clean_counts(geno1)
    geno2 <- clean_counts(geno2)
    
    geno_pool <- geno1 + geno2
    N <- sum(geno_pool)
    
    p1 <- get_p(geno1)
    p2 <- get_p(geno2)
    pbar <- get_p(geno_pool)
    
    g_pool_obs <- geno_pool / N
    g_pool_hwe <- geno_hwe(pbar)
    
    expected_counts <- g_pool_hwe * N
    chisq <- sum((geno_pool - expected_counts)^2 / expected_counts)
    pval  <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
    
    subtitle_text <- paste0(
      "N = ", N,
      " | p = ", round(pbar,3),
      " | χ² = ", round(chisq,2),
      " | P = ", signif(pval,3)
    )
    
  } else if (!is.null(p1) && !is.null(p2)) {
    
    pbar <- m*p1 + (1-m)*p2
    
    g1 <- geno_hwe(p1)
    g2 <- geno_hwe(p2)
    
    g_pool_obs <- m*g1 + (1-m)*g2
    g_pool_hwe <- geno_hwe(pbar)
    
    N <- N_theoretical
    
    obs_counts <- g_pool_obs * N
    expected_counts <- g_pool_hwe * N
    
    chisq <- sum((obs_counts - expected_counts)^2 / expected_counts)
    pval  <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
    
    subtitle_text <- paste0(
      "Theoretical N = ", N,
      " | p = ", round(pbar,3),
      " | χ² = ", round(chisq,2),
      " | P = ", signif(pval,3)
    )
    
  } else {
    stop("Provide either (p1,p2) OR (geno1,geno2)")
  }
  
  allele_df <- data.frame(
    pop = c("Pop 1","Pop 2","Pooled"),
    A = c(p1,p2,pbar),
    a = 1 - c(p1,p2,pbar)
  )
  
  allele_long <- tidyr::pivot_longer(
    allele_df,
    cols = c("A","a"),
    names_to = "allele",
    values_to = "freq"
  )
  
  geno_df <- data.frame(
    genotype = required,
    Observed = as.numeric(g_pool_obs),
    Expected = as.numeric(g_pool_hwe)
  )
  
  geno_long <- tidyr::pivot_longer(
    geno_df,
    cols = -genotype,
    names_to = "type",
    values_to = "freq"
  )
  
  pA <- ggplot2::ggplot(allele_long,
                        ggplot2::aes(pop, freq, fill = allele)) +
    ggplot2::geom_col(position = "dodge", color="black") +
    ggplot2::scale_fill_manual(values = qg_cols$allele) +
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::labs(title="Allele frequencies",
                  y="Frequency", x=NULL) +
    theme_qg()
  
  pG <- ggplot2::ggplot(geno_long,
                        ggplot2::aes(genotype, freq, fill = type)) +
    ggplot2::geom_col(position = "dodge", color="black") +
    ggplot2::scale_fill_manual(values = qg_cols$model) +
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::labs(
      title = "Pooled genotypes deviate from HWE",
      subtitle = subtitle_text,
      y = "Frequency", x = NULL
    ) +
    theme_qg()
  
  patchwork::wrap_plots(pA, pG, ncol = 2) +
    patchwork::plot_annotation(
      title = "Admixture generates Hardy–Weinberg deviation (Wahlund effect)"
    )
}

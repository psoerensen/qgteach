# ============================================================
# POPULATION GENETICS TEACHING PLOTS USING NEW FRAMEWORK
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


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

# 1) Colored dots
plot_pedigree_alleles(ped, alleles, mode = "dots")


## Inbreeding and Population Size

plot_inbreeding_small_pop(
  N_values = c(10, 50, 100, 500),
  generations = 100
)


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

plot_sequence(seqs,highlight_regions = list(c(12,15), c(25,27)))


## Gene Pool Concept

pop <- generate_population(N=16, p=0.7, seed=1)
pop <- generate_population(geno_counts = c(AA = 5, Aa = 8, aa = 3))

geno <- layout_population_grid(pop$geno)

p_ind  <- plot_population_individuals(geno, pop$geno_labels, base_size)
p_cnt  <- plot_population_counts(pop$geno_labels, base_size)

p_ind | p_cnt 


## Genetic Drift(N small vs large)

p_small <- plot_drift(N = 10, p0=0.05, generations=100)
p_large <- plot_drift(N = 200, p0=0.05, generations=100)

p_small | p_large


## Selection Trajectories

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

p_dir_traj | p_het_traj


## Plot Allele and Genotype Frequencies from Counts
pG <- plot_genotype_freq(p = 0.7)
pA <- plot_allele_freq(p = 0.7)
pA | pG

counts <- c(AA = 10, Aa = 7, aa = 10)
pG <- plot_genotype_freq(geno_counts = counts)
pA <- plot_allele_freq(geno_counts = counts)
pA | pG


## Hardy-Weinberg
plot_hwe_comparison(c(AA = 10, Aa = 7, aa = 50))


## Inbreeding
plot_inbreeding(p = 0.5, F = 0.4)


## Admixture (Wahlund Effect)
plot_admixture_wahlund(p1 = 0.2, p2 = 0.8, m=0.5)
plot_admixture_wahlund(
  geno1 = c(AA=40, Aa=20, aa=40),
  geno2 = c(AA=10, Aa=30, aa=60)
)

## Allele and Genotype Frequencies
pA <- plot_allele_freq(p = 0.7)
pG <- plot_genotype_freq(p = 0.7)

pA | pG


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


plot_heterozygosity_decay(
  N = 20,
  generations = 100,
  H0 = 0.5
)

plot_heterozygosity_decay(
  N = 200,
  generations = 100,
  H0 = 0.5
)

p_small <- plot_heterozygosity_decay(N = 20, generations = 100)
p_large <- plot_heterozygosity_decay(N = 200, generations = 100)

p_small | p_large

## Haplotype sweep plot
psweep <- make_sweep_plots(
  mode = "Realistic Hard Sweep",
  n_hap = 100,
  L = 201,
  sweep_freq = 0.9,
  width = 30,
  sort_haps = TRUE,
  seed = 1
)

psweep$hapPlot
psweep$gdPlot
psweep$ldPlot


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


library(ape)

# ---- Example comparison ----
set.seed(1)

tree_100 <- simulate_coalescent_tree(n = 8, Ne = 100)
tree_500 <- simulate_coalescent_tree(n = 8, Ne = 500)

y_max <- max(
  max(node.depth.edgelength(tree_100)),
  max(node.depth.edgelength(tree_500))
)

par(mfrow = c(1,2))
plot_coalescent_tree_fixed(tree_100, main = "Ne = 100", y_max = y_max)
plot_coalescent_tree_fixed(tree_500, main = "Ne = 500", y_max = y_max)
par(mfrow = c(1,1))

set.seed(1)
tree_base <- rcoal(8)

tree_100 <- tree_base
tree_500 <- tree_base

tree_100$edge.length <- tree_base$edge.length * (2 * 100)
tree_500$edge.length <- tree_base$edge.length * (2 * 500)

y_max <- max(max(node.depth.edgelength(tree_100)),
             max(node.depth.edgelength(tree_500)))

par(mfrow = c(1,2))
plot_coalescent_tree_fixed(tree_100, main = "Same genealogy, Ne = 100", y_max = y_max)
plot_coalescent_tree_fixed(tree_500, main = "Same genealogy, Ne = 500", y_max = y_max)
par(mfrow = c(1,1))


plot_fixation_probability(N = 50, p0 = 0.3)
plot_Fst_vs_migration(Ne = 100)
plot_mutation_selection_balance(mu = 1e-5, s = 0.1)
plot_site_frequency_spectrum(n = 20)


plot_population_genetics_concept_map()
plot_population_genetics_concept_map_minimal()



# install.packages("ggplot2")  # if needed
library(ggplot2)

set.seed(1)

# ---- parameters you can change ----
n_each <- 10
mix_red  <- 7
mix_blue <- 3
bucket_w <- 1.0
bucket_h <- 1.2

# helper to place balls inside a bucket area
make_bucket_points <- function(n, bucket, color, x0, x1, y0, y1) {
  data.frame(
    bucket = bucket,
    color  = color,
    x = runif(n, x0 + 0.08, x1 - 0.08),
    y = runif(n, y0 + 0.10, y1 - 0.10)
  )
}

# bucket coordinates in a common plotting space
buckets <- data.frame(
  bucket = c("Population A", "Population B", "Admixed"),
  xmin   = c(0, 1.3, 2.6),
  xmax   = c(1.0, 2.3, 3.6),
  ymin   = 0,
  ymax   = bucket_h
)

# points (balls)
pts <- rbind(
  make_bucket_points(n_each, "Population A", "red",
                     buckets$xmin[1], buckets$xmax[1], buckets$ymin[1], buckets$ymax[1]),
  make_bucket_points(n_each, "Population B", "blue",
                     buckets$xmin[2], buckets$xmax[2], buckets$ymin[2], buckets$ymax[2]),
  make_bucket_points(mix_red, "Admixed", "red",
                     buckets$xmin[3], buckets$xmax[3], buckets$ymin[3], buckets$ymax[3]),
  make_bucket_points(mix_blue, "Admixed", "blue",
                     buckets$xmin[3], buckets$xmax[3], buckets$ymin[3], buckets$ymax[3])
)

# color mapping
cols <- c(red = "#D94848", blue = "#2B6CB0")

ggplot() +
  # draw bucket outlines
  geom_rect(
    data = buckets,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "grey30", linewidth = 0.7
  ) +
  # bucket labels
  geom_text(
    data = buckets,
    aes(x = (xmin + xmax)/2, y = ymax + 0.12, label = bucket),
    fontface = "bold", size = 4
  ) +
  # balls
  geom_point(
    data = pts,
    aes(x = x, y = y, color = color),
    size = 5
  ) +
  scale_color_manual(values = cols, guide = "none") +
  coord_equal(xlim = c(-0.1, 3.7), ylim = c(-0.1, 1.45), expand = FALSE) +
  theme_void() +
  ggtitle("Illustration of admixture: two source populations and an admixed population") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))




plot_admixture_teaching <- function(
    group_sizes,
    alpha_list,
    group_names = NULL,
    colors = NULL,
    seed = 123,
    title = NULL
) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(MCMCpack)
  
  set.seed(seed)
  
  K <- length(alpha_list[[1]])
  
  if (is.null(group_names)) {
    group_names <- paste0("Group ", seq_along(group_sizes))
  }
  
  if (is.null(colors)) {
    colors <- scales::hue_pal()(K)
  }
  
  # Generate ancestry proportions
  props_list <- lapply(seq_along(group_sizes), function(i) {
    rdirichlet(group_sizes[i], alpha_list[[i]])
  })
  
  props <- do.call(rbind, props_list)
  
  df <- data.frame(
    individual = 1:nrow(props),
    group = rep(group_names, group_sizes),
    props
  )
  
  colnames(df)[3:(K+2)] <- paste0("Cluster", 1:K)
  
  df <- df %>% arrange(group)
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("Cluster"),
                 names_to = "ancestry",
                 values_to = "proportion")
  
  ggplot(df_long, aes(x = factor(individual),
                      y = proportion,
                      fill = ancestry)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Individuals",
         y = "Ancestry proportion",
         title = title) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_admixture_teaching(
  group_sizes = c(10,10,10,10),
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
  colors = c("#D94848","#2B6CB0","#2F9E44"),
  title = "ADMIXTURE-style barplot (K = 3)"
)



plot_migration_pulse <- function(
    n_per_pop = 40,
    generations = 6,
    pulse_generation = 2,
    migration_rate = 0.2,
    seed = 123
) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  set.seed(seed)
  
  # Initial ancestry
  popA <- rep(1, n_per_pop)
  popB <- rep(0, n_per_pop)
  
  results <- list()
  
  for (g in 0:(generations - 1)) {
    
    if (g == pulse_generation) {
      new_popA <- (1 - migration_rate) * popA + migration_rate * popB
      new_popB <- (1 - migration_rate) * popB + migration_rate * popA
      popA <- new_popA
      popB <- new_popB
    }
    
    df <- data.frame(
      generation = g,
      population = rep(c("Pop A", "Pop B"), each = n_per_pop),
      ancestry_A = c(popA, popB)
    )
    
    df$ancestry_B <- 1 - df$ancestry_A
    
    results[[g + 1]] <- df
  }
  
  df_all <- bind_rows(results)
  
  # Create individual index within each generation
  df_all <- df_all %>%
    group_by(generation) %>%
    mutate(individual = row_number()) %>%
    ungroup()
  
  df_long <- df_all %>%
    pivot_longer(cols = c(ancestry_A, ancestry_B),
                 names_to = "cluster",
                 values_to = "proportion")
  
  ggplot(df_long,
         aes(x = individual,
             y = proportion,
             fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = c("#D94848","#2B6CB0")) +
    facet_wrap(~ generation, ncol = 1) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = paste("Migration Pulse at Generation", pulse_generation),
      y = "Ancestry proportion",
      x = NULL
    )
}


install.packages("gganimate")
install.packages("gifski")

plot_migration_pulse_3pop <- function(
    n_per_pop = 30,
    generations = 8,
    pulse_generation = 2,
    migration_rate = 0.2,
    seed = 123,
    animate = TRUE
) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gganimate)
  
  set.seed(seed)
  
  # Initial ancestry vectors
  popA <- matrix(rep(c(1,0,0), each = n_per_pop), ncol = 3, byrow = TRUE)
  popB <- matrix(rep(c(0,1,0), each = n_per_pop), ncol = 3, byrow = TRUE)
  popC <- matrix(rep(c(0,0,1), each = n_per_pop), ncol = 3, byrow = TRUE)
  
  results <- list()
  
  for (g in 0:(generations - 1)) {
    
    if (g == pulse_generation) {
      
      # A <-> B exchange
      newA <- (1 - migration_rate) * popA +
        migration_rate * popB
      newB <- (1 - migration_rate) * popB +
        migration_rate * popA
      
      popA <- newA
      popB <- newB
      # popC unchanged
    }
    
    df <- data.frame(
      generation = g,
      population = rep(c("Pop A", "Pop B", "Pop C"),
                       each = n_per_pop),
      ancestry_A = c(popA[,1], popB[,1], popC[,1]),
      ancestry_B = c(popA[,2], popB[,2], popC[,2]),
      ancestry_C = c(popA[,3], popB[,3], popC[,3])
    )
    
    results[[g + 1]] <- df
  }
  
  df_all <- bind_rows(results)
  
  df_all <- df_all %>%
    group_by(generation) %>%
    mutate(individual = row_number()) %>%
    ungroup()
  
  df_long <- df_all %>%
    pivot_longer(cols = starts_with("ancestry"),
                 names_to = "cluster",
                 values_to = "proportion")
  
  p <- ggplot(df_long,
              aes(x = individual,
                  y = proportion,
                  fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = c(
      ancestry_A = "#D94848",
      ancestry_B = "#2B6CB0",
      ancestry_C = "#2F9E44"
    )) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(
      y = "Ancestry proportion",
      x = NULL,
      title = "Migration Pulse (3 Populations)\nGeneration: {closest_state}"
    )
  
  if (animate) {
    p <- p +
      transition_states(generation,
                        transition_length = 1,
                        state_length = 1) +
      ease_aes("linear")
  } else {
    p <- p + facet_wrap(~ generation, ncol = 1)
  }
  
  return(p)
}

plot_recombination_two_loci <- function() {
  
  library(ggplot2)
  
  # Data for chromosomes before recombination
  df_before <- data.frame(
    x = c(0, 10),
    y = c(3, 2),
    haplotype = c("AB", "ab")
  )
  
  # Data for chromosomes after recombination
  df_after <- data.frame(
    x = c(0, 10),
    y = c(0, -1),
    haplotype = c("Ab", "aB")
  )
  
  # Base plot
  p <- ggplot() +
    
    # Before recombination chromosomes
    geom_segment(aes(x = 0, xend = 10, y = 3, yend = 3), linewidth = 2) +
    geom_segment(aes(x = 0, xend = 10, y = 2, yend = 2), linewidth = 2) +
    
    # After recombination chromosomes
    geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0), linewidth = 2) +
    geom_segment(aes(x = 0, xend = 10, y = -1, yend = -1), linewidth = 2) +
    
    # Locus positions
    geom_point(aes(x = 3, y = 3), size = 4) +
    geom_point(aes(x = 7, y = 3), size = 4) +
    geom_point(aes(x = 3, y = 2), size = 4) +
    geom_point(aes(x = 7, y = 2), size = 4) +
    
    geom_point(aes(x = 3, y = 0), size = 4) +
    geom_point(aes(x = 7, y = 0), size = 4) +
    geom_point(aes(x = 3, y = -1), size = 4) +
    geom_point(aes(x = 7, y = -1), size = 4) +
    
    # Labels for alleles (before)
    annotate("text", x = 3, y = 3.3, label = "A", size = 6) +
    annotate("text", x = 7, y = 3.3, label = "B", size = 6) +
    annotate("text", x = 3, y = 2.3, label = "a", size = 6) +
    annotate("text", x = 7, y = 2.3, label = "b", size = 6) +
    
    # Labels for alleles (after)
    annotate("text", x = 3, y = 0.3, label = "A", size = 6) +
    annotate("text", x = 7, y = 0.3, label = "b", size = 6) +
    annotate("text", x = 3, y = -0.7, label = "a", size = 6) +
    annotate("text", x = 7, y = -0.7, label = "B", size = 6) +
    
    # Crossover indicator
    geom_segment(aes(x = 5, xend = 5, y = 1.6, yend = 3.4),
                 linetype = "dashed", linewidth = 1) +
    
    annotate("text", x = 5, y = 3.7, label = "Crossover", size = 5) +
    
    annotate("text", x = 5, y = 3.9,
             label = "Parental haplotypes (before recombination)",
             size = 5) +
    
    annotate("text", x = 5, y = 1,
             label = "Recombinant haplotypes (after recombination)",
             size = 5) +
    
    coord_cartesian(xlim = c(-1, 11), ylim = c(-2, 5)) +
    theme_void()
  
  return(p)
}

plot_recombination_two_loci()

plot_recombination_colored <- function() {
  
  library(ggplot2)
  
  # Crossover position
  x_cross <- 5
  
  ggplot() +
    
    # ---- BEFORE RECOMBINATION ----
  
  # Chromosome 1 (AB) - red
  geom_segment(aes(x = 0, xend = 10, y = 3, yend = 3),
               linewidth = 3, color = "#D94848") +
    
    # Chromosome 2 (ab) - blue
    geom_segment(aes(x = 0, xend = 10, y = 2, yend = 2),
                 linewidth = 3, color = "#2B6CB0") +
    
    # Loci markers (before)
    geom_point(aes(x = 3, y = 3), size = 4) +
    geom_point(aes(x = 7, y = 3), size = 4) +
    geom_point(aes(x = 3, y = 2), size = 4) +
    geom_point(aes(x = 7, y = 2), size = 4) +
    
    annotate("text", x = 3, y = 3.3, label = "A", size = 6) +
    annotate("text", x = 7, y = 3.3, label = "B", size = 6) +
    annotate("text", x = 3, y = 2.3, label = "a", size = 6) +
    annotate("text", x = 7, y = 2.3, label = "b", size = 6) +
    
    # Crossover line
    geom_segment(aes(x = x_cross, xend = x_cross,
                     y = 1.6, yend = 3.4),
                 linetype = "dashed", linewidth = 1) +
    
    annotate("text", x = x_cross, y = 3.6,
             label = "Crossover", size = 5) +
    
    # ---- AFTER RECOMBINATION ----
  
  # Recombinant 1 (Ab)
  geom_segment(aes(x = 0, xend = x_cross, y = 0, yend = 0),
               linewidth = 3, color = "#D94848") +
    geom_segment(aes(x = x_cross, xend = 10, y = 0, yend = 0),
                 linewidth = 3, color = "#2B6CB0") +
    
    # Recombinant 2 (aB)
    geom_segment(aes(x = 0, xend = x_cross, y = -1, yend = -1),
                 linewidth = 3, color = "#2B6CB0") +
    geom_segment(aes(x = x_cross, xend = 10, y = -1, yend = -1),
                 linewidth = 3, color = "#D94848") +
    
    # Loci markers (after)
    geom_point(aes(x = 3, y = 0), size = 4) +
    geom_point(aes(x = 7, y = 0), size = 4) +
    geom_point(aes(x = 3, y = -1), size = 4) +
    geom_point(aes(x = 7, y = -1), size = 4) +
    
    annotate("text", x = 3, y = 0.3, label = "A", size = 6) +
    annotate("text", x = 7, y = 0.3, label = "b", size = 6) +
    annotate("text", x = 3, y = -0.7, label = "a", size = 6) +
    annotate("text", x = 7, y = -0.7, label = "B", size = 6) +
    
    annotate("text", x = 5, y = 4.2,
             label = "Parental haplotypes (before recombination)",
             size = 5) +
    
    annotate("text", x = 5, y = 1,
             label = "Recombinant haplotypes (after recombination)",
             size = 5) +
    
    coord_cartesian(xlim = c(-1, 11), ylim = c(-2, 5)) +
    theme_void()
}


plot_recombination_colored()


plot_haplotype_inheritance <- function() {
  
  library(ggplot2)
  
  xA <- 3
  xB <- 7
  x_cross <- 5
  
  ggplot() +
    
    # =========================
  # MOTHER
  # =========================
  annotate("text", x = 5, y = 5.2, label = "Mother", size = 6, fontface = "bold") +
    
    # Maternal haplotypes (red)
    geom_segment(aes(x = 0, xend = 10, y = 4.5, yend = 4.5),
                 linewidth = 3, color = "#D94848") +
    geom_segment(aes(x = 0, xend = 10, y = 4, yend = 4),
                 linewidth = 3, color = "#D94848") +
    
    annotate("text", x = xA, y = 4.8, label = "A", size = 6) +
    annotate("text", x = xB, y = 4.8, label = "B", size = 6) +
    annotate("text", x = xA, y = 4.3, label = "a", size = 6) +
    annotate("text", x = xB, y = 4.3, label = "b", size = 6) +
    
    # =========================
  # FATHER
  # =========================
  annotate("text", x = 5, y = 3.3, label = "Father", size = 6, fontface = "bold") +
    
    geom_segment(aes(x = 0, xend = 10, y = 2.8, yend = 2.8),
                 linewidth = 3, color = "#2B6CB0") +
    geom_segment(aes(x = 0, xend = 10, y = 2.3, yend = 2.3),
                 linewidth = 3, color = "#2B6CB0") +
    
    annotate("text", x = xA, y = 3.1, label = "A", size = 6) +
    annotate("text", x = xB, y = 3.1, label = "b", size = 6) +
    annotate("text", x = xA, y = 2.6, label = "a", size = 6) +
    annotate("text", x = xB, y = 2.6, label = "B", size = 6) +
    
    # =========================
  # CASE 1: NO CROSSOVER
  # =========================
  annotate("text", x = 5, y = 1.8,
           label = "Case 1: No Crossover",
           size = 5, fontface = "bold") +
    
    # Offspring haplotypes (intact)
    geom_segment(aes(x = 0, xend = 10, y = 1.3, yend = 1.3),
                 linewidth = 3, color = "#D94848") +
    geom_segment(aes(x = 0, xend = 10, y = 0.8, yend = 0.8),
                 linewidth = 3, color = "#2B6CB0") +
    
    annotate("text", x = xA, y = 1.6, label = "A", size = 6) +
    annotate("text", x = xB, y = 1.6, label = "B", size = 6) +
    annotate("text", x = xA, y = 1.1, label = "a", size = 6) +
    annotate("text", x = xB, y = 1.1, label = "B", size = 6) +
    
    annotate("text", x = 5, y = 0.3,
             label = "Genotype: A a  B B",
             size = 5) +
    
    # =========================
  # CASE 2: CROSSOVER
  # =========================
  annotate("text", x = 5, y = -0.5,
           label = "Case 2: Crossover",
           size = 5, fontface = "bold") +
    
    # Crossover indicator
    geom_segment(aes(x = x_cross, xend = x_cross,
                     y = 2.1, yend = 4.9),
                 linetype = "dashed", linewidth = 1) +
    
    # Recombinant gametes
    geom_segment(aes(x = 0, xend = x_cross, y = -1, yend = -1),
                 linewidth = 3, color = "#D94848") +
    geom_segment(aes(x = x_cross, xend = 10, y = -1, yend = -1),
                 linewidth = 3, color = "#2B6CB0") +
    
    geom_segment(aes(x = 0, xend = x_cross, y = -1.5, yend = -1.5),
                 linewidth = 3, color = "#2B6CB0") +
    geom_segment(aes(x = x_cross, xend = 10, y = -1.5, yend = -1.5),
                 linewidth = 3, color = "#D94848") +
    
    annotate("text", x = xA, y = -0.7, label = "A", size = 6) +
    annotate("text", x = xB, y = -0.7, label = "b", size = 6) +
    annotate("text", x = xA, y = -1.2, label = "a", size = 6) +
    annotate("text", x = xB, y = -1.2, label = "B", size = 6) +
    
    coord_cartesian(xlim = c(-1, 11), ylim = c(-2.5, 6)) +
    theme_void()
}

plot_haplotype_inheritance()


draw_haplotype <- function(y,
                           alleles,
                           color_left,
                           color_right = NULL,
                           xA = 3,
                           xB = 7,
                           x_cross = 5,
                           allele_size = 5,
                           locus_height = 0.15) {
  
  layers <- list()
  
  # ---- Chromosome segments ----
  if (is.null(color_right)) {
    
    layers <- c(layers,
                list(
                  geom_segment(aes(x = 0, xend = 10, y = y, yend = y),
                               linewidth = 3,
                               color = color_left)
                )
    )
    
  } else {
    
    layers <- c(layers,
                list(
                  geom_segment(aes(x = 0, xend = x_cross, y = y, yend = y),
                               linewidth = 3,
                               color = color_left),
                  geom_segment(aes(x = x_cross, xend = 10, y = y, yend = y),
                               linewidth = 3,
                               color = color_right)
                )
    )
  }
  
  # ---- Locus tick marks ----
  layers <- c(layers,
              list(
                geom_segment(aes(x = xA, xend = xA,
                                 y = y - locus_height,
                                 yend = y + locus_height),
                             linewidth = 1),
                geom_segment(aes(x = xB, xend = xB,
                                 y = y - locus_height,
                                 yend = y + locus_height),
                             linewidth = 1)
              )
  )
  
  # ---- Allele labels ----
  layers <- c(layers,
              list(
                annotate("text",
                         x = xA,
                         y = y + 0.3,
                         label = alleles[1],
                         size = allele_size),
                annotate("text",
                         x = xB,
                         y = y + 0.3,
                         label = alleles[2],
                         size = allele_size)
              )
  )
  
  return(layers)
}

p1 <- ggplot() +
  draw_haplotype(1,
                 c("A","B"),
                 color_left = "#D94848") +
  coord_cartesian(xlim = c(0,10), ylim = c(0,2)) +
  theme_void()

p2 <- ggplot() +
  draw_haplotype(1,
                 c("A","b"),
                 color_left = "#D94848",
                 color_right = "#2B6CB0") +
  coord_cartesian(xlim = c(0,10), ylim = c(0,2)) +
  theme_void()

(p1/p2/p1/p1) | (p1/p2/p1/p1)

((p1 / p1 / plot_spacer() / p2 / p1) |
    (p2 / p2 / plot_spacer() / p2 / p1))

draw_haplotype(0.7, c("A","b"), "#D94848", "#2B6CB0")

library(ggplot2)
library(patchwork)

xA <- 3
xB <- 7
x_cross <- 5

p_case1 <- ggplot() +
  
  annotate("text", x = 5, y = 5, label = "Case 1: No Crossover",
           fontface = "bold", size = 6) +
  
  # Mother
  annotate("text", x = 5, y = 4.4, label = "Mother", size = 5) +
  draw_haplotype(4, c("A","B"), "#D94848") +
  draw_haplotype(3.5, c("a","b"), "#D94848") +
  
  # Father
  annotate("text", x = 5, y = 2.8, label = "Father", size = 5) +
  draw_haplotype(2.4, c("A","b"), "#2B6CB0") +
  draw_haplotype(1.9, c("a","B"), "#2B6CB0") +
  
  # Offspring example
  annotate("text", x = 5, y = 1.1, label = "Offspring",
           fontface = "bold", size = 5) +
  draw_haplotype(0.7, c("A","B"), "#D94848") +
  draw_haplotype(0.2, c("a","B"), "#2B6CB0") +
  
  annotate("text", x = xA, y = 4.2, label = "A", size = 5) +
  annotate("text", x = xB, y = 4.2, label = "B", size = 5) +
  annotate("text", x = xA, y = 3.7, label = "a", size = 5) +
  annotate("text", x = xB, y = 3.7, label = "b", size = 5) +
  
  annotate("text", x = xA, y = 2.6, label = "A", size = 5) +
  annotate("text", x = xB, y = 2.6, label = "b", size = 5) +
  annotate("text", x = xA, y = 2.1, label = "a", size = 5) +
  annotate("text", x = xB, y = 2.1, label = "B", size = 5) +
  
  annotate("text", x = xA, y = 0.9, label = "A", size = 5) +
  annotate("text", x = xB, y = 0.9, label = "B", size = 5) +
  annotate("text", x = xA, y = 0.4, label = "a", size = 5) +
  annotate("text", x = xB, y = 0.4, label = "B", size = 5) +
  
  coord_cartesian(xlim = c(-1, 11), ylim = c(-0.5, 5.5)) +
  theme_void()

p_case2 <- ggplot() +
  
  annotate("text", x = 5, y = 5, label = "Case 2: Crossover",
           fontface = "bold", size = 6) +
  
  # Mother
  annotate("text", x = 5, y = 4.4, label = "Mother", size = 5) +
  draw_haplotype(4, c("A","B"), "#D94848") +
  draw_haplotype(3.5, c("a","b"), "#D94848") +
  
  # Father
  annotate("text", x = 5, y = 2.8, label = "Father", size = 5) +
  draw_haplotype(2.4, c("A","b"), "#2B6CB0") +
  draw_haplotype(1.9, c("a","B"), "#2B6CB0") +
  
  # Crossover line
  geom_segment(aes(x = x_cross, xend = x_cross,
                   y = 3.2, yend = 4.8),
               linetype = "dashed", linewidth = 1) +
  
  # Recombinant gametes
  annotate("text", x = 5, y = 1.1, label = "Recombinant Gametes",
           fontface = "bold", size = 5) +
  draw_haplotype(0.7, c("A","b"),
                 "#D94848", "#2B6CB0") +
  draw_haplotype(0.2, c("a","B"),
                 "#2B6CB0", "#D94848") +
  
  annotate("text", x = xA, y = 0.9, label = "A", size = 5) +
  annotate("text", x = xB, y = 0.9, label = "b", size = 5) +
  annotate("text", x = xA, y = 0.4, label = "a", size = 5) +
  annotate("text", x = xB, y = 0.4, label = "B", size = 5) +
  
  coord_cartesian(xlim = c(-1, 11), ylim = c(-0.5, 5.5)) +
  theme_void()



# # ============================================================
# # BIOLOGICAL TEACHING PLOT FRAMEWORK
# # ============================================================
# 
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # ============================================================
# # 1) BIOLOGICAL COLOUR SYSTEM
# # ============================================================
# 
# bio_cols <- list(
#   
#   allele = c(
#     A = "#D55E00",   # Vermillion (Okabe-Ito)
#     a = "#0072B2"    # Blue
#   ),
#   
#   genotype = c(
#     AA = "#D55E00",
#     Aa = "#8C6BB1",
#     aa = "#0072B2"
#   ),
#   
#   stage = c(
#     "Before" = "grey70",
#     "After"  = "#D55E00"
#   ),
#   
#   model = c(
#     "Observed" = "grey60",
#     "Expected (HWE)" = "grey85",
#     "HWE" = "grey70",
#     "Inbreeding" = "#D55E00"
#   ),
#   
#   ancestry = c(
#     "Ancestry 1" = "#D55E00",
#     "Ancestry 2" = "#0072B2"
#   )
# )
# 
# # ============================================================
# # 2) UNIFIED BIO THEME
# # ============================================================
# 
# theme_bio <- function(base_size = 18,
#                       legend_position = "top") {
#   
#   theme_minimal(base_size = base_size) +
#     theme(
#       plot.title = element_text(face = "plain"),
#       plot.subtitle = element_text(color = "grey30"),
#       panel.grid.minor = element_blank(),
#       legend.position = legend_position,
#       legend.title = element_blank()
#     )
# }
# 
# # ============================================================
# # 3) GENERIC BIO BAR PLOT
# # ============================================================
# 
# plot_bio_bar <- function(data,
#                          x,
#                          y,
#                          fill,
#                          palette,
#                          ylim = NULL,
#                          title = NULL,
#                          subtitle = NULL,
#                          ylab = "Frequency",
#                          base_size = 18) {
#   
#   g <- ggplot(data, aes({{x}}, {{y}}, fill = {{fill}})) +
#     geom_col(position = position_dodge(width = 0.7),
#              width = 0.6,
#              color = "black") +
#     scale_fill_manual(values = palette) +
#     labs(
#       title = title,
#       subtitle = subtitle,
#       x = NULL,
#       y = ylab
#     ) +
#     theme_bio(base_size)
#   
#   if (!is.null(ylim)) {
#     g <- g + coord_cartesian(ylim = ylim)
#   }
#   
#   g
# }
# 
# # ============================================================
# # 4) WRIGHT–FISHER DRIFT SIMULATION
# # ============================================================
# 
# simulate_drift <- function(N,
#                            generations = 100,
#                            replicates = 20,
#                            p0 = 0.7) {
#   
#   results <- vector("list", replicates)
#   
#   for (r in seq_len(replicates)) {
#     
#     p <- numeric(generations + 1)
#     p[1] <- p0
#     
#     for (t in seq_len(generations)) {
#       p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
#     }
#     
#     results[[r]] <- data.frame(
#       generation = 0:generations,
#       frequency  = p,
#       replicate  = r
#     )
#   }
#   
#   bind_rows(results)
# }
# 
# plot_drift <- function(N,
#                        generations = 100,
#                        replicates = 20,
#                        p0 = 0.7,
#                        show_mean = TRUE,
#                        base_size = 18) {
#   
#   df <- simulate_drift(N, generations, replicates, p0)
#   
#   g <- ggplot(df,
#               aes(generation, frequency,
#                   group = replicate)) +
#     geom_line(alpha = 0.4) +
#     geom_hline(yintercept = p0,
#                linetype = "dashed",
#                color = "grey40") +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = paste("Genetic Drift (N =", N, ")"),
#       x = "Generation",
#       y = "Allele frequency (p)"
#     ) +
#     theme_bio(base_size)
#   
#   if (show_mean) {
#     g <- g +
#       stat_summary(aes(group = 1),
#                    fun = mean,
#                    geom = "line",
#                    linewidth = 1.5,
#                    color = bio_cols$allele["A"])
#   }
#   
#   g
# }
# 
# # ============================================================
# # 5) SELECTION MODEL
# # ============================================================
# 
# step_selection <- function(p, wAA, wAa, waa) {
#   q <- 1 - p
#   wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
#   (p^2*wAA + p*q*wAa) / wbar
# }
# 
# simulate_selection <- function(p0,
#                                T,
#                                wAA,
#                                wAa,
#                                waa) {
#   
#   p <- numeric(T + 1)
#   p[1] <- p0
#   
#   for (t in seq_len(T)) {
#     p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
#   }
#   
#   data.frame(
#     generation = 0:T,
#     p = p
#   )
# }
# 
# plot_selection_trajectory <- function(p0,
#                                       T,
#                                       wAA,
#                                       wAa,
#                                       waa,
#                                       title = "Selection",
#                                       base_size = 18) {
#   
#   df <- simulate_selection(p0, T, wAA, wAa, waa)
#   
#   ggplot(df, aes(generation, p)) +
#     geom_line(linewidth = 1.5,
#               color = bio_cols$allele["A"]) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = title,
#       x = "Generation",
#       y = "Allele frequency (p)"
#     ) +
#     theme_bio(base_size)
# }
# 
# # ============================================================
# # 6) HARDY–WEINBERG COMPARISON
# # ============================================================
# 
# plot_hwe_comparison <- function(geno_counts,
#                                 base_size = 18) {
#   
#   # ----------------------------
#   # Input checks
#   # ----------------------------
#   
#   if (is.null(names(geno_counts))) {
#     stop("geno_counts must be named, e.g. c(AA=10, Aa=7, aa=10)")
#   }
#   
#   required <- c("AA", "Aa", "aa")
#   geno_counts[setdiff(required, names(geno_counts))] <- 0
#   geno_counts <- geno_counts[required]
#   
#   N <- sum(geno_counts)
#   
#   # ----------------------------
#   # Allele frequency
#   # ----------------------------
#   
#   nA <- 2 * geno_counts["AA"] + geno_counts["Aa"]
#   p_hat <- as.numeric(nA) / (2 * N)
#   q_hat <- 1 - p_hat
#   
#   # ----------------------------
#   # Expected counts
#   # ----------------------------
#   
#   expected <- c(
#     AA = p_hat^2,
#     Aa = 2*p_hat*q_hat,
#     aa = q_hat^2
#   ) * N
#   
#   # ----------------------------
#   # Chi-square test
#   # ----------------------------
#   
#   chisq <- sum((geno_counts - expected)^2 / expected)
#   pval  <- pchisq(chisq, df = 1, lower.tail = FALSE)
#   
#   # ----------------------------
#   # Prepare plotting data
#   # ----------------------------
#   
#   df <- data.frame(
#     genotype = rep(required, 2),
#     type = rep(c("Observed", "Expected (HWE)"), each = 3),
#     n = c(geno_counts, expected)
#   )
#   
#   # ----------------------------
#   # Color scheme
#   # ----------------------------
#   
#   hwe_cols <- c(
#     "Observed" = "#2C3E50",
#     "Expected (HWE)" = "#D5DBDB"
#   )
#   
#   # ----------------------------
#   # Plot
#   # ----------------------------
#   
#   ggplot2::ggplot(df,
#                   ggplot2::aes(genotype, n, fill = type)) +
#     
#     ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
#                       width = 0.6,
#                       color = "black") +
#     
#     ggplot2::scale_fill_manual(values = hwe_cols) +
#     
#     ggplot2::labs(
#       title = "Hardy–Weinberg Equilibrium Test",
#       subtitle = paste0(
#         "N = ", N,
#         "  |  Estimated p = ", round(p_hat, 3),
#         "  |  χ² = ", round(chisq, 3),
#         "  |  P-value = ", signif(pval, 3)
#       ),
#       x = NULL,
#       y = "Genotype count",
#       fill = ""
#     ) +
#     
#     ggplot2::coord_cartesian(ylim = c(0, max(df$n) * 1.15)) +
#     
#     ggplot2::theme_minimal(base_size = base_size) +
#     ggplot2::theme(
#       panel.grid.minor = ggplot2::element_blank(),
#       legend.position = "top",
#       plot.title = ggplot2::element_text(face = "bold")
#     )
# }
# 
# # ============================================================
# # END OF FRAMEWORK
# # ============================================================
# 
# 
# 
# 
# # ============================================================
# # BIOLOGICAL TEACHING PLOT FRAMEWORK (EXTENDED)
# # ============================================================
# 
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(patchwork)
# 
# # ============================================================
# # 1) BIOLOGICAL COLOUR SYSTEM
# # ============================================================
# 
# bio_cols <- list(
#   
#   allele = c(
#     A = "#D55E00",
#     a = "#0072B2"
#   ),
#   
#   genotype = c(
#     AA = "#D55E00",
#     Aa = "#8C6BB1",
#     aa = "#0072B2"
#   ),
#   
#   stage = c(
#     "Gen 0"  = "grey70",
#     "After"  = "#D55E00",
#     "Before" = "grey70"
#   ),
#   
#   model = c(
#     "Observed" = "grey60",
#     "Expected (HWE)" = "grey85",
#     "HWE" = "grey70",
#     "Inbreeding" = "#D55E00"
#   ),
#   
#   ancestry = c(
#     "Ancestry 1" = "#D55E00",
#     "Ancestry 2" = "#0072B2"
#   )
# )
# 
# # ============================================================
# # 2) UNIFIED BIO THEME
# # ============================================================
# 
# theme_bio <- function(base_size = 18,
#                       legend_position = "top") {
#   
#   theme_minimal(base_size = base_size) +
#     theme(
#       plot.title = element_text(face = "bold"),
#       plot.subtitle = element_text(color = "grey30"),
#       panel.grid.minor = element_blank(),
#       legend.position = legend_position,
#       legend.title = element_blank()
#     )
# }
# 
# # ============================================================
# # 3) GENOTYPE / ALLELE FREQUENCY HELPERS
# # ============================================================
# 
# geno_freq <- function(p) {
#   q <- 1 - p
#   tibble(
#     genotype = c("AA","Aa","aa"),
#     freq = c(p^2, 2*p*q, q^2)
#   )
# }
# 
# allele_freq <- function(p) {
#   tibble(
#     allele = c("A","a"),
#     freq = c(p, 1-p)
#   )
# }
# 
# plot_genotype_freq <- function(p = NULL,
#                                geno_counts = NULL,
#                                stage_label = NULL,
#                                base_size = 18) {
#   
#   if (!is.null(geno_counts)) {
#     
#     if (is.null(names(geno_counts)))
#       stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")
#     
#     required <- c("AA","Aa","aa")
#     geno_counts[setdiff(required, names(geno_counts))] <- 0
#     geno_counts <- geno_counts[required]
#     
#     N <- sum(geno_counts)
#     
#     df <- data.frame(
#       genotype = required,
#       freq = as.numeric(geno_counts) / N
#     )
#     
#   } else if (!is.null(p)) {
#     
#     q <- 1 - p
#     
#     df <- data.frame(
#       genotype = c("AA","Aa","aa"),
#       freq = c(p^2, 2*p*q, q^2)
#     )
#     
#   } else {
#     stop("Provide either p OR geno_counts")
#   }
#   
#   # Title logic
#   title_text <- if (is.null(stage_label)) {
#     "Genotype frequencies"
#   } else {
#     paste0("Genotype frequencies (", stage_label, ")")
#   }
#   
#   ggplot(df, aes(genotype, freq, fill = genotype)) +
#     geom_col(width = 0.6, color = "black") +
#     scale_fill_manual(values = bio_cols$genotype) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = title_text,
#       y = "Frequency",
#       x = NULL
#     ) +
#     theme_bio(base_size)
# }
# 
# plot_allele_freq <- function(p = NULL,
#                              geno_counts = NULL,
#                              stage_label = NULL,
#                              base_size = 18) {
#   
#   if (!is.null(geno_counts)) {
#     
#     if (is.null(names(geno_counts)))
#       stop("geno_counts must be named: c(AA=..., Aa=..., aa=...)")
#     
#     required <- c("AA","Aa","aa")
#     geno_counts[setdiff(required, names(geno_counts))] <- 0
#     geno_counts <- geno_counts[required]
#     
#     N <- sum(geno_counts)
#     
#     nA <- 2*geno_counts["AA"] + geno_counts["Aa"]
#     p_hat <- as.numeric(nA) / (2*N)
#     q_hat <- 1 - p_hat
#     
#     df <- data.frame(
#       allele = c("A","a"),
#       freq = c(p_hat, q_hat)
#     )
#     
#   } else if (!is.null(p)) {
#     
#     df <- data.frame(
#       allele = c("A","a"),
#       freq = c(p, 1-p)
#     )
#     
#   } else {
#     stop("Provide either p OR geno_counts")
#   }
#   
#   title_text <- if (is.null(stage_label)) {
#     "Allele frequencies"
#   } else {
#     paste0("Allele frequencies (", stage_label, ")")
#   }
#   
#   ggplot(df, aes(allele, freq, fill = allele)) +
#     geom_col(width = 0.6, color = "black") +
#     scale_fill_manual(values = bio_cols$allele) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = title_text,
#       y = "Frequency",
#       x = NULL
#     ) +
#     theme_bio(base_size)
# }
# 
# # ============================================================
# # 4) INBREEDING WRAPPER
# # ============================================================
# 
# plot_inbreeding <- function(p, F, base_size = 18) {
#   
#   q <- 1 - p
#   
#   geno_hwe <- tibble(
#     genotype = c("AA","Aa","aa"),
#     freq = c(p^2, 2*p*q, q^2),
#     model = "HWE"
#   )
#   
#   geno_F <- tibble(
#     genotype = c("AA","Aa","aa"),
#     freq = c(
#       p^2 + F*p*q,
#       2*p*q*(1 - F),
#       q^2 + F*p*q
#     ),
#     model = "Inbreeding"
#   )
#   
#   geno_long <- bind_rows(geno_hwe, geno_F)
#   
#   allele_df <- tibble(
#     allele = rep(c("A","a"), 2),
#     freq = rep(c(p, q), 2),
#     model = rep(c("HWE","Inbreeding"), each = 2)
#   )
#   
#   p1 <- ggplot(geno_long,
#                aes(genotype, freq, fill = model)) +
#     geom_col(position = position_dodge(width = 0.7),
#              width = 0.6,
#              color = "black") +
#     scale_fill_manual(values = bio_cols$model) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(title = "Genotype frequencies",
#          y = "Frequency", x = NULL) +
#     theme_bio(base_size)
#   
#   p2 <- ggplot(allele_df,
#                aes(allele, freq, fill = model)) +
#     geom_col(position = position_dodge(width = 0.7),
#              width = 0.6,
#              color = "black") +
#     scale_fill_manual(values = bio_cols$model) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(title = "Allele frequencies (unchanged)",
#          y = "Frequency", x = NULL) +
#     theme_bio(base_size)
#   
#   p1 | p2 +
#     plot_annotation(
#       title = paste("Inbreeding (F =", round(F,2), ")"),
#       theme = theme(plot.title = element_text(face="bold"))
#     )
# }
# 
# # ============================================================
# # 5) ADMIXTURE / WAHLUND EFFECT WRAPPER
# # ============================================================
# 
# plot_admixture_wahlund <- function(p1 = NULL,
#                                    p2 = NULL,
#                                    geno1 = NULL,
#                                    geno2 = NULL,
#                                    m = 0.5,
#                                    N_theoretical = 1000,
#                                    base_size = 18) {
#   
#   required <- c("AA","Aa","aa")
#   
#   clean_counts <- function(g) {
#     g[setdiff(required, names(g))] <- 0
#     g[required]
#   }
#   
#   get_p <- function(g) {
#     N <- sum(g)
#     nA <- 2*g["AA"] + g["Aa"]
#     as.numeric(nA)/(2*N)
#   }
#   
#   geno_hwe_fun <- function(p){
#     q <- 1 - p
#     c(AA = p^2, Aa = 2*p*q, aa = q^2)
#   }
#   
#   # -------------------------------------------------
#   # Case 1: genotype counts provided
#   # -------------------------------------------------
#   
#   if (!is.null(geno1) && !is.null(geno2)) {
#     
#     geno1 <- clean_counts(geno1)
#     geno2 <- clean_counts(geno2)
#     
#     geno_pool <- geno1 + geno2
#     N <- sum(geno_pool)
#     
#     p1 <- get_p(geno1)
#     p2 <- get_p(geno2)
#     pbar <- get_p(geno_pool)
#     
#     g_pool_obs <- geno_pool / N
#     g_pool_hwe <- geno_hwe_fun(pbar)
#     
#     expected_counts <- g_pool_hwe * N
#     chisq <- sum((geno_pool - expected_counts)^2 / expected_counts)
#     pval  <- pchisq(chisq, df = 1, lower.tail = FALSE)
#     
#     subtitle_text <- paste0(
#       "N = ", N,
#       " | Estimated p = ", round(pbar,3),
#       " | χ² = ", round(chisq,2),
#       " | P-value = ", signif(pval,3)
#     )
#     
#     # -------------------------------------------------
#     # Case 2: allele frequencies provided
#     # -------------------------------------------------
#     
#   } else if (!is.null(p1) && !is.null(p2)) {
#     
#     # pooled allele frequency
#     pbar <- m*p1 + (1-m)*p2
#     
#     g1 <- geno_hwe_fun(p1)
#     g2 <- geno_hwe_fun(p2)
#     
#     g_pool_obs <- m*g1 + (1-m)*g2
#     g_pool_hwe <- geno_hwe_fun(pbar)
#     
#     # Use theoretical N for chi-square scaling
#     N <- N_theoretical
#     
#     obs_counts <- g_pool_obs * N
#     expected_counts <- g_pool_hwe * N
#     
#     chisq <- sum((obs_counts - expected_counts)^2 / expected_counts)
#     pval  <- pchisq(chisq, df = 1, lower.tail = FALSE)
#     
#     subtitle_text <- paste0(
#       "Theoretical N = ", N,
#       " | Pooled p = ", round(pbar,3),
#       " | χ² = ", round(chisq,2),
#       " | P-value = ", signif(pval,3)
#     )    
#   } else {
#     stop("Provide either (p1,p2) OR (geno1,geno2)")
#   }
#   
#   # -------------------------------------------------
#   # Data frames
#   # -------------------------------------------------
#   
#   allele_df <- tibble::tibble(
#     pop = c("Pop 1","Pop 2","Pooled"),
#     A = c(p1,p2,pbar),
#     a = 1 - A
#   ) |>
#     tidyr::pivot_longer(cols = c(A,a),
#                         names_to = "allele",
#                         values_to = "freq")
#   
#   geno_df <- tibble::tibble(
#     genotype = required,
#     Observed = as.numeric(g_pool_obs),
#     `Expected (HWE)` = as.numeric(g_pool_hwe)
#   ) |>
#     tidyr::pivot_longer(cols = -genotype,
#                         names_to = "type",
#                         values_to = "freq")
#   
#   # -------------------------------------------------
#   # Plots
#   # -------------------------------------------------
#   
#   pA_plot <- ggplot2::ggplot(allele_df,
#                              ggplot2::aes(pop, freq, fill = allele)) +
#     ggplot2::geom_col(position = ggplot2::position_dodge(0.7),
#                       width = 0.6,
#                       color="black") +
#     ggplot2::scale_fill_manual(values = bio_cols$allele) +
#     ggplot2::coord_cartesian(ylim=c(0,1)) +
#     ggplot2::labs(title="Allele frequencies",
#                   y="Frequency", x=NULL) +
#     theme_bio(base_size)
#   
#   pG_plot <- ggplot2::ggplot(geno_df,
#                              ggplot2::aes(genotype, freq, fill = type)) +
#     ggplot2::geom_col(position = ggplot2::position_dodge(0.7),
#                       width = 0.6,
#                       color="black") +
#     ggplot2::scale_fill_manual(values = bio_cols$model) +
#     ggplot2::coord_cartesian(ylim=c(0,1)) +
#     ggplot2::labs(
#       title="Pooled genotypes deviate from HWE",
#       subtitle = subtitle_text,
#       y="Frequency", x=NULL
#     ) +
#     theme_bio(base_size)
#   
#   pA_plot | pG_plot +
#     patchwork::plot_annotation(
#       title="Admixture generates Hardy–Weinberg deviation (Wahlund effect)"
#     )
# }
# 
# # ============================================================
# # END OF EXTENDED FRAMEWORK
# # ============================================================
# 
# 
# 
# generate_population <- function(N = NULL,
#                                 p = NULL,
#                                 geno_counts = NULL,
#                                 seed = NULL) {
#   
#   if (!is.null(seed)) set.seed(seed)
#   
#   # ------------------------------------------------------------
#   # MODE 1: Genotype counts provided
#   # ------------------------------------------------------------
#   
#   if (!is.null(geno_counts)) {
#     
#     # Ensure named vector
#     stopifnot(all(c("AA","Aa","aa") %in% names(geno_counts)))
#     
#     N <- sum(geno_counts)
#     
#     # Expand genotypes
#     genotypes <- rep(names(geno_counts), geno_counts)
#     
#     # Convert to allele vector
#     alleles <- unlist(strsplit(genotypes, ""))
#     
#     # Shuffle individuals
#     genotypes <- sample(genotypes)
#     alleles <- unlist(strsplit(genotypes, ""))
#     
#   } else {
#     
#     # ------------------------------------------------------------
#     # MODE 2: Frequency p provided
#     # ------------------------------------------------------------
#     
#     stopifnot(!is.null(N), !is.null(p))
#     
#     n_alleles <- 2 * N
#     
#     alleles <- sample(c(
#       rep("A", round(p * n_alleles)),
#       rep("a", n_alleles - round(p * n_alleles))
#     ))
#     
#     # Build temporary genotypes
#     genotypes <- paste0(
#       alleles[seq(1, length(alleles), 2)],
#       alleles[seq(2, length(alleles), 2)]
#     )
#     
#     genotypes <- ifelse(genotypes %in% c("Aa","aA"), "Aa", genotypes)
#   }
#   
#   # ------------------------------------------------------------
#   # Build tidy structure
#   # ------------------------------------------------------------
#   
#   geno_labels <- tibble(
#     id = 1:N,
#     genotype = genotypes
#   )
#   
#   geno <- geno_labels %>%
#     mutate(copy = 1) %>%
#     tidyr::uncount(2, .id = "copy") %>%
#     group_by(id) %>%
#     mutate(allele = strsplit(genotype[1], "")[[1]][copy]) %>%
#     ungroup()
#   
#   list(
#     geno = geno,
#     geno_labels = geno_labels
#   )
# }
# 
# 
# layout_population_grid <- function(geno) {
#   
#   N <- length(unique(geno$id))
#   grid_size <- ceiling(sqrt(N))
#   
#   geno$row <- ceiling(geno$id / grid_size)
#   geno$col <- (geno$id - 1) %% grid_size + 1
#   
#   geno$x <- geno$col + ifelse(geno$copy == 1, -0.10, 0.10)
#   geno$y <- -geno$row
#   
#   geno
# }
# 
# plot_population_individuals <- function(geno,
#                                         geno_labels,
#                                         base_size = 20) {
#   
#   label_df <- geno_labels %>%
#     mutate(
#       row = ceiling(id / ceiling(sqrt(n()))),
#       col = (id - 1) %% ceiling(sqrt(n())) + 1,
#       x = col,
#       y = -row + 0.30
#     )
#   
#   ggplot() +
#     geom_point(data = geno,
#                aes(x, y, fill = allele),
#                shape = 21,
#                size = 6,
#                color = "black",
#                stroke = 0.7) +
#     geom_text(data = label_df,
#               aes(x, y, label = genotype),
#               size = 4,
#               color = "grey20") +
#     scale_fill_manual(values = bio_cols$allele) +
#     coord_equal() +
#     labs(title = "Diploid individuals") +
#     theme_void(base_size = base_size) +
#     theme(legend.position = "right")
# }
# 
# plot_population_counts <- function(geno_labels,
#                                    base_size = 20) {
#   
#   geno_counts <- geno_labels %>%
#     count(genotype)
#   
#   plot_bio_bar(
#     geno_counts,
#     genotype,
#     n,
#     genotype,
#     palette = bio_cols$genotype,
#     ylim = c(0, max(geno_counts$n) * 1.1),
#     title = "Genotype counts",
#     ylab = "Count",
#     base_size = base_size
#   )
# }
# 
# 
# 
# plot_sequence <- function(seqs,
#                           highlight_regions = NULL,
#                           show_seg_sites = TRUE,
#                           base_size = 16,
#                           title = "Sequence variation across homologous chromosomes") {
#   
#   # ------------------------------------------------------------
#   # Basic checks
#   # ------------------------------------------------------------
#   
#   stopifnot(is.character(seqs))
#   stopifnot(length(unique(nchar(seqs))) == 1)
#   
#   n_seq <- length(seqs)
#   seq_length <- nchar(seqs[1])
#   
#   # ------------------------------------------------------------
#   # Convert to tidy format
#   # ------------------------------------------------------------
#   
#   df <- tibble(
#     id = factor(seq_along(seqs)),
#     seq = seqs
#   ) %>%
#     mutate(
#       pos  = purrr::map(seq, ~ seq_along(strsplit(.x, "")[[1]])),
#       base = purrr::map(seq, ~ strsplit(.x, "")[[1]])
#     ) %>%
#     tidyr::unnest(c(pos, base))
#   
#   # ------------------------------------------------------------
#   # Detect segregating sites
#   # ------------------------------------------------------------
#   
#   seg_sites <- df %>%
#     group_by(pos) %>%
#     summarise(n = n_distinct(base), .groups = "drop") %>%
#     filter(n > 1) %>%
#     pull(pos)
#   
#   df <- df %>%
#     mutate(is_seg = pos %in% seg_sites)
#   
#   # ------------------------------------------------------------
#   # Base color palette (color-blind friendly)
#   # ------------------------------------------------------------
#   
#   base_cols <- c(
#     "A" = "#D55E00",
#     "T" = "#009E73",
#     "C" = "#0072B2",
#     "G" = "#CC79A7",
#     "-" = "grey60"
#   )
#   
#   # ------------------------------------------------------------
#   # Base plot
#   # ------------------------------------------------------------
#   
#   p <- ggplot(df, aes(pos, forcats::fct_rev(id))) +
#     
#     geom_text(
#       aes(label = base,
#           color = base,
#           fontface = ifelse(show_seg_sites & is_seg, "bold", "plain")),
#       family = "mono",
#       size = 4.5
#     ) +
#     
#     scale_color_manual(values = base_cols) +
#     
#     scale_x_continuous(
#       breaks = seq(0, seq_length, 5),
#       limits = c(0.5, seq_length + 0.5),
#       expand = c(0, 0)
#     ) +
#     
#     scale_y_discrete(
#       limits = rev(levels(df$id)),
#       expand = expansion(mult = c(0, 0))
#     ) +
#     
#     coord_fixed(ratio = 0.6, clip = "off") +
#     
#     labs(
#       x = "Genomic position",
#       y = NULL,
#       title = title
#     ) +
#     
#     theme_bio(base_size = base_size) +
#     theme(
#       legend.position = "none",
#       panel.grid = element_blank(),
#       axis.text.y = element_text(size = base_size * 0.7)
#     )
#   
#   # ------------------------------------------------------------
#   # Optional highlighted regions
#   # ------------------------------------------------------------
#   
#   if (!is.null(highlight_regions)) {
#     
#     for (region in highlight_regions) {
#       
#       p <- p +
#         annotate("rect",
#                  xmin = region[1] - 0.5,
#                  xmax = region[2] + 0.5,
#                  ymin = 0.5,
#                  ymax = n_seq + 0.5,
#                  fill = NA,
#                  color = "black",
#                  linewidth = 0.7)
#     }
#   }
#   print(p)
#   invisible(p)
# }
# 
# 
# plot_inbreeding_small_pop <- function(
#     N_values = c(20, 50, 200),
#     generations = 100,
#     base_size = 18
# ) {
#   
#   library(dplyr)
#   library(ggplot2)
#   
#   df <- expand.grid(
#     generation = 0:generations,
#     N = N_values
#   ) %>%
#     arrange(N, generation) %>%
#     mutate(
#       F = 1 - (1 - 1/(2*N))^generation,
#       N_label = paste0("N = ", N)
#     )
#   
#   ggplot(df, aes(generation, F, color = N_label)) +
#     geom_line(linewidth = 1.4) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = "Inbreeding increases faster in small populations",
#       x = "Generation",
#       y = "Inbreeding coefficient (F)",
#       color = NULL
#     ) +
#     theme_bio(base_size = base_size)
# }
# 
# 
# plot_dominant_vs_recessive_selection <- function(
#     p0 = 0.1,
#     generations = 100,
#     s = 0.2,
#     base_size = 18
# ) {
#   
#   library(dplyr)
#   library(ggplot2)
#   
#   # -----------------------------------------
#   # Selection recursion
#   # -----------------------------------------
#   
#   step_selection <- function(p, wAA, wAa, waa) {
#     q <- 1 - p
#     wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
#     (p^2*wAA + p*q*wAa) / wbar
#   }
#   
#   simulate <- function(p0, generations, wAA, wAa, waa) {
#     p <- numeric(generations + 1)
#     p[1] <- p0
#     for (t in 1:generations)
#       p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
#     p
#   }
#   
#   # -----------------------------------------
#   # Case 1: Favored dominant allele A
#   # -----------------------------------------
#   
#   p_dom <- simulate(
#     p0,
#     generations,
#     wAA = 1 + s,
#     wAa = 1 + s,
#     waa = 1
#   )
#   
#   # -----------------------------------------
#   # Case 2: Favored recessive allele A
#   # -----------------------------------------
#   
#   p_rec <- simulate(
#     p0,
#     generations,
#     wAA = 1 + s,
#     wAa = 1,
#     waa = 1
#   )
#   
#   df <- data.frame(
#     generation = rep(0:generations, 2),
#     p = c(p_dom, p_rec),
#     type = rep(
#       c("Dominant allele",
#         "Recessive allele"),
#       each = generations + 1
#     )
#   )
#   
#   # -----------------------------------------
#   # Plot
#   # -----------------------------------------
#   
#   ggplot(df, aes(generation, p, color = type)) +
#     geom_line(linewidth = 1.4) +
#     coord_cartesian(ylim = c(0,1)) +
#     scale_color_manual(
#       values = c(
#         "Dominant allele" = "#D55E00",  # red/orange
#         "Recessive allele" = "#0072B2"  # blue
#       )
#     ) +
#     labs(
#       title = "Change in allele frequency under natural selection",
#       x = "Generation",
#       y = "Allele frequency (p)",
#       color = NULL
#     ) +
#     theme_bio(base_size = base_size)
# }
# 
# 
# 
# library(ggplot2)
# library(dplyr)
# 
# plot_ibd_pedigree <- function() {
#   
#   # ------------------------------------------------------------
#   # Node positions
#   # ------------------------------------------------------------
#   
#   nodes <- tibble(
#     id = c("A", "B", "C", "I"),
#     x  = c(0, -2, 2, 0),
#     y  = c(4, 2, 2, 0)
#   )
#   
#   # ------------------------------------------------------------
#   # Edges (parent-child connections)
#   # ------------------------------------------------------------
#   
#   edges <- tibble(
#     x    = c(0, 0, -2, 2),
#     y    = c(4, 4, 2, 2),
#     xend = c(-2, 2, 0, 0),
#     yend = c(2, 2, 0, 0)
#   )
#   
#   # ------------------------------------------------------------
#   # Plot
#   # ------------------------------------------------------------
#   
#   ggplot() +
#     
#     # Parent-child lines
#     geom_segment(data = edges,
#                  aes(x = x, y = y,
#                      xend = xend, yend = yend),
#                  linewidth = 1) +
#     
#     # Individuals
#     geom_point(data = nodes,
#                aes(x = x, y = y),
#                size = 12,
#                shape = 21,
#                fill = "white",
#                color = "black",
#                stroke = 1.2) +
#     
#     # Labels
#     geom_text(data = nodes,
#               aes(x = x, y = y, label = id),
#               size = 6) +
#     
#     coord_fixed() +
#     xlim(-3, 3) +
#     ylim(-1, 5) +
#     theme_void(base_size = 18) +
#     labs(title = "Pedigree illustrating identity by descent")
# }
# 
# 
# plot_ibd_with_alleles <- function() {
#   
#   nodes <- tibble(
#     id = c("A", "B", "C", "I"),
#     x  = c(0, -2, 2, 0),
#     y  = c(4, 2, 2, 0)
#   )
#   
#   edges <- tibble(
#     x    = c(0, 0, -2, 2),
#     y    = c(4, 4, 2, 2),
#     xend = c(-2, 2, 0, 0),
#     yend = c(2, 2, 0, 0)
#   )
#   
#   # Allele positions inside each circle
#   alleles <- tibble(
#     id = rep(c("A","B","C","I"), each = 2),
#     dx = rep(c(-0.3, 0.3), 4),
#     dy = 0,
#     allele_color = c(
#       "#D55E00", "#0072B2",   # A
#       "#0072B2", "#D55E00",   # B
#       "#0072B2", "#0072B2",   # C
#       "#D55E00", "#D55E00"    # I (identical by descent)
#     )
#   ) %>%
#     left_join(nodes, by = "id") %>%
#     mutate(ax = x + dx,
#            ay = y + dy)
#   
#   ggplot() +
#     
#     geom_segment(data = edges,
#                  aes(x, y, xend = xend, yend = yend),
#                  linewidth = 1) +
#     
#     geom_point(data = nodes,
#                aes(x, y),
#                size = 12,
#                shape = 21,
#                fill = "white",
#                color = "black",
#                stroke = 1.2) +
#     
#     geom_point(data = alleles,
#                aes(ax, ay),
#                size = 5,
#                shape = 21,
#                fill = alleles$allele_color,
#                color = "black",
#                stroke = 0.7) +
#     
#     geom_text(data = nodes,
#               aes(x, y + 0.8, label = id),
#               size = 5) +
#     
#     coord_fixed() +
#     xlim(-3, 3) +
#     ylim(-1, 5) +
#     theme_void(base_size = 18) +
#     labs(title = "Identity by descent in a pedigree")
# }
# 
# 
# plot_pedigree <- function(pedigree,
#                           highlight_path = NULL,
#                           base_size = 18) {
#   
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   
#   # ---------------------------------------------
#   # Compute generation depth recursively
#   # ---------------------------------------------
#   
#   pedigree$generation <- NA
#   
#   # Founders = generation 1
#   pedigree$generation[is.na(pedigree$father) &
#                         is.na(pedigree$mother)] <- 1
#   
#   # Simple propagation (works for shallow pedigrees)
#   repeat {
#     updated <- FALSE
#     
#     for (i in seq_len(nrow(pedigree))) {
#       if (is.na(pedigree$generation[i])) {
#         
#         parents <- pedigree %>%
#           filter(id %in% c(pedigree$father[i],
#                            pedigree$mother[i]))
#         
#         if (all(!is.na(parents$generation))) {
#           pedigree$generation[i] <- max(parents$generation) + 1
#           updated <- TRUE
#         }
#       }
#     }
#     
#     if (!updated) break
#   }
#   
#   # ---------------------------------------------
#   # Assign x positions within generations
#   # ---------------------------------------------
#   
#   nodes <- pedigree %>%
#     arrange(generation, id) %>%
#     group_by(generation) %>%
#     mutate(
#       x = seq_along(id),
#       y = -generation
#     ) %>%
#     ungroup() %>%
#     select(id, x, y)
#   
#   # ---------------------------------------------
#   # Build edges
#   # ---------------------------------------------
#   
#   edges <- pedigree %>%
#     pivot_longer(cols = c(father, mother),
#                  names_to = "parent_type",
#                  values_to = "parent") %>%
#     filter(!is.na(parent)) %>%
#     select(parent, child = id)
#   
#   edge_coords <- edges %>%
#     left_join(nodes %>% rename(parent = id,
#                                x = x,
#                                y = y),
#               by = "parent") %>%
#     left_join(nodes %>% rename(child = id,
#                                xend = x,
#                                yend = y),
#               by = "child")
#   
#   # ---------------------------------------------
#   # Highlight path
#   # ---------------------------------------------
#   
#   if (!is.null(highlight_path)) {
#     highlight_edges <- tibble(
#       parent = highlight_path[-length(highlight_path)],
#       child  = highlight_path[-1]
#     )
#     
#     edge_coords <- edge_coords %>%
#       mutate(is_highlight =
#                paste(parent, child) %in%
#                paste(highlight_edges$parent,
#                      highlight_edges$child))
#   } else {
#     edge_coords$is_highlight <- FALSE
#   }
#   
#   # ---------------------------------------------
#   # Plot
#   # ---------------------------------------------
#   
#   ggplot() +
#     geom_segment(data = edge_coords,
#                  aes(x = x, y = y,
#                      xend = xend, yend = yend,
#                      linewidth = is_highlight),
#                  color = "black") +
#     scale_linewidth_manual(values = c("TRUE" = 2,
#                                       "FALSE" = 0.8),
#                            guide = "none") +
#     geom_point(data = nodes,
#                aes(x, y),
#                size = 14,
#                shape = 21,
#                fill = "white",
#                color = "black",
#                stroke = 1.2) +
#     geom_text(data = nodes,
#               aes(x, y + 0.5, label = id),
#               size = base_size * 0.3) +
#     coord_fixed() +
#     theme_void(base_size = base_size)
# }
# 
# 
# plot_pedigree_alleles <- function(ped,
#                                   alleles = NULL,
#                                   mode = c("dots", "letters", "both"),
#                                   allele_cols = c("A" = "#D55E00",
#                                                   "a" = "#0072B2"),
#                                   dx_mult = 0.25,
#                                   dy_mult = -0.35,   # negative = upward in kinship2 plot coords
#                                   cex_dots = 1.4,
#                                   cex_letters = 0.9,
#                                   letters_col = "black",
#                                   letters_font = 2,
#                                   dot_border_col = "black",
#                                   dot_border_lwd = 0.8,
#                                   mar = c(4, 4, 4, 4),
#                                   xpd = NA,
#                                   ...) {
#   
#   mode <- match.arg(mode)
#   
#   # Basic checks
#   if (!inherits(ped, "pedigree")) {
#     stop("`ped` must be a kinship2::pedigree object.")
#   }
#   
#   if (!is.null(alleles)) {
#     # Expect a named list: alleles[["id"]] = c("A","a")
#     if (!is.list(alleles)) stop("`alleles` must be a list or NULL.")
#     # Soft check: ensure each provided entry has length 2
#     bad <- names(alleles)[vapply(alleles, length, integer(1)) != 2]
#     if (length(bad) > 0) {
#       stop("Each element of `alleles` must be a length-2 vector. Problem IDs: ",
#            paste(bad, collapse = ", "))
#     }
#   }
#   
#   # Set plotting params (temporarily)
#   op <- par(no.readonly = TRUE)
#   on.exit(par(op), add = TRUE)
#   
#   par(mar = mar)
#   par(xpd = xpd)
#   
#   # Plot pedigree; capture coordinates returned by plot.pedigree()
#   coords <- plot(ped, ...)
#   
#   # If no alleles supplied, we're done
#   if (is.null(alleles)) {
#     return(invisible(coords))
#   }
#   
#   # Extract positions
#   x <- coords$x
#   y <- coords$y
#   ids <- as.character(ped$id)
#   
#   # Offsets based on node box
#   dx <- coords$boxw * dx_mult
#   dy <- coords$boxh * dy_mult
#   
#   # Helper: safe allele fetch
#   get_pair <- function(id) {
#     pair <- alleles[[id]]
#     if (is.null(pair)) return(c(NA_character_, NA_character_))
#     as.character(pair)
#   }
#   
#   for (i in seq_along(ids)) {
#     pair <- get_pair(ids[i])
#     a1 <- pair[1]; a2 <- pair[2]
#     
#     # Skip if missing
#     if (is.na(a1) || is.na(a2)) next
#     
#     # --- Dots overlay ---
#     if (mode %in% c("dots", "both")) {
#       points(x[i] - dx, y[i] + dy,
#              pch = 21,
#              bg  = allele_cols[a1],
#              col = dot_border_col,
#              lwd = dot_border_lwd,
#              cex = cex_dots)
#       
#       points(x[i] + dx, y[i] + dy,
#              pch = 21,
#              bg  = allele_cols[a2],
#              col = dot_border_col,
#              lwd = dot_border_lwd,
#              cex = cex_dots)
#     }
#     
#     # --- Letters overlay ---
#     if (mode %in% c("letters", "both")) {
#       text(x[i] - dx, y[i] + dy,
#            labels = a1,
#            cex = cex_letters,
#            col = letters_col,
#            font = letters_font)
#       
#       text(x[i] + dx, y[i] + dy,
#            labels = a2,
#            cex = cex_letters,
#            col = letters_col,
#            font = letters_font)
#     }
#   }
#   
#   invisible(coords)
# }
# 
# 
# 
# # Selective sweep
# 
# library(ggplot2)
# 
# # =========================================================
# # TEXTBOOK HARD SWEEP (clean, dramatic)
# # =========================================================
# generate_textbook_sweep <- function(
#     n_hap = 100,
#     L = 201,
#     sweep_freq = 1.0,
#     core_width = 30,
#     seed = NULL
# ) {
#   if (!is.null(seed)) set.seed(seed)
#   stopifnot(L %% 2 == 1)
#   
#   center <- (L + 1) / 2
#   pos <- seq_len(L)
#   dist <- abs(pos - center)
#   
#   H <- matrix(0, nrow = n_hap, ncol = L)
#   
#   core <- dist <= core_width
#   is_sweep <- rbinom(n_hap, 1, sweep_freq) == 1
#   
#   for (i in seq_len(n_hap)) {
#     H[i, ] <- rbinom(L, 1, 0.5)
#     if (is_sweep[i]) {
#       H[i, core] <- 1
#       H[i, center] <- 1
#     } else {
#       H[i, center] <- 0
#     }
#   }
#   
#   list(H = H, center = center, pos = pos, dist = dist)
# }
# 
# # =========================================================
# # REALISTIC HARD SWEEP (noisy)
# # =========================================================
# generate_realistic_sweep <- function(
#     n_hap = 100,
#     L = 201,
#     sweep_freq = 0.9,
#     rec_scale = 30,
#     seed = NULL
# ) {
#   if (!is.null(seed)) set.seed(seed)
#   stopifnot(L %% 2 == 1)
#   
#   center <- (L + 1) / 2
#   pos <- seq_len(L)
#   dist <- abs(pos - center)
#   
#   base_p <- runif(L, 0.2, 0.8)
#   copy_prob <- exp(-dist / rec_scale)
#   copy_prob[center] <- 1
#   
#   H <- matrix(0, nrow = n_hap, ncol = L)
#   is_sweep <- rbinom(n_hap, 1, sweep_freq) == 1
#   
#   for (i in seq_len(n_hap)) {
#     H[i, ] <- rbinom(L, 1, base_p)
#     if (is_sweep[i]) {
#       match <- rbinom(L, 1, copy_prob) == 1
#       H[i, match] <- 1
#       H[i, center] <- 1
#     } else {
#       H[i, center] <- 0
#     }
#   }
#   
#   list(H = H, center = center, pos = pos, dist = dist)
# }
# 
# # =========================================================
# # STATISTICS (same as Shiny app)
# # =========================================================
# gene_diversity <- function(H) {
#   p <- colMeans(H)
#   2 * p * (1 - p)
# }
# 
# ld_r2_with_center <- function(H, center) {
#   x <- H[, center]
#   px <- mean(x)
#   vx <- px * (1 - px)
#   r2 <- numeric(ncol(H))
#   
#   for (j in seq_len(ncol(H))) {
#     y <- H[, j]
#     py <- mean(y)
#     vy <- py * (1 - py)
#     if (vx == 0 || vy == 0) {
#       r2[j] <- NA_real_
#     } else {
#       cov_xy <- mean((x - px) * (y - py))
#       r <- cov_xy / sqrt(vx * vy)
#       r2[j] <- r^2
#     }
#   }
#   r2
# }
# 
# # =========================================================
# # EXACT SAME PLOTS AS SHINY (extracted)
# # =========================================================
# make_sweep_plots <- function(
#     mode = c("Textbook Hard Sweep", "Realistic Hard Sweep"),
#     n_hap = 100,
#     L = 201,
#     sweep_freq = 1.0,
#     width = 30,          # core_width OR rec_scale
#     sort_haps = TRUE,
#     seed = NULL
# ) {
#   mode <- match.arg(mode)
#   
#   L <- as.integer(L)
#   if (L %% 2 == 0) L <- L + 1
#   
#   sweep_data <- if (mode == "Textbook Hard Sweep") {
#     generate_textbook_sweep(
#       n_hap = n_hap, L = L, sweep_freq = sweep_freq,
#       core_width = width, seed = seed
#     )
#   } else {
#     generate_realistic_sweep(
#       n_hap = n_hap, L = L, sweep_freq = sweep_freq,
#       rec_scale = width, seed = seed
#     )
#   }
#   
#   H <- sweep_data$H
#   center <- sweep_data$center
#   
#   # Optional sorting to make sweep block visually clean (same as Shiny)
#   if (sort_haps) {
#     H <- H[order(H[, center], decreasing = TRUE), , drop = FALSE]
#   }
#   
#   # ---- Haplotypes plot (same as Shiny) ----
#   df_hap <- data.frame(
#     Hap = rep(seq_len(nrow(H)), each = ncol(H)),
#     Pos = rep(seq_len(ncol(H)), times = nrow(H)),
#     Allele = as.vector(H)
#   )
#   df_hap$Allele <- as.vector(t(H))  # same final assignment as in Shiny app
#   
#   hapPlot <- ggplot(df_hap, aes(x = Pos, y = Hap, fill = factor(Allele))) +
#     geom_tile() +
#     geom_vline(xintercept = center, linetype = "dashed") +
#     theme_minimal() +
#     labs(fill = "Allele",
#          title = "Haplotype structure around selected site")
#   
#   # ---- Gene diversity plot (same as Shiny) ----
#   dist <- abs(seq_len(ncol(H)) - center)
#   gd <- gene_diversity(H)
#   df_gd <- data.frame(Distance = dist, GD = gd)
#   
#   gdPlot <- ggplot(df_gd, aes(x = Distance, y = GD)) +
#     geom_line(linewidth = 1.2) +
#     theme_minimal() +
#     ylim(0, 0.5) +
#     labs(title = "Diversity trough around selected site",
#          y = "2p(1-p)")
#   
#   # ---- LD plot (same as Shiny) ----
#   r2 <- ld_r2_with_center(H, center)
#   df_ld <- data.frame(Distance = dist, r2 = r2)
#   
#   ldPlot <- ggplot(df_ld, aes(x = Distance, y = r2)) +
#     geom_line(linewidth = 1.2) +
#     theme_minimal() +
#     ylim(0, 1) +
#     labs(title = "LD peak near selected site",
#          y = expression(r^2))
#   
#   list(hapPlot = hapPlot, gdPlot = gdPlot, ldPlot = ldPlot)
# }
# 
# 
# plot_mutation_dynamics <- function(p0 = 0.9,
#                                    mu = 0.001,
#                                    nu = 0.0001,
#                                    generations = 200,
#                                    base_size = 18) {
#   
#   p <- numeric(generations + 1)
#   p[1] <- p0
#   
#   for(t in 1:generations){
#     p[t+1] <- p[t]*(1-mu) + (1-p[t])*nu
#   }
#   
#   df <- tibble(
#     gen = 0:generations,
#     p = p,
#     H = 2*p*(1-p)
#   )
#   
#   p_eq <- nu / (mu + nu)
#   
#   ggplot(df, aes(gen, p)) +
#     geom_line(linewidth = 1.5, color = "#0072B2") +
#     geom_hline(yintercept = p_eq,
#                linetype = "dashed") +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = "Mutation drives allele frequency toward equilibrium",
#       subtitle = paste0("Equilibrium p* = ", round(p_eq,4)),
#       x = "Generation",
#       y = "Allele frequency (p)"
#     ) +
#     theme_bio(base_size)
# }
# 
# plot_migration_dynamics <- function(p1 = 0.2,
#                                     p2 = 0.8,
#                                     m = 0.05,
#                                     generations = 100,
#                                     base_size = 18) {
#   
#   p1_traj <- numeric(generations + 1)
#   p2_traj <- numeric(generations + 1)
#   
#   p1_traj[1] <- p1
#   p2_traj[1] <- p2
#   
#   for(t in 1:generations){
#     p1_traj[t+1] <- (1-m)*p1_traj[t] + m*p2_traj[t]
#     p2_traj[t+1] <- (1-m)*p2_traj[t] + m*p1_traj[t]
#   }
#   
#   df <- tibble(
#     gen = rep(0:generations, 2),
#     p = c(p1_traj, p2_traj),
#     pop = rep(c("Population 1", "Population 2"),
#               each = generations + 1)
#   )
#   
#   ggplot(df, aes(gen, p, color = pop)) +
#     geom_line(linewidth = 1.5) +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = "Migration causes allele frequencies to converge",
#       subtitle = paste0("Migration rate m = ", m),
#       x = "Generation",
#       y = "Allele frequency (p)"
#     ) +
#     theme_bio(base_size)
# }
# 
# plot_effective_population_size <- function(Nm = 10,
#                                            Nf = 90,
#                                            base_size = 18) {
#   
#   N <- Nm + Nf
#   Ne <- (4*Nm*Nf)/(Nm + Nf)
#   
#   df <- tibble(
#     Type = c("Census size (N)", "Effective size (Ne)"),
#     Size = c(N, Ne)
#   )
#   
#   ggplot(df, aes(Type, Size, fill = Type)) +
#     geom_col(width = 0.6, color = "black") +
#     labs(
#       title = "Effective population size can be much smaller than N",
#       subtitle = paste0("Nm = ", Nm, ", Nf = ", Nf),
#       y = "Population size",
#       x = NULL
#     ) +
#     theme_bio(base_size) +
#     theme(legend.position = "none")
# }
# 
# 
# plot_heterozygosity_decay <- function(N = 50,
#                                       generations = 100,
#                                       H0 = 0.5,
#                                       base_size = 18) {
#   
#   H <- numeric(generations + 1)
#   H[1] <- H0
#   
#   for(t in 1:generations){
#     H[t+1] <- H[t] * (1 - 1/(2*N))
#   }
#   
#   df <- tibble(
#     gen = 0:generations,
#     H = H
#   )
#   
#   ggplot(df, aes(gen, H)) +
#     geom_line(linewidth = 1.5, color = "#D55E00") +
#     labs(
#       title = "Heterozygosity declines under genetic drift",
#       subtitle = paste0("Population size N = ", N),
#       x = "Generation",
#       y = "Expected heterozygosity (H)"
#     ) +
#     theme_bio(base_size)
# }
# 
# 
# plot_molecular_clock <- function(mu = 1e-8,
#                                  L = 1e6,
#                                  generations = 300,
#                                  stochastic = FALSE,
#                                  base_size = 18,
#                                  seed = 1) {
#   
#   if (stochastic) set.seed(seed)
#   
#   t <- 0:generations
#   
#   if (!stochastic) {
#     
#     divergence <- 2 * mu * L * t
#     
#   } else {
#     
#     new_mut <- rpois(generations, lambda = 2 * mu * L)
#     divergence <- c(0, cumsum(new_mut))
#   }
#   
#   df <- tibble(
#     generation = t,
#     divergence = divergence
#   )
#   
#   ggplot(df, aes(generation, divergence)) +
#     geom_line(linewidth = 1.5, color = "#0072B2") +
#     labs(
#       title = "Molecular clock: divergence increases over time",
#       subtitle = paste0("μ = ", mu, ", L = ", format(L, scientific=TRUE)),
#       x = "Time (generations)",
#       y = "Number of substitutions"
#     ) +
#     theme_bio(base_size)
# }
# 
# 
# plot_domestication_bottleneck <- function(
#     N_ancestral = 1000,
#     N_bottleneck = 20,
#     N_post = 200,
#     bottleneck_start = 40,
#     bottleneck_duration = 20,
#     generations = 120,
#     H0 = 0.5,
#     base_size = 18
# ) {
#   
#   H <- numeric(generations + 1)
#   H[1] <- H0
#   
#   for(t in 1:generations){
#     
#     if (t < bottleneck_start) {
#       N_current <- N_ancestral
#     } else if (t < bottleneck_start + bottleneck_duration) {
#       N_current <- N_bottleneck
#     } else {
#       N_current <- N_post
#     }
#     
#     H[t+1] <- H[t] * (1 - 1/(2*N_current))
#   }
#   
#   df <- tibble(
#     generation = 0:generations,
#     heterozygosity = H
#   )
#   
#   ggplot(df, aes(generation, heterozygosity)) +
#     geom_line(linewidth = 1.5, color = "#D55E00") +
#     geom_vline(xintercept = bottleneck_start,
#                linetype = "dashed") +
#     geom_vline(xintercept = bottleneck_start + bottleneck_duration,
#                linetype = "dashed") +
#     labs(
#       title = "Domestication bottleneck reduces genetic diversity",
#       subtitle = paste0(
#         "Ancestral N=", N_ancestral,
#         " → Bottleneck N=", N_bottleneck
#       ),
#       x = "Generation",
#       y = "Expected heterozygosity"
#     ) +
#     theme_bio(base_size)
# }
# 
# plot_coalescent_tree <- function(n = 10,
#                                  Ne = 100,
#                                  base_size = 16,
#                                  seed = 1,
#                                  max_height = NULL) {
#   
#   library(ape)
#   set.seed(seed)
#   
#   tree <- rcoal(n)
#   tree$edge.length <- tree$edge.length * 2 * Ne
#   
#   tree_height <- max(node.depth.edgelength(tree))
#   
#   # Determine common scaling height
#   if (is.null(max_height)) {
#     max_height <- tree_height
#   }
#   
#   # Scale tree to desired height
#   scaling_factor <- max_height / tree_height
#   tree$edge.length <- tree$edge.length * scaling_factor
#   
#   plot(tree,
#        direction = "upwards",
#        show.tip.label = TRUE,
#        cex = 0.9)
#   
#   title(main = paste0("Coalescent tree (Ne = ", Ne, ")"))
#   mtext("Time (relative scale)",
#         side = 2,
#         line = 2)
# }
# 
# plot_bottleneck_coalescent <- function(n = 10,
#                                        Ne_ancestral = 1000,
#                                        Ne_recent = 20,
#                                        base_size = 16,
#                                        seed = 1) {
#   
#   library(ape)
#   set.seed(seed)
#   
#   tree <- rcoal(n)
#   
#   # Artificially compress recent branches
#   tree$edge.length <- tree$edge.length * 2 * Ne_ancestral
#   
#   # Shorten terminal branches (recent time)
#   terminal_edges <- which(tree$edge[,2] <= n)
#   tree$edge.length[terminal_edges] <-
#     tree$edge.length[terminal_edges] *
#     (Ne_recent / Ne_ancestral)
#   
#   plot(tree,
#        direction = "upwards",
#        show.tip.label = TRUE,
#        cex = 0.9)
#   
#   title(main = paste0("Bottleneck: Ne ",
#                       Ne_recent, " → ",
#                       Ne_ancestral),
#         cex.main = base_size * 0.07)
#   
#   mtext("Time (generations)",
#         side = 2,
#         line = 2,
#         cex = base_size * 0.06)
# }
# 
# 
# # 1) Simulate + scale tree in generations
# simulate_coalescent_tree <- function(n = 8, Ne = 100, seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#   tr <- rcoal(n)
#   tr$edge.length <- tr$edge.length * (2 * Ne)  # coalescent units -> generations
#   tr
# }
# 
# # 2) Plot tree with a fixed y-axis (time axis when direction="upwards")
# plot_coalescent_tree_fixed <- function(tree, main = "", y_max = NULL) {
#   if (is.null(y_max)) y_max <- max(node.depth.edgelength(tree))
#   
#   plot(tree,
#        direction = "upwards",
#        show.tip.label = TRUE,
#        cex = 0.9,
#        y.lim = c(0, y_max))   # <-- THIS is the key
#   
#   title(main = main)
#   mtext("Time (generations)", side = 2, line = 2)
# }
# 
# 
# # ============================================================
# # FINAL CORE POPULATION GENETICS FUNCTIONS
# # Fixation probability
# # FST vs migration
# # Mutation–selection balance
# # Site frequency spectrum
# # ============================================================
# 
# library(ggplot2)
# library(dplyr)
# 
# # ------------------------------------------------------------
# # 1️⃣ Fixation Probability Under Drift
# # ------------------------------------------------------------
# 
# plot_fixation_probability <- function(N = 50,
#                                       p0 = 0.3,
#                                       n_sim = 200,
#                                       base_size = 18,
#                                       seed = 1) {
#   
#   set.seed(seed)
#   
#   fix_count <- 0
#   
#   for (sim in 1:n_sim) {
#     
#     p <- p0
#     
#     while (p > 0 & p < 1) {
#       p <- rbinom(1, 2*N, p) / (2*N)
#     }
#     
#     if (p == 1) fix_count <- fix_count + 1
#   }
#   
#   prob_fix <- fix_count / n_sim
#   
#   df <- data.frame(
#     Outcome = c("Fixation", "Loss"),
#     Probability = c(prob_fix, 1 - prob_fix)
#   )
#   
#   ggplot(df, aes(Outcome, Probability, fill = Outcome)) +
#     geom_col(width = 0.6, color = "black") +
#     coord_cartesian(ylim = c(0,1)) +
#     labs(
#       title = "Neutral fixation probability",
#       subtitle = paste0("Start p = ", p0,
#                         " | Estimated P(fix) = ",
#                         round(prob_fix,3),
#                         " | Theoretical = ", p0),
#       y = "Probability",
#       x = NULL
#     ) +
#     theme_bio(base_size) +
#     theme(legend.position = "none")
# }
# 
# 
# # ------------------------------------------------------------
# # 2️⃣ FST vs Migration (Island Model)
# # ------------------------------------------------------------
# 
# plot_Fst_vs_migration <- function(Ne = 100,
#                                   m_values = seq(0.001, 0.5, length.out = 200),
#                                   base_size = 18) {
#   
#   Fst <- 1 / (4 * Ne * m_values + 1)
#   
#   df <- data.frame(
#     m = m_values,
#     Fst = Fst
#   )
#   
#   ggplot(df, aes(m, Fst)) +
#     geom_line(linewidth = 1.4, color = "#0072B2") +
#     labs(
#       title = "Population differentiation decreases with migration",
#       subtitle = paste0("Ne = ", Ne,
#                         " | FST = 1 / (4Nem + 1)"),
#       x = "Migration rate (m)",
#       y = expression(F[ST])
#     ) +
#     theme_bio(base_size)
# }
# 
# 
# # ------------------------------------------------------------
# # 3️⃣ Mutation–Selection Balance (Recessive deleterious)
# # ------------------------------------------------------------
# 
# plot_mutation_selection_balance <- function(mu = 1e-5,
#                                             s = 0.1,
#                                             base_size = 18) {
#   
#   q_hat <- sqrt(mu / s)
#   
#   q_vals <- seq(0, 0.05, length.out = 500)
#   selection_pressure <- mu - s * q_vals^2
#   
#   df <- data.frame(
#     q = q_vals,
#     balance = selection_pressure
#   )
#   
#   ggplot(df, aes(q, balance)) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_line(linewidth = 1.4, color = "#D55E00") +
#     geom_vline(xintercept = q_hat,
#                linetype = "dotted",
#                linewidth = 1) +
#     labs(
#       title = "Mutation–selection balance (recessive allele)",
#       subtitle = paste0("Equilibrium q̂ = sqrt(μ/s) = ",
#                         signif(q_hat,3)),
#       x = "Allele frequency (q)",
#       y = "Net change (mutation − selection)"
#     ) +
#     theme_bio(base_size)
# }
# 
# 
# # ------------------------------------------------------------
# # 4️⃣ Site Frequency Spectrum (Neutral expectation)
# # ------------------------------------------------------------
# 
# plot_site_frequency_spectrum <- function(n = 20,
#                                          base_size = 18) {
#   
#   k <- 1:(n-1)
#   SFS <- 1 / k
#   
#   df <- data.frame(
#     frequency_class = k,
#     proportion = SFS / sum(SFS)
#   )
#   
#   ggplot(df, aes(frequency_class, proportion)) +
#     geom_col(width = 0.7,
#              fill = "#4C78A8",
#              color = "black") +
#     labs(
#       title = "Neutral site frequency spectrum",
#       subtitle = "Expected SFS ∝ 1/k",
#       x = "Derived allele count (k)",
#       y = "Proportion of segregating sites"
#     ) +
#     theme_bio(base_size)
# }
# 
# 
# # ============================================================
# # MASTER POPULATION GENETICS CONCEPT MAP
# # ============================================================
# plot_population_genetics_concept_map <- function(base_size = 20) {
#   
#   library(ggplot2)
#   library(dplyr)
#   
#   nodes <- tibble::tribble(
#     ~label,                        ~x,  ~y,
#     
#     "Allele frequencies",          -4,  9,
#     "Genotype frequencies\n(HWE)",  4,  9,
#     
#     "Genetic drift",               -8,  6,
#     "Selection",                   -3,  6,
#     "Mutation",                     3,  6,
#     "Migration",                    8,  6,
#     
#     "Effective population size (Ne)", -8, 3.5,
#     "Inbreeding (F)",               -3, 3.5,
#     "Population structure (FST)",    3, 3.5,
#     "Wahlund effect",                8, 3.5,
#     
#     # Genealogy layer
#     "Coalescent tree",              -4, 0,
#     "Molecular clock",               0, 0,
#     "Site frequency spectrum",       6, 0,   # moved right
#     
#     # Genomics layer
#     "Selective sweeps",             -3, -3,
#     "Linkage disequilibrium",        6, -3    # moved right
#     )
#   
#   arrows <- tibble::tribble(
#     ~x, ~y, ~xend, ~yend,
#     
#     -4,9, -8,6,
#     -4,9, -3,6,
#     4,9,   3,6,
#     4,9,   8,6,
#     
#     -8,6, -8,3.5,
#     -8,3.5, -4,1,
#     
#     -3,6, -3,3.5,
#     
#     3,6,  0,1,
#     
#     8,6,  3,3.5,
#     3,3.5, 4,1,
#     
#     -4,1, -3,-2,
#     4,1,   3,-2
#   )
#   
#   ggplot() +
#     
#     # subtle horizontal layer guides (very light)
#     geom_hline(yintercept = c(7.5, 5, 2.5, 0),
#                color = "grey92",
#                linewidth = 0.6) +
#     
#     geom_segment(data = arrows,
#                  aes(x = x, y = y,
#                      xend = xend, yend = yend),
#                  arrow = arrow(length = unit(0.28, "cm")),
#                  linewidth = 0.5,
#                  color = "grey55") +
#     
#     geom_label(data = nodes,
#                aes(x, y, label = label),
#                fill = "white",
#                color = "black",
#                label.size = 0.4,
#                size = 5,
#                label.padding = unit(0.6, "lines"),
#                label.r = unit(0.2, "lines")) +
#     
#     coord_cartesian(
#       xlim = c(-11, 13),   # slightly wider right margin
#       ylim = c(-4, 10),
#       clip = "off"
#     ) +
#     
#     labs(
#       title = "Conceptual Structure of Population Genetics"
#     ) +
#     
#     theme_void(base_size = base_size) +
#     theme(
#       plot.title = element_text(
#         face = "bold",
#         size = base_size * 1.3,
#         hjust = 0.5
#       ),
#       plot.margin = margin(40, 60, 40, 60)
#     )
# }

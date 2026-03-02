## The Gene Pool Concept

library(ggplot2)

# -------------------------------------------------
# Parameters
# -------------------------------------------------
set.seed(1)

N <- 25
p <- 0.7
n_alleles <- 2 * N

# -------------------------------------------------
# Generate allele pool
# -------------------------------------------------
alleles <- sample(c(
  rep("A", round(p * n_alleles)),
  rep("a", n_alleles - round(p * n_alleles))
))

# Random circular layout
theta  <- runif(n_alleles, 0, 2*pi)
radius <- sqrt(runif(n_alleles))

df <- data.frame(
  x = radius * cos(theta),
  y = radius * sin(theta),
  allele = alleles
)

# -------------------------------------------------
# Plot
# -------------------------------------------------
ggplot(df, aes(x, y, fill = allele)) +
  geom_point(shape = 21, size = 6,
             color = "black", stroke = 0.7) +
  scale_fill_manual(values = c(
    "A" = "#D55E00",   # warm orange-red
    "a" = "#0072B2"    # deep blue
  )) +
  coord_equal() +
  labs(
    title = "The Gene Pool",
    subtitle = paste0(
      "Diploid population: N = ", N,
      "   |   Total alleles = ", n_alleles,
      "   |   p = ", round(p, 2)
    ),
    fill = "Allele"
  ) +
  theme_void(base_size = 20) +
  theme(legend.position = "right")


library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(1)

# -------------------------------------------------
# Parameters
# -------------------------------------------------
N <- 25
p <- 0.7
n_alleles <- 2 * N

# -------------------------------------------------
# Soft colour palette
# -------------------------------------------------
col_allele <- c(
  "A" = "#C26D38",
  "a" = "#4C78A8"
)

col_genotype <- c(
  "AA" = "#C26D38",
  "Aa" = "#8C6BB1",
  "aa" = "#4C78A8"
)

# -------------------------------------------------
# Generate allele pool
# -------------------------------------------------
alleles <- sample(c(
  rep("A", round(p * n_alleles)),
  rep("a", n_alleles - round(p * n_alleles))
))

geno <- tibble(
  id = rep(1:N, each = 2),
  copy = rep(c(1, 2), N),
  allele = alleles
)

geno_labels <- geno %>%
  group_by(id) %>%
  summarise(genotype = paste(sort(allele), collapse = "")) %>%
  mutate(genotype = ifelse(genotype == "aA", "Aa", genotype))

geno <- left_join(geno, geno_labels, by = "id")

# -------------------------------------------------
# Layout grid (5 × 5)
# -------------------------------------------------
geno$row <- ceiling(geno$id / 5)
geno$col <- (geno$id - 1) %% 5 + 1

geno$x <- geno$col + ifelse(geno$copy == 1, -0.10, 0.10)
geno$y <- -geno$row

label_df <- geno_labels %>%
  mutate(
    row = ceiling(id / 5),
    col = (id - 1) %% 5 + 1,
    x = col,
    y = -row + 0.30   # slightly closer to alleles
  )

# -------------------------------------------------
# Panel 1: Individuals
# -------------------------------------------------
p_ind <- ggplot() +
  geom_point(data = geno,
             aes(x, y, fill = allele),
             shape = 21,
             size = 6,
             color = "black",
             stroke = 0.7) +
  geom_text(data = label_df,
            aes(x, y, label = genotype),
            size = 4.2,              # smaller
            fontface = "plain",      # not bold
            color = "grey20") +      # softer than pure black
  scale_fill_manual(values = col_allele) +
  coord_equal() +
  labs(title = "Diploid individuals",
       fill = "Allele") +
  theme_void(base_size = 20) +
  theme(legend.position = "right")

# -------------------------------------------------
# Panel 2: Genotype counts
# -------------------------------------------------
geno_counts <- geno_labels %>%
  count(genotype)

p_counts <- ggplot(geno_counts,
                   aes(genotype, n, fill = genotype)) +
  geom_col(width = 0.6,
           color = "black") +
  scale_fill_manual(values = col_genotype) +
  labs(title = "Genotype counts",
       x = NULL,
       y = "Count") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

# -------------------------------------------------
# Combine
# -------------------------------------------------
p_ind | p_counts +
  plot_annotation(
    title = "From Gene Pool to Genotype Counts",
    theme = theme(plot.title = element_text(size = 24, face = "bold"))
  )



## Genetic Drift Depends on Population Size

# --- Libraries (once per document ideally) ---
library(ggplot2)
library(dplyr)

# --- Global theme (define once in setup chunk ideally) ---
theme_slide <- theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.minor = element_blank()
  )

# --- Wright–Fisher simulation ---
wf_drift <- function(N, p0, T, replicates) {
  results <- data.frame()
  
  for (r in 1:replicates) {
    p <- numeric(T + 1)
    p[1] <- p0
    
    for (t in 1:T) {
      p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
    }
    
    results <- rbind(results,
                     data.frame(
                       generation = 0:T,
                       p = p,
                       replicate = r,
                       N = paste0("N = ", N)
                     ))
  }
  results
}

set.seed(1)

T <- 100
p0 <- 0.7
replicates <- 20

df <- bind_rows(
  wf_drift(20,  p0, T, replicates),
  wf_drift(200, p0, T, replicates)
)

# --- Plot ---
ggplot(df, aes(generation, p, group = replicate)) +
  geom_line(alpha = 0.35, linewidth = 0.8) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.6,
               color = "#D55E00") +
  geom_hline(yintercept = p0,
             linetype = "dashed",
             linewidth = 1) +
  facet_wrap(~N) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Generation",
    y = "Allele frequency (p)"
  ) +
  theme_slide


## Selection Changes Allele Frequencies

library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------------------------------
# Global theme (define once in setup ideally)
# -------------------------------------------------
theme_slide <- theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# -------------------------------------------------
# Colour system (consistent biology)
# -------------------------------------------------
col_genotype <- c(
  "AA" = "#D55E00",
  "Aa" = "#7B3294",
  "aa" = "#0072B2"
)

# -------------------------------------------------
# Selection model
# -------------------------------------------------
step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  (p^2*wAA + p*q*wAa) / wbar
}

simulate_selection <- function(p0, T, wAA, wAa, waa) {
  p <- numeric(T + 1)
  p[1] <- p0
  for (t in 1:T)
    p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
  data.frame(gen = 0:T, p = p)
}

geno_freq <- function(p) {
  q <- 1 - p
  data.frame(
    genotype = c("AA","Aa","aa"),
    freq = c(p^2, 2*p*q, q^2)
  )
}

# -------------------------------------------------
# Parameters
# -------------------------------------------------
p0 <- 0.7
T  <- 30

# Directional selection
dir <- simulate_selection(p0, T, wAA=1.0, wAa=0.9, waa=0.6)

# Heterozygote advantage
het <- simulate_selection(p0, T, wAA=0.8, wAa=1.0, waa=0.8)

# -------------------------------------------------
# Allele trajectories
# -------------------------------------------------
p_dir_traj <- ggplot(dir, aes(gen, p)) +
  geom_line(linewidth = 1.6, color = "#D55E00") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Directional selection",
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_slide

p_het_traj <- ggplot(het, aes(gen, p)) +
  geom_line(linewidth = 1.6, color = "#0072B2") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Heterozygote advantage",
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_slide

# -------------------------------------------------
# Genotype frequencies: Gen 0 vs Gen T
# -------------------------------------------------
geno_dir <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(dir$p[T+1]) %>% mutate(generation = paste0("Gen ", T))
)

geno_het <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(het$p[T+1]) %>% mutate(generation = paste0("Gen ", T))
)

p_dir_geno <- ggplot(geno_dir,
                     aes(genotype, freq, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  scale_fill_manual(values = col_genotype) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL,
       y = "Genotype frequency",
       fill = "") +
  theme_slide

p_het_geno <- ggplot(geno_het,
                     aes(genotype, freq, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  scale_fill_manual(values = col_genotype) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL,
       y = "Genotype frequency",
       fill = "") +
  theme_slide

# -------------------------------------------------
# Combine layout
# -------------------------------------------------
(p_dir_traj | p_het_traj) /
  (p_dir_geno | p_het_geno)



## Hardy–Weinberg: Observed vs Expected


library(ggplot2)
library(dplyr)
library(tidyr)

# -------------------------------------------------
# Observed genotype counts
# -------------------------------------------------
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count = c(40, 20, 40)
)

N <- sum(obs$count)

# -------------------------------------------------
# Estimate allele frequency
# -------------------------------------------------
p_hat <- with(obs,
              (2*count[genotype=="AA"] +
                 count[genotype=="Aa"]) / (2*N))
q_hat <- 1 - p_hat

# -------------------------------------------------
# Combine observed + expected (HWE)
# -------------------------------------------------
df <- obs %>%
  mutate(expected = c(p_hat^2, 2*p_hat*q_hat, q_hat^2) * N) %>%
  pivot_longer(c(count, expected),
               names_to = "type",
               values_to = "n") %>%
  mutate(type = recode(type,
                       count = "Observed",
                       expected = "Expected (HWE)"))

# -------------------------------------------------
# Plot
# -------------------------------------------------
ggplot(df, aes(genotype, n, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  scale_fill_manual(values = c(
    "Observed" = "grey60",
    "Expected (HWE)" = "grey85"
  )) +
  labs(
    title = "Hardy–Weinberg: Observed vs Expected",
    subtitle = paste0("N = ", N,
                      "   |   p̂ = ", round(p_hat, 3),
                      "   |   q̂ = ", round(q_hat, 3)),
    x = NULL,
    y = "Genotype count",
    fill = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")


## Inbreeding Redistributes Genotypes (Not Alleles)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Example observed counts (heterozygote deficit)
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count    = c(40, 20, 40)
)

N <- sum(obs$count)

# Allele frequency
nA <- 2*obs$count[obs$genotype=="AA"] + obs$count[obs$genotype=="Aa"]
p_hat <- as.numeric(nA) / (2*N)
q_hat <- 1 - p_hat

# Heterozygosity and F estimate
H_obs <- obs$count[obs$genotype=="Aa"] / N
H_exp <- 2*p_hat*q_hat
F_hat <- 1 - H_obs / H_exp

# ---- Genotype frequencies ----
geno_hwe <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2, 2*p_hat*q_hat, q_hat^2),
  model = "HWE"
)

geno_F <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2 + F_hat*p_hat*q_hat,
           2*p_hat*q_hat*(1 - F_hat),
           q_hat^2 + F_hat*p_hat*q_hat),
  model = "Inbreeding"
)

geno_long <- bind_rows(geno_hwe, geno_F)

# ---- Allele frequencies (identical) ----
allele_df <- tibble(
  allele = rep(c("A","a"), 2),
  freq = rep(c(p_hat, q_hat), 2),
  model = rep(c("HWE","Inbreeding"), each = 2)
)

# ---- Plot 1: Genotypes ----
p1 <- ggplot(geno_long,
             aes(genotype, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Genotype frequencies",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Plot 2: Alleles ----
p2 <- ggplot(allele_df,
             aes(allele, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Allele frequencies (unchanged)",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

# ---- Combine ----
p1 | p2 +
  plot_annotation(
    title = "Inbreeding redistributes genotypes without changing allele frequency",
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )

## Allele and Genotype Frequencies

library(ggplot2)
library(patchwork)

# Parameters
N <- 200
p <- 0.7
set.seed(1)

# ---------------------------
# Simulate Gene Pool
# ---------------------------
n_alleles <- 2 * N
n_A <- round(p * n_alleles)

alleles <- c(rep("A", n_A),
             rep("a", n_alleles - n_A))
alleles <- sample(alleles)

# ---------------------------
# Form Genotypes
# ---------------------------
genotypes <- paste0(alleles[seq(1, n_alleles, 2)],
                    alleles[seq(2, n_alleles, 2)])

genotypes <- ifelse(genotypes %in% c("Aa", "aA"), "Aa", genotypes)

df_allele  <- data.frame(allele = alleles)
df_geno    <- data.frame(genotype = genotypes)

# ---------------------------
# Allele Counts
# ---------------------------
allele_counts <- as.data.frame(table(df_allele$allele))
colnames(allele_counts) <- c("allele", "count")

p1 <- ggplot(allele_counts,
             aes(x = allele, y = count, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "A" = "#C26D38",
    "a" = "#4C78A8"
  )) +
  labs(title = "Allele Counts (Gene Pool)",
       subtitle = paste("Total alleles = 2N =", 2*N),
       x = "Allele", y = "Count") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Allele Frequencies
# ---------------------------
allele_freq <- as.data.frame(prop.table(table(df_allele$allele)))
colnames(allele_freq) <- c("allele", "frequency")

p2 <- ggplot(allele_freq,
             aes(x = allele, y = frequency, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "A" = "#C26D38",
    "a" = "#4C78A8"
  )) +
  ylim(0,1) +
  labs(title = "Allele Frequencies",
       subtitle = paste("p =", round(allele_freq$frequency[1],2),
                        "| q =", round(allele_freq$frequency[2],2)),
       x = "Allele", y = "Frequency") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Genotype Counts
# ---------------------------
geno_counts <- as.data.frame(table(df_geno$genotype))
colnames(geno_counts) <- c("genotype", "count")

p3 <- ggplot(geno_counts,
             aes(x = genotype, y = count, fill = genotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "AA" = "#C26D38",
    "Aa" = "#8C6BB1",
    "aa" = "#4C78A8"
  )) +
  labs(title = "Genotype Counts",
       subtitle = paste("Total individuals = N =", N),
       x = "Genotype", y = "Count") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Genotype Frequencies
# ---------------------------
geno_freq <- as.data.frame(prop.table(table(df_geno$genotype)))
colnames(geno_freq) <- c("genotype", "frequency")

p4 <- ggplot(geno_freq,
             aes(x = genotype, y = frequency, fill = genotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c(
    "AA" = "#C26D38",
    "Aa" = "#8C6BB1",
    "aa" = "#4C78A8"
  )) +
  ylim(0,1) +
  labs(title = "Genotype Frequencies",
       subtitle = expression(AA == p^2 ~ "," ~ Aa == 2*p*q ~ "," ~ aa == q^2),
       x = "Genotype", y = "Frequency") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Combine 2×2 Layout
# ---------------------------
(p1 | p2) /
  (p3 | p4)



## Genetic Drift Across Replicate Populations

set.seed(1)

replicates <- 20
N <- 50
generations <- 100
p0 <- 0.7

results <- data.frame()

for (r in 1:replicates) {
  p <- numeric(generations + 1)
  p[1] <- p0
  
  for (t in 1:generations) {
    p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
  }
  
  temp <- data.frame(
    generation = 0:generations,
    frequency = p,
    replicate = as.factor(r)
  )
  
  results <- rbind(results, temp)
}

ggplot(results, aes(generation, frequency, group = replicate)) +
  geom_line(alpha = 0.6) +
  ylim(0,1) +
  labs(
    title = "Genetic Drift Across Replicate Populations",
    subtitle = paste("Wright–Fisher Model | N =", N),
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")



## Genetic Drift Depends on Population Size

library(ggplot2)

set.seed(1)

# Parameters
generations <- 100
replicates  <- 20
p0 <- 0.7

simulate_drift <- function(N, generations, replicates, p0) {
  results <- data.frame()
  
  for (r in 1:replicates) {
    p <- numeric(generations + 1)
    p[1] <- p0
    
    for (t in 1:generations) {
      p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
    }
    
    temp <- data.frame(
      generation = 0:generations,
      frequency  = p,
      replicate  = as.factor(r),
      N = paste0("N = ", N)
    )
    
    results <- rbind(results, temp)
  }
  
  return(results)
}

# Simulate both population sizes
res_small <- simulate_drift(20, generations, replicates, p0)
res_large <- simulate_drift(200, generations, replicates, p0)

results <- rbind(res_small, res_large)

# Plot
ggplot(results,
       aes(generation, frequency, group = replicate)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = p0,
             linetype = "dashed",
             color = "red") +
  ylim(0,1) +
  facet_wrap(~N) +
  labs(
    title = "Genetic Drift Depends on Population Size",
    subtitle = "Wright–Fisher Model | Same starting frequency",
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")



















##########################################################
# Selection 1-generation change in genotype frequencies
##########################################################

theme_slide <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

library(ggplot2)
library(patchwork)

p <- 0.7
q <- 1 - p

wAA <- 1.0
wAa <- 0.9
waa <- 0.6

geno_before <- data.frame(
  genotype = c("AA","Aa","aa"),
  freq = c(p^2, 2*p*q, q^2)
)

wbar <- sum(geno_before$freq * c(wAA, wAa, waa))

p_prime <- (p^2 * wAA + p*q * wAa) / wbar

alleles_long <- data.frame(
  allele = rep(c("A","a"), 2),
  stage  = rep(c("Before","After"), each = 2),
  freq   = c(p, q, p_prime, 1 - p_prime)
)


##########################################################
# Selection 1 vs 10 generation changes in genotype/allele frequencies
##########################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Parameters
p0 <- 0.7
wAA <- 1.0
wAa <- 0.9
waa <- 0.6
T <- 10

# Selection recursion
step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  (p^2*wAA + p*q*wAa) / wbar
}

# Iterate allele frequency
p <- numeric(T + 1)
p[1] <- p0
for (t in 1:T) {
  p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
}

p10 <- p[T + 1]

# Function to get genotype frequencies
geno_freq <- function(p) {
  q <- 1 - p
  c(AA = p^2, Aa = 2*p*q, aa = q^2)
}

df_geno <- rbind(
  data.frame(generation = "Gen 0",
             genotype = names(geno_freq(p0)),
             freq = geno_freq(p0)),
  data.frame(generation = "Gen 10",
             genotype = names(geno_freq(p10)),
             freq = geno_freq(p10))
)

df_allele <- data.frame(
  generation = rep(c("Gen 0","Gen 10"), each = 2),
  allele = rep(c("A","a"), 2),
  freq = c(p0, 1 - p0, p10, 1 - p10)
)

# Genotype plot
p1 <- ggplot(df_geno, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  ylim(0,1) +
  labs(title = "Genotype Frequencies",
       y = "Frequency", x = NULL, fill = "") +
  theme_minimal(base_size = 14)

# Allele plot
p2 <- ggplot(df_allele, aes(allele, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  ylim(0,1) +
  labs(title = "Allele Frequencies",
       y = "Frequency", x = NULL, fill = "") +
  theme_minimal(base_size = 14)

p1 | p2


##########################################################
# Directional selection vs heterozygote advantage
##########################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Helper functions ----
step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  (p^2*wAA + p*q*wAa) / wbar
}

simulate_selection <- function(p0, T, wAA, wAa, waa) {
  p <- numeric(T + 1)
  p[1] <- p0
  for (t in 1:T) p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
  tibble(gen = 0:T, p = p)
}

geno_freq <- function(p) {
  q <- 1 - p
  tibble(genotype = c("AA","Aa","aa"),
         freq = c(p^2, 2*p*q, q^2))
}

# ---- Parameters ----
p0 <- 0.7
T  <- 30

dir <- list(name = "Directional selection",
            wAA = 1.0, wAa = 0.9, waa = 0.6)

het <- list(name = "Heterozygote advantage",
            wAA = 0.8, wAa = 1.0, waa = 0.8)

# ---- Simulate ----
df_dir <- simulate_selection(p0, T, dir$wAA, dir$wAa, dir$waa)
df_het <- simulate_selection(p0, T, het$wAA, het$wAa, het$waa)

geno_dir <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(df_dir$p[T + 1]) %>% mutate(generation = paste0("Gen ", T))
)

geno_het <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(df_het$p[T + 1]) %>% mutate(generation = paste0("Gen ", T))
)

# ---- Plots ----
p_traj_dir <- ggplot(df_dir, aes(gen, p)) +
  geom_line(linewidth = 1.3) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = dir$name,
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_traj_het <- ggplot(df_het, aes(gen, p)) +
  geom_line(linewidth = 1.3) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = het$name,
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_geno_dir <- ggplot(geno_dir, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Genotype frequency", fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_geno_het <- ggplot(geno_het, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Genotype frequency", fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Combine (no subtitle) ----
(p_traj_dir | p_traj_het) /
  (p_geno_dir | p_geno_het) +
  plot_annotation(
    title = "Directional selection vs heterozygote advantage",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )

# -------------------------------------------------
# Panel 1: Genotype frequencies (neutral soft grey)
# -------------------------------------------------
p1 <- ggplot(geno_before, aes(genotype, freq)) +
  geom_col(width = 0.6, fill = "grey70", color="black") +
  ylim(0,1) +
  labs(title = "Genotypes") +
  theme_slide

# -------------------------------------------------
# Panel 2: Fitness (soft blue)
# -------------------------------------------------
fitness_df <- data.frame(
  genotype = c("AA","Aa","aa"),
  fitness  = c(wAA, wAa, waa)
)

p2 <- ggplot(fitness_df, aes(genotype, fitness)) +
  geom_col(width = 0.6, fill = "#4C78A8", color="black") +
  labs(title = "Fitness (Selection)") +
  theme_slide

# -------------------------------------------------
# Panel 3: Allele frequency change
# -------------------------------------------------
p3 <- ggplot(alleles_long, aes(allele, freq, fill = stage)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color="black") +
  ylim(0,1) +
  scale_fill_manual(values = c(
    "Before" = "grey75",
    "After"  = "#C26D38"
  )) +
  labs(title = "Allele Frequency Change", fill = "") +
  theme_slide

(p1 | p2) / p3


# HWE

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Example data (heterozygote deficit) ----
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count = c(40, 20, 40)
)

N <- sum(obs$count)

# Estimate allele frequency
nA <- 2*obs$count[obs$genotype=="AA"] + obs$count[obs$genotype=="Aa"]
p_hat <- as.numeric(nA) / (2*N)
q_hat <- 1 - p_hat

# Expected under HWE
exp <- tibble(
  genotype = c("AA","Aa","aa"),
  expected = c(p_hat^2, 2*p_hat*q_hat, q_hat^2) * N
)

df <- obs %>%
  left_join(exp, by="genotype")

# ---- Panel 1: Observed vs Expected ----
df_long <- df %>%
  pivot_longer(cols = c(count, expected),
               names_to = "type",
               values_to = "n") %>%
  mutate(type = recode(type,
                       count = "Observed",
                       expected = "Expected (HWE)"))

p1 <- ggplot(df_long,
             aes(genotype, n, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  labs(title = "Observed vs Expected (HWE)",
       x = NULL,
       y = "Genotype count",
       fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Panel 2: Chi-square contributions ----
chi_df <- df %>%
  mutate(contrib = (count - expected)^2 / expected)

X2 <- sum(chi_df$contrib)
pval <- pchisq(X2, df = 1, lower.tail = FALSE)

p2 <- ggplot(chi_df,
             aes(genotype, contrib)) +
  geom_col(width = 0.6,
           fill = "grey70",
           color = "black") +
  labs(title = paste0("χ² contributions  (χ² = ",
                      round(X2,2),
                      ",  p = ",
                      signif(pval,3),
                      ")"),
       x = NULL,
       y = "Contribution to χ²") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Combine ----
p1 | p2


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Example observed counts (heterozygote deficit)
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count    = c(40, 20, 40)
)

N <- sum(obs$count)

# Allele frequency
nA <- 2*obs$count[obs$genotype=="AA"] + obs$count[obs$genotype=="Aa"]
p_hat <- as.numeric(nA) / (2*N)
q_hat <- 1 - p_hat

# Heterozygosity and F estimate
H_obs <- obs$count[obs$genotype=="Aa"] / N
H_exp <- 2*p_hat*q_hat
F_hat <- 1 - H_obs / H_exp

# ---- Genotype frequencies ----
geno_hwe <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2, 2*p_hat*q_hat, q_hat^2),
  model = "HWE"
)

geno_F <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2 + F_hat*p_hat*q_hat,
           2*p_hat*q_hat*(1 - F_hat),
           q_hat^2 + F_hat*p_hat*q_hat),
  model = "Inbreeding"
)

geno_long <- bind_rows(geno_hwe, geno_F)

# ---- Allele frequencies (identical) ----
allele_df <- tibble(
  allele = rep(c("A","a"), 2),
  freq = rep(c(p_hat, q_hat), 2),
  model = rep(c("HWE","Inbreeding"), each = 2)
)

# ---- Plot 1: Genotypes ----
p1 <- ggplot(geno_long,
             aes(genotype, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Genotype frequencies",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Plot 2: Alleles ----
p2 <- ggplot(allele_df,
             aes(allele, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Allele frequencies (unchanged)",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

# ---- Combine ----
p1 | p2 +
  plot_annotation(
    title = "Inbreeding redistributes genotypes without changing allele frequency",
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )



# Admixture

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Parameters ----
m  <- 0.5   # fraction from Pop1 (so Pop2 is 1-m)
p1 <- 0.2   # allele A frequency in Pop1
p2 <- 0.8   # allele A frequency in Pop2

# ---- Helper: HWE genotype freqs ----
geno_hwe <- function(p){
  q <- 1 - p
  c(AA = p^2, Aa = 2*p*q, aa = q^2)
}

# ---- Pooled allele frequency ----
pbar <- m*p1 + (1-m)*p2

# ---- Within-pop HWE (assume each pop is in HWE) ----
g1 <- geno_hwe(p1)
g2 <- geno_hwe(p2)

# ---- Pooled (admixed) genotype frequencies: mixture of the two populations ----
g_pool_obs <- m*g1 + (1-m)*g2

# ---- HWE expectation if you (incorrectly) treat pooled sample as one random-mating pop ----
g_pool_hwe <- geno_hwe(pbar)

# ---- Data for plotting ----
allele_df <- tibble(
  pop = c("Pop 1", "Pop 2", "Pooled"),
  pA  = c(p1, p2, pbar),
  pa  = 1 - pA
) %>%
  pivot_longer(cols = c(pA, pa), names_to = "allele", values_to = "freq") %>%
  mutate(allele = recode(allele, pA = "A", pa = "a"))

geno_df <- tibble(
  genotype = c("AA","Aa","aa"),
  Observed_pooled = as.numeric(g_pool_obs),
  Expected_HWE    = as.numeric(g_pool_hwe)
) %>%
  pivot_longer(cols = c(Observed_pooled, Expected_HWE),
               names_to = "type", values_to = "freq") %>%
  mutate(type = recode(type,
                       Observed_pooled = "Observed",
                       Expected_HWE    = "HWE Expectation"))

# ---- Panel 1: Allele frequencies ----
p1_plot <- ggplot(allele_df, aes(pop, freq, fill = allele)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "#d73027", "a" = "#4575b4")) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Allele frequencies",
       x = NULL, y = "Frequency") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Panel 2: Genotype frequencies (Wahlund effect) ----
p2_plot <- ggplot(geno_df, aes(genotype, freq, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Pooled genotypes deviate from HWE (Wahlund effect)",
       subtitle = paste0(
         "50% from each population   |   Overall allele frequency = ",
         round(pbar, 2)),
       x = NULL, y = "Frequency", fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

p1_plot | p2_plot


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(1)

# -------------------------
# Parameters you can tweak
# -------------------------
n_ind <- 80              # individuals
K <- 2

p1 <- 0.2                # allele A frequency in ancestry 1
p2 <- 0.8                # allele A frequency in ancestry 2

# admixture proportions per individual (q = ancestry-1 proportion)
# smaller "conc" => more extreme bars; larger => more similar individuals
mean_q <- 0.5
conc   <- 6
alpha  <- mean_q * conc
beta   <- (1 - mean_q) * conc

# -------------------------
# Simulate ADMIXTURE Q-matrix (K=2)
# -------------------------
q1 <- rbeta(n_ind, alpha, beta)     # ancestry 1 proportion
q2 <- 1 - q1

Q_long <- tibble(
  id = 1:n_ind,
  anc1 = q1,
  anc2 = q2
) %>%
  pivot_longer(cols = c(anc1, anc2), names_to = "ancestry", values_to = "q") %>%
  mutate(ancestry = recode(ancestry, anc1 = "Ancestry 1", anc2 = "Ancestry 2"))

# order individuals by ancestry-1 proportion (common in STRUCTURE plots)
id_order <- tibble(id = 1:n_ind, q1 = q1) %>% arrange(desc(q1)) %>% pull(id)
Q_long$id <- factor(Q_long$id, levels = id_order)

# -------------------------
# Simulate genotypes under local ancestry
# (each allele chooses ancestry ~ Bernoulli(q1); then allele ~ Bernoulli(p_ancestry))
# -------------------------
draw_allele <- function(q1_i) {
  anc <- rbinom(1, 1, q1_i)              # 1 => ancestry 1, 0 => ancestry 2
  pA  <- ifelse(anc == 1, p1, p2)
  rbinom(1, 1, pA)                       # 1 => A, 0 => a
}

geno <- sapply(q1, function(qi) {
  a1 <- draw_allele(qi)
  a2 <- draw_allele(qi)
  # genotype coded as number of A alleles: 0,1,2
  a1 + a2
})

# Observed pooled genotype counts
obs_counts <- tibble(
  genotype = c("AA","Aa","aa"),
  count = c(sum(geno == 2), sum(geno == 1), sum(geno == 0))
)

N <- sum(obs_counts$count)

# Pooled allele frequency estimate from pooled genotypes
p_hat <- (2*obs_counts$count[obs_counts$genotype=="AA"] +
            obs_counts$count[obs_counts$genotype=="Aa"]) / (2*N)
q_hat <- 1 - p_hat

# Expected under HWE using pooled p_hat
exp_counts <- tibble(
  genotype = c("AA","Aa","aa"),
  expected = c(p_hat^2, 2*p_hat*q_hat, q_hat^2) * N
)

df_hwe <- obs_counts %>%
  left_join(exp_counts, by="genotype") %>%
  pivot_longer(cols = c(count, expected), names_to = "type", values_to = "n") %>%
  mutate(type = recode(type, count = "Observed (pooled)", expected = "Expected (HWE)"))

# Panel 1
p_admix <- ggplot(Q_long, aes(x = id, y = q, fill = ancestry)) +
  geom_col(width = 1) +
  labs(
    title = "Admixture proportions",
    x = NULL, y = "Ancestry proportion", fill = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

# Panel 2
p_hwe <- ggplot(df_hwe, aes(genotype, n, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  labs(
    title = "Pooled sample deviates from HWE",
    subtitle = paste0("Allele frequency = ", round(p_hat, 3)),
    x = NULL, y = "Genotype count", fill = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

p_admix | p_hwe +
  plot_annotation(
    title = "Admixture can generate Hardy–Weinberg deviation (Wahlund effect)",
    theme = theme(plot.title = element_text(size = 22, face = "bold"))
  )


# Selective sweep

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
# STATISTICS (same as Shiny app)
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
# EXACT SAME PLOTS AS SHINY (extracted)
# =========================================================
make_sweep_plots <- function(
    mode = c("Textbook Hard Sweep", "Realistic Hard Sweep"),
    n_hap = 100,
    L = 201,
    sweep_freq = 1.0,
    width = 30,          # core_width OR rec_scale
    sort_haps = TRUE,
    seed = NULL
) {
  mode <- match.arg(mode)
  
  L <- as.integer(L)
  if (L %% 2 == 0) L <- L + 1
  
  sweep_data <- if (mode == "Textbook Hard Sweep") {
    generate_textbook_sweep(
      n_hap = n_hap, L = L, sweep_freq = sweep_freq,
      core_width = width, seed = seed
    )
  } else {
    generate_realistic_sweep(
      n_hap = n_hap, L = L, sweep_freq = sweep_freq,
      rec_scale = width, seed = seed
    )
  }
  
  H <- sweep_data$H
  center <- sweep_data$center
  
  # Optional sorting to make sweep block visually clean (same as Shiny)
  if (sort_haps) {
    H <- H[order(H[, center], decreasing = TRUE), , drop = FALSE]
  }
  
  # ---- Haplotypes plot (same as Shiny) ----
  df_hap <- data.frame(
    Hap = rep(seq_len(nrow(H)), each = ncol(H)),
    Pos = rep(seq_len(ncol(H)), times = nrow(H)),
    Allele = as.vector(H)
  )
  df_hap$Allele <- as.vector(t(H))  # same final assignment as in Shiny app
  
  hapPlot <- ggplot(df_hap, aes(x = Pos, y = Hap, fill = factor(Allele))) +
    geom_tile() +
    geom_vline(xintercept = center, linetype = "dashed") +
    theme_minimal() +
    labs(fill = "Allele",
         title = "Haplotype structure around selected site")
  
  # ---- Gene diversity plot (same as Shiny) ----
  dist <- abs(seq_len(ncol(H)) - center)
  gd <- gene_diversity(H)
  df_gd <- data.frame(Distance = dist, GD = gd)
  
  gdPlot <- ggplot(df_gd, aes(x = Distance, y = GD)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    ylim(0, 0.5) +
    labs(title = "Diversity trough around selected site",
         y = "2p(1-p)")
  
  # ---- LD plot (same as Shiny) ----
  r2 <- ld_r2_with_center(H, center)
  df_ld <- data.frame(Distance = dist, r2 = r2)
  
  ldPlot <- ggplot(df_ld, aes(x = Distance, y = r2)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    ylim(0, 1) +
    labs(title = "LD peak near selected site",
         y = expression(r^2))
  
  list(hapPlot = hapPlot, gdPlot = gdPlot, ldPlot = ldPlot)
}

# =========================================================
# EXAMPLE: run once and print plots
# =========================================================
plots <- make_sweep_plots(
  mode = "Realistic Hard Sweep",
  n_hap = 100,
  L = 201,
  sweep_freq = 0.9,
  width = 30,
  sort_haps = TRUE,
  seed = 1
)

plots$hapPlot
plots$gdPlot
plots$ldPlot



library(tidyverse)

# ------------------------------------------------------------
# Aligned sequences
# ------------------------------------------------------------

seqs <- c(
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GCCATCGCTCC--TACTTAGAGAG---TTAGTC",
  "GCCATCGCTCC----CTTAGAGAG---TTAGTC",
  "GCCATCGCTCC----CTTAGAGAG---TTAGTC",
  "GCCATCGCTCCGTTACGTAGAGAG---CTTAGTC"
)

df <- tibble(
  id = factor(1:7),
  seq = seqs
) %>%
  mutate(pos = map(seq, ~ seq_along(strsplit(.x, "")[[1]])),
         base = map(seq, ~ strsplit(.x, "")[[1]])) %>%
  unnest(c(pos, base))

# SNP detection
snp_pos <- df %>%
  filter(base != "-") %>%
  group_by(pos) %>%
  summarise(n = n_distinct(base), .groups="drop") %>%
  filter(n > 1) %>%
  pull(pos)

# Regions
indel_region  <- 12:16
micro_region  <- 21:28

# ------------------------------------------------------------
# Stronger, high-contrast base colors
# ------------------------------------------------------------

base_cols <- c(
  "A" = "#D55E00",  # strong orange-red
  "T" = "#009E73",  # strong green
  "C" = "#0072B2",  # strong blue
  "G" = "#CC79A7",  # magenta
  "-" = "grey50"
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ggplot(df, aes(pos, fct_rev(id))) +
  
  # subtle background grid
  geom_tile(fill = "grey96", color = NA, height = 0.85) +
  
  # region highlights
  annotate("rect",
           xmin = min(indel_region)-0.5,
           xmax = max(indel_region)+0.5,
           ymin = 0.5, ymax = 7.5,
           fill = "grey88", alpha = 0.5) +
  
  annotate("rect",
           xmin = min(micro_region)-0.5,
           xmax = max(micro_region)+0.5,
           ymin = 0.5, ymax = 7.5,
           fill = "grey92", alpha = 0.5) +
  
  # letters (larger + bold monospace)
  geom_text(aes(label = base, color = base),
            family = "mono",
            fontface = "bold",
            size = 5.5) +
  
  scale_color_manual(values = base_cols) +
  
  scale_x_continuous(
    breaks = seq(0, 40, 5),
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  
  labs(
    x = "Genomic position",
    y = NULL,
    title = "Sequence variation across homologous chromosomes"
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  
  # SNP markers (slightly thicker)
  geom_segment(
    data = tibble(pos = snp_pos),
    aes(x = pos, xend = pos, y = 7.6, yend = 7.9),
    inherit.aes = FALSE,
    linewidth = 0.8
  )


library(tidyverse)

# ------------------------------------------------------------
# Aligned sequences
# ------------------------------------------------------------

seqs <- c(
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGCATCGCGCCGTTACGTAGAGAGAGGTGAATC",
  "GGTATCGCTCC--TACTTAGAGAG---TTAGTC",  # C->T transition near start
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",  # same indel block
  "GGTATCGCTCC----CTTAGAGAG---TTAGTC",
  "GGTATCGCTCC--TACTTAGAGAG---CTTAGT"  # add "--" to match length/structure
)

nchar(seqs)
stopifnot(length(unique(nchar(seqs))) == 1)
df <- tibble(
  id = factor(1:7),
  seq = seqs
) %>%
  mutate(pos = map(seq, ~ seq_along(strsplit(.x, "")[[1]])),
         base = map(seq, ~ strsplit(.x, "")[[1]])) %>%
  unnest(c(pos, base))

# ------------------------------------------------------------
# Detect segregating (polymorphic) positions
# ------------------------------------------------------------

seg_sites <- df %>%
  group_by(pos) %>%
  summarise(n = n_distinct(base), .groups="drop") %>%
  filter(n > 1) %>%
  pull(pos)

df <- df %>%
  mutate(is_seg = pos %in% seg_sites)

# ------------------------------------------------------------
# Publication-friendly palette
# ------------------------------------------------------------

base_cols <- c(
  "A" = "#D55E00",
  "T" = "#009E73",
  "C" = "#0072B2",
  "G" = "#CC79A7",
  "-" = "grey60"
)

n_pos <- max(df$pos)

ggplot(df, aes(pos, fct_rev(id))) +
  
  geom_text(aes(label = base,
                color = base,
                fontface = ifelse(is_seg, "bold", "plain")),
            family = "mono",
            size = 4.5) +
  
  scale_color_manual(values = base_cols) +
  
  scale_x_continuous(
    breaks = seq(0, n_pos, 5),
    limits = c(0.5, n_pos + 0.5),
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    limits = rev(levels(df$id)),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_fixed(ratio = 0.6, clip = "off") +
  
  labs(
    x = "Genomic position",
    y = NULL,
    title = "Sequence variation across homologous chromosomes"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(12, 12, 12, 12)
  )

n_id <- length(unique(df$id))

ggplot(df, aes(pos, fct_rev(id))) +
  
  # Microsatellite outline box
  annotate("rect",
           xmin = 12 - 0.5,
           xmax = 15 + 0.5,
           ymin = 0.5,
           ymax = n_id + 0.5,
           fill = NA,
           color = "black",
           linewidth = 0.7) +
  
  # Microsatellite outline box
  annotate("rect",
           xmin = 25 - 0.5,
           xmax = 27 + 0.5,
           ymin = 0.5,
           ymax = n_id + 0.5,
           fill = NA,
           color = "black",
           linewidth = 0.7) +
  
  
  geom_text(aes(label = base,
                color = base,
                fontface = ifelse(is_seg, "bold", "plain")),
            family = "mono",
            size = 4.5) +
  
  scale_color_manual(values = base_cols) +
  
  scale_x_continuous(
    breaks = seq(0, n_pos, 5),
    limits = c(0.5, n_pos + 0.5),
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    limits = rev(levels(df$id)),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_fixed(ratio = 0.6, clip = "off") +
  
  labs(
    x = "Genomic position",
    y = NULL,
    title = "Sequence variation across homologous chromosomes"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(12, 12, 12, 12)
  )








##########################################################
# Gene Pool
##########################################################

library(ggplot2)

# Parameters
N <- 25
p <- 0.7
set.seed(1)

n_alleles <- 2 * N
n_A <- round(p * n_alleles)
alleles <- c(rep("A", n_A), rep("a", n_alleles - n_A))
alleles <- sample(alleles)

# Random circular layout
theta <- runif(n_alleles, 0, 2*pi)
radius <- sqrt(runif(n_alleles, 0, 1))

df <- data.frame(
  x = radius * cos(theta),
  y = radius * sin(theta),
  allele = alleles
)

ggplot(df, aes(x, y, fill = allele)) +
  geom_point(shape = 21, size = 6, color = "black", stroke = 0.7) +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  coord_equal() +
  labs(
    title = "The Gene Pool",
    subtitle = paste("Diploid population: N =", N,
                     "| 2N =", 2*N,
                     "| p =", round(n_A/(2*N),2)),
    fill = "Allele"
  ) +
  theme_void(base_size = 18) +
  theme(legend.position = "right")


##########################################################
# Allele and Genotype Frequencies
##########################################################


library(ggplot2)
library(patchwork)

# Parameters
N <- 200
p <- 0.7
set.seed(1)

# ---------------------------
# Simulate Gene Pool
# ---------------------------
n_alleles <- 2 * N
n_A <- round(p * n_alleles)

alleles <- c(rep("A", n_A),
             rep("a", n_alleles - n_A))
alleles <- sample(alleles)

# ---------------------------
# Form Genotypes
# ---------------------------
genotypes <- paste0(alleles[seq(1, n_alleles, 2)],
                    alleles[seq(2, n_alleles, 2)])

# Standardize heterozygotes
genotypes <- ifelse(genotypes %in% c("Aa", "aA"), "Aa", genotypes)

# Data frames
df_allele  <- data.frame(allele = alleles)
df_geno    <- data.frame(genotype = genotypes)

# ---------------------------
# Allele Counts
# ---------------------------
allele_counts <- as.data.frame(table(df_allele$allele))
colnames(allele_counts) <- c("allele", "count")

p1 <- ggplot(allele_counts,
             aes(x = allele, y = count, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  labs(title = "Allele Counts (Gene Pool)",
       subtitle = paste("Total alleles = 2N =", 2*N),
       x = "Allele", y = "Count") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Allele Frequencies
# ---------------------------
allele_freq <- as.data.frame(prop.table(table(df_allele$allele)))
colnames(allele_freq) <- c("allele", "frequency")

p2 <- ggplot(allele_freq,
             aes(x = allele, y = frequency, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  ylim(0,1) +
  labs(title = "Allele Frequencies",
       subtitle = paste("p =", round(allele_freq$frequency[1],2),
                        "| q =", round(allele_freq$frequency[2],2)),
       x = "Allele", y = "Frequency") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Genotype Counts
# ---------------------------
geno_counts <- as.data.frame(table(df_geno$genotype))
colnames(geno_counts) <- c("genotype", "count")

p3 <- ggplot(geno_counts,
             aes(x = genotype, y = count, fill = genotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("AA" = "darkred",
                               "Aa" = "purple",
                               "aa" = "darkblue")) +
  labs(title = "Genotype Counts",
       subtitle = paste("Total individuals = N =", N),
       x = "Genotype", y = "Count") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Genotype Frequencies
# ---------------------------
geno_freq <- as.data.frame(prop.table(table(df_geno$genotype)))
colnames(geno_freq) <- c("genotype", "frequency")

p4 <- ggplot(geno_freq,
             aes(x = genotype, y = frequency, fill = genotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("AA" = "darkred",
                               "Aa" = "purple",
                               "aa" = "darkblue")) +
  ylim(0,1) +
  labs(title = "Genotype Frequencies",
       subtitle = expression(AA == p^2 ~ "," ~ Aa == 2*p*q ~ "," ~ aa == q^2),
       x = "Genotype", y = "Frequency") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# ---------------------------
# Combine 2×2 Layout
# ---------------------------
(p1 | p2) /
  (p3 | p4)

##########################################################
# Drift
##########################################################

library(ggplot2)

set.seed(1)

# Parameters
N <- 50              # population size
generations <- 100
p0 <- 0.7            # initial allele frequency

p <- numeric(generations + 1)
p[1] <- p0

for (t in 1:generations) {
  p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
}

df <- data.frame(
  generation = 0:generations,
  frequency = p
)

ggplot(df, aes(generation, frequency)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = p0, linetype = "dashed", color = "red") +
  ylim(0,1) +
  labs(
    title = "Genetic Drift (Wright–Fisher Model)",
    subtitle = paste("Population size N =", N,
                     "| Initial p =", p0),
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal(base_size = 16)


##########################################################
# Drift multiple replicates
##########################################################

set.seed(1)

replicates <- 20
N <- 50
generations <- 100
p0 <- 0.7

results <- data.frame()

for (r in 1:replicates) {
  p <- numeric(generations + 1)
  p[1] <- p0
  
  for (t in 1:generations) {
    p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
  }
  
  temp <- data.frame(
    generation = 0:generations,
    frequency = p,
    replicate = as.factor(r)
  )
  
  results <- rbind(results, temp)
}

ggplot(results, aes(generation, frequency, group = replicate)) +
  geom_line(alpha = 0.6) +
  ylim(0,1) +
  labs(
    title = "Genetic Drift Across Replicate Populations",
    subtitle = paste("Wright–Fisher Model | N =", N),
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


##########################################################
# Drift with N=20 or N=200
##########################################################

library(ggplot2)

set.seed(1)

# Parameters
generations <- 100
replicates  <- 20
p0 <- 0.7

simulate_drift <- function(N, generations, replicates, p0) {
  results <- data.frame()
  
  for (r in 1:replicates) {
    p <- numeric(generations + 1)
    p[1] <- p0
    
    for (t in 1:generations) {
      p[t + 1] <- rbinom(1, 2*N, p[t]) / (2*N)
    }
    
    temp <- data.frame(
      generation = 0:generations,
      frequency  = p,
      replicate  = as.factor(r),
      N = paste0("N = ", N)
    )
    
    results <- rbind(results, temp)
  }
  
  return(results)
}

# Simulate both population sizes
res_small <- simulate_drift(20, generations, replicates, p0)
res_large <- simulate_drift(200, generations, replicates, p0)

results <- rbind(res_small, res_large)

# Plot
ggplot(results,
       aes(generation, frequency, group = replicate)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = p0,
             linetype = "dashed",
             color = "red") +
  ylim(0,1) +
  facet_wrap(~N) +
  labs(
    title = "Genetic Drift Depends on Population Size",
    subtitle = "Wright–Fisher Model | Same starting frequency",
    x = "Generation",
    y = "Allele Frequency (p)"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

##########################################################
# Selection 1-generation change in genotype frequencies
##########################################################

theme_slide <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    panel.grid.minor = element_blank()
  )
library(ggplot2)
library(patchwork)

p <- 0.7
q <- 1 - p

wAA <- 1.0
wAa <- 0.9
waa <- 0.6

geno_before <- data.frame(
  genotype = c("AA","Aa","aa"),
  freq = c(p^2, 2*p*q, q^2)
)

wbar <- sum(geno_before$freq * c(wAA, wAa, waa))

p_prime <- (p^2 * wAA + p*q * wAa) / wbar

alleles_long <- data.frame(
  allele = rep(c("A","a"), 2),
  stage  = rep(c("Before","After"), each = 2),
  freq   = c(p, q, p_prime, 1 - p_prime)
)

# Panel 1
p1 <- ggplot(geno_before, aes(genotype, freq)) +
  geom_col(width = 0.6, fill = "grey60", color="black") +
  ylim(0,1) +
  labs(title = "Genotypes") +
  theme_slide

# Panel 2
fitness_df <- data.frame(
  genotype = c("AA","Aa","aa"),
  fitness  = c(wAA, wAa, waa)
)

p2 <- ggplot(fitness_df, aes(genotype, fitness)) +
  geom_col(width = 0.6, fill = "steelblue", color="black") +
  labs(title = "Fitness (Selection)") +
  theme_slide

# Panel 3
p3 <- ggplot(alleles_long, aes(allele, freq, fill = stage)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color="black") +
  ylim(0,1) +
  scale_fill_manual(values = c("grey70","firebrick")) +
  labs(title = "Allele Frequency Change", fill = "") +
  theme_slide

(p1 | p2) / p3


##########################################################
# Selection 1 vs 10 generation changes in genotype/allele frequencies
##########################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Parameters
p0 <- 0.7
wAA <- 1.0
wAa <- 0.9
waa <- 0.6
T <- 10

# Selection recursion
step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  (p^2*wAA + p*q*wAa) / wbar
}

# Iterate allele frequency
p <- numeric(T + 1)
p[1] <- p0
for (t in 1:T) {
  p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
}

p10 <- p[T + 1]

# Function to get genotype frequencies
geno_freq <- function(p) {
  q <- 1 - p
  c(AA = p^2, Aa = 2*p*q, aa = q^2)
}

df_geno <- rbind(
  data.frame(generation = "Gen 0",
             genotype = names(geno_freq(p0)),
             freq = geno_freq(p0)),
  data.frame(generation = "Gen 10",
             genotype = names(geno_freq(p10)),
             freq = geno_freq(p10))
)

df_allele <- data.frame(
  generation = rep(c("Gen 0","Gen 10"), each = 2),
  allele = rep(c("A","a"), 2),
  freq = c(p0, 1 - p0, p10, 1 - p10)
)

# Genotype plot
p1 <- ggplot(df_geno, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  ylim(0,1) +
  labs(title = "Genotype Frequencies",
       y = "Frequency", x = NULL, fill = "") +
  theme_minimal(base_size = 14)

# Allele plot
p2 <- ggplot(df_allele, aes(allele, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  ylim(0,1) +
  labs(title = "Allele Frequencies",
       y = "Frequency", x = NULL, fill = "") +
  theme_minimal(base_size = 14)

p1 | p2


##########################################################
# Directional selection vs heterozygote advantage
##########################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Helper functions ----
step_selection <- function(p, wAA, wAa, waa) {
  q <- 1 - p
  wbar <- p^2*wAA + 2*p*q*wAa + q^2*waa
  (p^2*wAA + p*q*wAa) / wbar
}

simulate_selection <- function(p0, T, wAA, wAa, waa) {
  p <- numeric(T + 1)
  p[1] <- p0
  for (t in 1:T) p[t + 1] <- step_selection(p[t], wAA, wAa, waa)
  tibble(gen = 0:T, p = p)
}

geno_freq <- function(p) {
  q <- 1 - p
  tibble(genotype = c("AA","Aa","aa"),
         freq = c(p^2, 2*p*q, q^2))
}

# ---- Parameters ----
p0 <- 0.7
T  <- 30

dir <- list(name = "Directional selection",
            wAA = 1.0, wAa = 0.9, waa = 0.6)

het <- list(name = "Heterozygote advantage",
            wAA = 0.8, wAa = 1.0, waa = 0.8)

# ---- Simulate ----
df_dir <- simulate_selection(p0, T, dir$wAA, dir$wAa, dir$waa)
df_het <- simulate_selection(p0, T, het$wAA, het$wAa, het$waa)

geno_dir <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(df_dir$p[T + 1]) %>% mutate(generation = paste0("Gen ", T))
)

geno_het <- bind_rows(
  geno_freq(p0) %>% mutate(generation = "Gen 0"),
  geno_freq(df_het$p[T + 1]) %>% mutate(generation = paste0("Gen ", T))
)

# ---- Plots ----
p_traj_dir <- ggplot(df_dir, aes(gen, p)) +
  geom_line(linewidth = 1.3) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = dir$name,
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_traj_het <- ggplot(df_het, aes(gen, p)) +
  geom_line(linewidth = 1.3) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = het$name,
       x = "Generation",
       y = "Allele frequency (p)") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_geno_dir <- ggplot(geno_dir, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Genotype frequency", fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

p_geno_het <- ggplot(geno_het, aes(genotype, freq, fill = generation)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = NULL, y = "Genotype frequency", fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Combine (no subtitle) ----
(p_traj_dir | p_traj_het) /
  (p_geno_dir | p_geno_het) +
  plot_annotation(
    title = "Directional selection vs heterozygote advantage",
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )


# HWE

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Example data (heterozygote deficit) ----
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count = c(40, 20, 40)
)

N <- sum(obs$count)

# Estimate allele frequency
nA <- 2*obs$count[obs$genotype=="AA"] + obs$count[obs$genotype=="Aa"]
p_hat <- as.numeric(nA) / (2*N)
q_hat <- 1 - p_hat

# Expected under HWE
exp <- tibble(
  genotype = c("AA","Aa","aa"),
  expected = c(p_hat^2, 2*p_hat*q_hat, q_hat^2) * N
)

df <- obs %>%
  left_join(exp, by="genotype")

# ---- Panel 1: Observed vs Expected ----
df_long <- df %>%
  pivot_longer(cols = c(count, expected),
               names_to = "type",
               values_to = "n") %>%
  mutate(type = recode(type,
                       count = "Observed",
                       expected = "Expected (HWE)"))

p1 <- ggplot(df_long,
             aes(genotype, n, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  labs(title = "Observed vs Expected (HWE)",
       x = NULL,
       y = "Genotype count",
       fill = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Panel 2: Chi-square contributions ----
chi_df <- df %>%
  mutate(contrib = (count - expected)^2 / expected)

X2 <- sum(chi_df$contrib)
pval <- pchisq(X2, df = 1, lower.tail = FALSE)

p2 <- ggplot(chi_df,
             aes(genotype, contrib)) +
  geom_col(width = 0.6,
           fill = "grey70",
           color = "black") +
  labs(title = paste0("χ² contributions  (χ² = ",
                      round(X2,2),
                      ",  p = ",
                      signif(pval,3),
                      ")"),
       x = NULL,
       y = "Contribution to χ²") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

# ---- Combine ----
p1 | p2


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Example observed counts (heterozygote deficit)
obs <- tibble(
  genotype = c("AA","Aa","aa"),
  count    = c(40, 20, 40)
)

N <- sum(obs$count)

# Allele frequency
nA <- 2*obs$count[obs$genotype=="AA"] + obs$count[obs$genotype=="Aa"]
p_hat <- as.numeric(nA) / (2*N)
q_hat <- 1 - p_hat

# Heterozygosity and F estimate
H_obs <- obs$count[obs$genotype=="Aa"] / N
H_exp <- 2*p_hat*q_hat
F_hat <- 1 - H_obs / H_exp

# ---- Genotype frequencies ----
geno_hwe <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2, 2*p_hat*q_hat, q_hat^2),
  model = "HWE"
)

geno_F <- tibble(
  genotype = c("AA","Aa","aa"),
  freq = c(p_hat^2 + F_hat*p_hat*q_hat,
           2*p_hat*q_hat*(1 - F_hat),
           q_hat^2 + F_hat*p_hat*q_hat),
  model = "Inbreeding"
)

geno_long <- bind_rows(geno_hwe, geno_F)

# ---- Allele frequencies (identical) ----
allele_df <- tibble(
  allele = rep(c("A","a"), 2),
  freq = rep(c(p_hat, q_hat), 2),
  model = rep(c("HWE","Inbreeding"), each = 2)
)

# ---- Plot 1: Genotypes ----
p1 <- ggplot(geno_long,
             aes(genotype, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Genotype frequencies",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Plot 2: Alleles ----
p2 <- ggplot(allele_df,
             aes(allele, freq, fill = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Allele frequencies (unchanged)",
       x = NULL,
       y = "Frequency",
       fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

# ---- Combine ----
p1 | p2 +
  plot_annotation(
    title = "Inbreeding redistributes genotypes without changing allele frequency",
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )


# Admixture

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Parameters ----
m  <- 0.5   # fraction from Pop1 (so Pop2 is 1-m)
p1 <- 0.2   # allele A frequency in Pop1
p2 <- 0.8   # allele A frequency in Pop2

# ---- Helper: HWE genotype freqs ----
geno_hwe <- function(p){
  q <- 1 - p
  c(AA = p^2, Aa = 2*p*q, aa = q^2)
}

# ---- Pooled allele frequency ----
pbar <- m*p1 + (1-m)*p2

# ---- Within-pop HWE (assume each pop is in HWE) ----
g1 <- geno_hwe(p1)
g2 <- geno_hwe(p2)

# ---- Pooled (admixed) genotype frequencies: mixture of the two populations ----
g_pool_obs <- m*g1 + (1-m)*g2

# ---- HWE expectation if you (incorrectly) treat pooled sample as one random-mating pop ----
g_pool_hwe <- geno_hwe(pbar)

# ---- Data for plotting ----
allele_df <- tibble(
  pop = c("Pop 1", "Pop 2", "Pooled"),
  pA  = c(p1, p2, pbar),
  pa  = 1 - pA
) %>%
  pivot_longer(cols = c(pA, pa), names_to = "allele", values_to = "freq") %>%
  mutate(allele = recode(allele, pA = "A", pa = "a"))

geno_df <- tibble(
  genotype = c("AA","Aa","aa"),
  Observed_pooled = as.numeric(g_pool_obs),
  Expected_HWE    = as.numeric(g_pool_hwe)
) %>%
  pivot_longer(cols = c(Observed_pooled, Expected_HWE),
               names_to = "type", values_to = "freq") %>%
  mutate(type = recode(type,
                       Observed_pooled = "Observed (pooled)",
                       Expected_HWE    = "Expected (HWE from pooled p̄)"))

# ---- Panel 1: Allele frequencies ----
p1_plot <- ggplot(allele_df, aes(pop, freq, fill = allele)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "#d73027", "a" = "#4575b4")) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Allele frequencies",
       x = NULL, y = "Frequency") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# ---- Panel 2: Genotype frequencies (Wahlund effect) ----
p2_plot <- ggplot(geno_df, aes(genotype, freq, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "Pooled genotypes deviate from HWE (Wahlund effect)",
       subtitle = paste0("m = ", m, " | p̄ = ", round(pbar, 2)),
       x = NULL, y = "Frequency", fill = "") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

p1_plot | p2_plot


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(1)

# -------------------------
# Parameters you can tweak
# -------------------------
n_ind <- 80              # individuals
K <- 2

p1 <- 0.2                # allele A frequency in ancestry 1
p2 <- 0.8                # allele A frequency in ancestry 2

# admixture proportions per individual (q = ancestry-1 proportion)
# smaller "conc" => more extreme bars; larger => more similar individuals
mean_q <- 0.5
conc   <- 6
alpha  <- mean_q * conc
beta   <- (1 - mean_q) * conc

# -------------------------
# Simulate ADMIXTURE Q-matrix (K=2)
# -------------------------
q1 <- rbeta(n_ind, alpha, beta)     # ancestry 1 proportion
q2 <- 1 - q1

Q_long <- tibble(
  id = 1:n_ind,
  anc1 = q1,
  anc2 = q2
) %>%
  pivot_longer(cols = c(anc1, anc2), names_to = "ancestry", values_to = "q") %>%
  mutate(ancestry = recode(ancestry, anc1 = "Ancestry 1", anc2 = "Ancestry 2"))

# order individuals by ancestry-1 proportion (common in STRUCTURE plots)
id_order <- tibble(id = 1:n_ind, q1 = q1) %>% arrange(desc(q1)) %>% pull(id)
Q_long$id <- factor(Q_long$id, levels = id_order)

# -------------------------
# Simulate genotypes under local ancestry
# (each allele chooses ancestry ~ Bernoulli(q1); then allele ~ Bernoulli(p_ancestry))
# -------------------------
draw_allele <- function(q1_i) {
  anc <- rbinom(1, 1, q1_i)              # 1 => ancestry 1, 0 => ancestry 2
  pA  <- ifelse(anc == 1, p1, p2)
  rbinom(1, 1, pA)                       # 1 => A, 0 => a
}

geno <- sapply(q1, function(qi) {
  a1 <- draw_allele(qi)
  a2 <- draw_allele(qi)
  # genotype coded as number of A alleles: 0,1,2
  a1 + a2
})

# Observed pooled genotype counts
obs_counts <- tibble(
  genotype = c("AA","Aa","aa"),
  count = c(sum(geno == 2), sum(geno == 1), sum(geno == 0))
)

N <- sum(obs_counts$count)

# Pooled allele frequency estimate from pooled genotypes
p_hat <- (2*obs_counts$count[obs_counts$genotype=="AA"] +
            obs_counts$count[obs_counts$genotype=="Aa"]) / (2*N)
q_hat <- 1 - p_hat

# Expected under HWE using pooled p_hat
exp_counts <- tibble(
  genotype = c("AA","Aa","aa"),
  expected = c(p_hat^2, 2*p_hat*q_hat, q_hat^2) * N
)

df_hwe <- obs_counts %>%
  left_join(exp_counts, by="genotype") %>%
  pivot_longer(cols = c(count, expected), names_to = "type", values_to = "n") %>%
  mutate(type = recode(type, count = "Observed (pooled)", expected = "Expected (HWE)"))

# -------------------------
# Panel 1: ADMIXTURE / STRUCTURE barplot
# -------------------------
p_admix <- ggplot(Q_long, aes(x = id, y = q, fill = ancestry)) +
  geom_col(width = 1) +
  labs(
    title = "Admixture proportions per individual (K = 2)",
    x = NULL, y = "Ancestry proportion", fill = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

# -------------------------
# Panel 2: Pooled HWE deviation
# -------------------------
p_hwe <- ggplot(df_hwe, aes(genotype, n, fill = type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +
  labs(
    title = "Pooling admixed individuals can deviate from HWE",
    subtitle = paste0("Pooled p̂ = ", round(p_hat, 3), " (ancestry allele freqs: p1=", p1, ", p2=", p2, ")"),
    x = NULL, y = "Genotype count", fill = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

p_admix | p_hwe


# Selective sweep

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
# STATISTICS (same as Shiny app)
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
# EXACT SAME PLOTS AS SHINY (extracted)
# =========================================================
make_sweep_plots <- function(
    mode = c("Textbook Hard Sweep", "Realistic Hard Sweep"),
    n_hap = 100,
    L = 201,
    sweep_freq = 1.0,
    width = 30,          # core_width OR rec_scale
    sort_haps = TRUE,
    seed = NULL
) {
  mode <- match.arg(mode)
  
  L <- as.integer(L)
  if (L %% 2 == 0) L <- L + 1
  
  sweep_data <- if (mode == "Textbook Hard Sweep") {
    generate_textbook_sweep(
      n_hap = n_hap, L = L, sweep_freq = sweep_freq,
      core_width = width, seed = seed
    )
  } else {
    generate_realistic_sweep(
      n_hap = n_hap, L = L, sweep_freq = sweep_freq,
      rec_scale = width, seed = seed
    )
  }
  
  H <- sweep_data$H
  center <- sweep_data$center
  
  # Optional sorting to make sweep block visually clean (same as Shiny)
  if (sort_haps) {
    H <- H[order(H[, center], decreasing = TRUE), , drop = FALSE]
  }
  
  # ---- Haplotypes plot (same as Shiny) ----
  df_hap <- data.frame(
    Hap = rep(seq_len(nrow(H)), each = ncol(H)),
    Pos = rep(seq_len(ncol(H)), times = nrow(H)),
    Allele = as.vector(H)
  )
  df_hap$Allele <- as.vector(t(H))  # same final assignment as in Shiny app
  
  hapPlot <- ggplot(df_hap, aes(x = Pos, y = Hap, fill = factor(Allele))) +
    geom_tile() +
    geom_vline(xintercept = center, linetype = "dashed") +
    theme_minimal() +
    labs(fill = "Allele",
         title = "Haplotype structure around selected site")
  
  # ---- Gene diversity plot (same as Shiny) ----
  dist <- abs(seq_len(ncol(H)) - center)
  gd <- gene_diversity(H)
  df_gd <- data.frame(Distance = dist, GD = gd)
  
  gdPlot <- ggplot(df_gd, aes(x = Distance, y = GD)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    ylim(0, 0.5) +
    labs(title = "Diversity trough around selected site",
         y = "2p(1-p)")
  
  # ---- LD plot (same as Shiny) ----
  r2 <- ld_r2_with_center(H, center)
  df_ld <- data.frame(Distance = dist, r2 = r2)
  
  ldPlot <- ggplot(df_ld, aes(x = Distance, y = r2)) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    ylim(0, 1) +
    labs(title = "LD peak near selected site",
         y = expression(r^2))
  
  list(hapPlot = hapPlot, gdPlot = gdPlot, ldPlot = ldPlot)
}

# =========================================================
# EXAMPLE: run once and print plots
# =========================================================
plots <- make_sweep_plots(
  mode = "Realistic Hard Sweep",
  n_hap = 100,
  L = 201,
  sweep_freq = 0.9,
  width = 30,
  sort_haps = TRUE,
  seed = 1
)

plots$hapPlot
plots$gdPlot
plots$ldPlot




library(ggplot2)

# Parameters
N <- 200          # number of individuals
p <- 0.7          # allele frequency of A
set.seed(1)

# Gene pool (2N alleles)
n_alleles <- 2 * N
n_A <- round(p * n_alleles)
alleles <- c(rep("A", n_A), rep("a", n_alleles - n_A))
alleles <- sample(alleles)

# Pair alleles into diploid individuals
genotypes <- paste0(alleles[seq(1, n_alleles, 2)],
                    alleles[seq(2, n_alleles, 2)])

# Standardize heterozygotes
genotypes <- ifelse(genotypes %in% c("Aa", "aA"), "Aa", genotypes)

df_geno <- data.frame(genotype = genotypes)
df_allele <- data.frame(allele = alleles)


allele_freq <- as.data.frame(prop.table(table(df_allele$allele)))
colnames(allele_freq) <- c("allele", "frequency")

ggplot(allele_freq, aes(x = allele, y = frequency, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  ylim(0, 1) +
  labs(
    title = "Allele Frequencies in the Gene Pool",
    subtitle = paste("p =", round(allele_freq$frequency[allele_freq$allele=="A"], 2),
                     "| q =", round(allele_freq$frequency[allele_freq$allele=="a"], 2)),
    x = "Allele",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 16)


geno_freq <- as.data.frame(prop.table(table(df_geno$genotype)))
colnames(geno_freq) <- c("genotype", "frequency")

ggplot(geno_freq, aes(x = genotype, y = frequency, fill = genotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("AA" = "darkred",
                               "Aa" = "purple",
                               "aa" = "darkblue")) +
  ylim(0, 1) +
  labs(
    title = "Genotype Frequencies in the Population",
    subtitle = "Genotypes formed by random pairing of alleles",
    x = "Genotype",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 16)


library(ggplot2)
library(patchwork)

# Parameters
N <- 200
p <- 0.7
set.seed(1)

# Gene pool
n_alleles <- 2 * N
n_A <- round(p * n_alleles)
alleles <- c(rep("A", n_A), rep("a", n_alleles - n_A))
alleles <- sample(alleles)

df_allele <- data.frame(allele = alleles)

# ---- Counts ----
allele_counts <- as.data.frame(table(df_allele$allele))
colnames(allele_counts) <- c("allele", "count")

# ---- Frequencies ----
allele_freq <- as.data.frame(prop.table(table(df_allele$allele)))
colnames(allele_freq) <- c("allele", "frequency")

# ---- Plot 1: Counts ----
p_counts <- ggplot(allele_counts, aes(x = allele, y = count, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  labs(
    title = "Allele Counts (Gene Pool)",
    subtitle = paste("Total alleles = 2N =", 2*N),
    x = "Allele",
    y = "Count"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

# ---- Plot 2: Frequencies ----
p_freq <- ggplot(allele_freq, aes(x = allele, y = frequency, fill = allele)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  ylim(0, 1) +
  labs(
    title = "Allele Frequencies",
    subtitle = paste("p =", round(allele_freq$frequency[allele_freq$allele=="A"],2),
                     "| q =", round(allele_freq$frequency[allele_freq$allele=="a"],2)),
    x = "Allele",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

# ---- Combine side-by-side ----
p_counts | p_freq


library(ggplot2)
#install.packages("gganimate")
library(gganimate)

# Parameters
N <- 20
gens <- 40
p0 <- 0.5
set.seed(2)

n_alleles <- 2*N
p <- numeric(gens)
p[1] <- p0

# Simulate drift
for (t in 2:gens) {
  p[t] <- rbinom(1, n_alleles, p[t-1]) / n_alleles
}

df <- data.frame(
  generation = 1:gens,
  p = p
)

ggplot(df, aes(generation, p)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, color = "red") +
  ylim(0,1) +
  labs(
    title = "Genetic Drift",
    subtitle = "Generation: {closest_state}",
    y = "Allele frequency (p)"
  ) +
  theme_minimal(base_size = 16) +
  transition_states(generation, transition_length = 1, state_length = 1)

animate(last_plot(), width = 900, height = 500, fps = 5)

library(ggplot2)
library(dplyr)

# --------------------------
# Parameters
# --------------------------
N <- 8                 # number of diploid individuals
p <- 0.6               # allele A frequency
set.seed(1)

# Total alleles
n_alleles <- 2 * N
n_A <- round(p * n_alleles)
n_a <- n_alleles - n_A

# Create allele vector
alleles <- c(rep("A", n_A), rep("a", n_a))
alleles <- sample(alleles)   # shuffle

# Assign alleles to individuals (2 per individual)
df <- data.frame(
  individual = rep(1:N, each = 2),
  allele_copy = rep(1:2, times = N),
  allele = alleles
)

# Plot
ggplot(df, aes(x = allele_copy, y = individual, fill = allele)) +
  geom_point(shape = 21, size = 8, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("A" = "red", "a" = "blue")) +
  scale_y_reverse() +
  labs(
    title = "The Gene Pool (Diploid Population)",
    subtitle = paste("N =", N, "individuals  |  2N =", 2*N, "alleles  |  p =", round(n_A/(2*N),2)),
    x = "Allele copy",
    y = "Individual",
    fill = "Allele"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Fitness Distribution Analysis by Mutation Type
# This script analyzes and visualizes fitness distributions for different mutation types

library(dplyr)
library(ggplot2)

#' Plot fitness density distributions by mutation type
#'
#' @param fit_data Fitness data for all variants
#' @param wt_data Wild-type fitness data
#' @param syn_data Synonymous mutation fitness data
#' @return ggplot object showing density distributions

plot_fitness_density <- function(fit_data, wt_data, syn_data) {
  
  # WT variants
  df_wt <- wt_data %>%
    mutate(mut_type = "WT")
  
  # Synonymous mutations
  df_syn <- syn_data %>%
    mutate(mut_type = "Synonymous")
  
  # Stop mutations (STOP or readthrough)
  df_stop <- fit_data %>%
    filter(STOP == TRUE | STOP_readthrough == TRUE) %>%
    mutate(mut_type = "Stop")
  
  # Missense mutations (Nham_aa > 0 and not indel/STOP)
  df_mis <- fit_data %>%
    filter(Nham_aa > 0 & indel == FALSE & STOP == FALSE) %>%
    mutate(mut_type = "Missense")
  
  # Combine all data
  df_plot <- bind_rows(df_wt, df_syn, df_stop, df_mis) %>%
    select(mut_type, fitness)
  
  # Create plot
  p <- ggplot(df_plot, aes(x = fitness, fill = mut_type, color = mut_type)) +
    geom_density(alpha = 0.4) +  # Semi-transparent fill
    scale_fill_manual(values = c(
      "WT" = "black",
      "Synonymous" = "#1B38A6",
      "Missense" = "#F1DD10",
      "Stop" = "#F4270C"
    )) +
    scale_color_manual(values = c(
      "WT" = "black",
      "Synonymous" = "#1B38A6",
      "Missense" = "#F1DD10",
      "Stop" = "#F4270C"
    )) +
    labs(x = "Fitness", y = "Density", title = "Fitness Density Distribution by Mutation Type") +
    theme_classic(base_size = 8) +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8)
    )
  
  return(p)
}

# =====================
# Analysis for Different Assays and Blocks
# =====================

# Define wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# K13 Analysis
load("path/to/K13_block1.RData")
p_k13_b1 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K13_block1_density.pdf", p_k13_b1, width = 4, height = 4, dpi = 300)

load("path/to/K13_block2.RData")
p_k13_b2 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K13_block2_density.pdf", p_k13_b2, width = 4, height = 4, dpi = 300)

load("path/to/K13_block3.RData")
p_k13_b3 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K13_block3_density.pdf", p_k13_b3, width = 4, height = 4, dpi = 300)

# K19 Analysis
load("path/to/K19_block1.RData")
p_k19_b1 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K19_block1_density.pdf", p_k19_b1, width = 4, height = 4, dpi = 300)

load("path/to/K19_block2.RData")
p_k19_b2 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K19_block2_density.pdf", p_k19_b2, width = 4, height = 4, dpi = 300)

load("path/to/K19_block3.RData")
p_k19_b3 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/K19_block3_density.pdf", p_k19_b3, width = 4, height = 4, dpi = 300)

# RAF1 Analysis
load("path/to/RAF1_block1.RData")
p_raf1_b1 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/RAF1_block1_density.pdf", p_raf1_b1, width = 4, height = 4, dpi = 300)

load("path/to/RAF1_block2.RData")
p_raf1_b2 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/RAF1_block2_density.pdf", p_raf1_b2, width = 4, height = 4, dpi = 300)

load("path/to/RAF1_block3.RData")
p_raf1_b3 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/RAF1_block3_density.pdf", p_raf1_b3, width = 4, height = 4, dpi = 300)

# Abundance Analysis
load("path/to/Abundance_block1.RData")
p_abundance_b1 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/Abundance_block1_density.pdf", p_abundance_b1, width = 4, height = 4, dpi = 300)

load("path/to/Abundance_block2.RData")
p_abundance_b2 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/Abundance_block2_density.pdf", p_abundance_b2, width = 4, height = 4, dpi = 300)

load("path/to/Abundance_block3.RData")
p_abundance_b3 <- plot_fitness_density(
  fit_data = all_variants %>% filter(sigma < 0.5),
  wt_data = wildtype %>% filter(sigma < 0.5),
  syn_data = synonymous %>% filter(sigma < 0.5)
)
ggsave("path/to/Abundance_block3_density.pdf", p_abundance_b3, width = 4, height = 4, dpi = 300)
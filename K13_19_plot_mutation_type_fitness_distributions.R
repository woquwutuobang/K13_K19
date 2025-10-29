library(ggplot2)
library(dplyr)
library(patchwork)
library(wlab.block)

# ============================================================
# Function: plot_assay_fitness_density
# Description:
#   Generate density plots of normalized fitness distributions for 
#   synonymous, missense, and stop mutations across experimental blocks
#   for any given assay type (abundance, binding, etc.).
#   The top panel shows the global distribution across all blocks,
#   while the bottom panel shows block-specific distributions.
# 
# Parameters:
#   assay_type: Name of the assay (e.g., "abundance", "RAF1", "SOS1", etc.)
#   block1: Path to block1 fitness data file
#   block2: Path to block2 fitness data file  
#   block3: Path to block3 fitness data file
#   output_file: Optional path to save the plot (if NULL, only displays plot)
#
# Returns:
#   Combined ggplot object with two panels
# ============================================================

plot_fitness_density <- function(assay_type, block1, block2, block3, output_file = NULL) {
  
  # Load normalized fitness data from three experimental blocks
  nor_fit <- nor_fitness(block1 = block1, block2 = block2, block3 = block3)
  
  # Classify mutation types based on sequence changes
  nor_fit_classified <- nor_fit %>%
    mutate(
      mut_type = case_when(
        # Synonymous: amino acid unchanged but nucleotide changed
        Nham_aa == 0 & Nham_nt > 0 ~ "Synonymous",
        # Stop: introduces stop codon
        STOP == TRUE | STOP_readthrough == TRUE ~ "Stop",
        # Missense: amino acid changed, not indel, not stop
        Nham_aa > 0 & indel == FALSE & STOP == FALSE & STOP_readthrough == FALSE ~ "Missense"
      )
    ) %>%
    filter(!is.na(mut_type))  # Remove unclassified mutations
  
  # Display mutation type distribution
  cat("Mutation type distribution for", assay_type, ":\n")
  print(table(nor_fit_classified$mut_type))
  
  # Prepare data for plotting
  nor_fit_plot <- nor_fit_classified %>%
    mutate(mut_type = factor(mut_type, levels = c("Synonymous", "Missense", "Stop")))
  
  # Display distribution for plotting data
  cat("Mutation types for plotting", assay_type, ":\n")
  print(table(nor_fit_plot$mut_type))
  
  # Create global density plot (all blocks combined)
  p_overall <- ggplot(nor_fit_plot, aes(x = nor_fitness, color = mut_type)) +
    geom_density(size = 1) +  # Line plot instead of fill
    scale_color_manual(
      values = c(
        "Synonymous" = "#09B636", 
        "Missense" = "#F4AD0C",
        "Stop" = "#FF6A56"
      ),
      name = "Mutation Type"
    ) +
    labs(
      title = paste(toupper(assay_type), "- Fitness Distribution by Mutation Type"),
      x = "Normalized Fitness",
      y = "Density"
    ) +
    # Set x-axis limits
    xlim(-1.5, 0.5) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 8),
      text = element_text(size = 8),
      axis.text = element_text(size = 8)
    )
  
  # Create block-specific density plots
  p_by_block <- ggplot(nor_fit_plot, aes(x = nor_fitness, color = mut_type)) +
    geom_density(size = 1) +  # Line plot instead of fill
    facet_wrap(~ block, ncol = 4) +
    scale_color_manual(
      values = c(
        "Synonymous" = "#09B636",
        "Missense" = "#F4AD0C", 
        "Stop" = "#FF6A56"
      )
    ) +
    labs(
      title = paste(toupper(assay_type), "- Fitness Distribution by Block"),
      x = "Normalized Fitness", 
      y = "Density"
    ) +
    # Set x-axis limits
    xlim(-1.5, 0.5) +
    theme_classic() +
    theme(
      legend.position = "none",  # Remove legend from this panel
      plot.title = element_text(hjust = 0.5, size = 8),
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8)
    )
  
  # Combine plots with shared legend
  combined_plot <- p_overall / p_by_block + 
    plot_layout(heights = c(1, 2), guides = "collect") &
    theme(legend.position = "bottom")
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot, width = 4, height = 6, units = "in")
    cat("Plot saved to:", output_file, "\n")
  }
  
  # Return the combined plot object
  return(combined_plot)
}



### Abundance

plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "Abundance",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/Abundance_normalized_density.pdf"
)



### K13
plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "K13",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/K13_normalized_density.pdf"
)



### K19
plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "K19",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/K19_normalized_density.pdf"
)

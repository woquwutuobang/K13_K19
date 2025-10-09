# ================================
# Fitness Distribution Analysis by Mutation Type
# ================================

library(dplyr)
library(ggplot2)
library(patchwork)

#' Plot fitness density distributions by mutation type
#'
#' @param fit_data Fitness data for all variants
#' @param wt_data Wild-type fitness data
#' @param syn_data Synonymous mutation fitness data
#' @param assay_name Assay name (for title annotation)
#' @param block_name Block name (for subtitle annotation)
#' @return ggplot object showing density distributions
plot_fitness_density <- function(fit_data, wt_data, syn_data, assay_name, block_name) {
  
  # WT
  df_wt <- wt_data %>% mutate(mut_type = "WT")
  
  # Synonymous
  df_syn <- syn_data %>% mutate(mut_type = "Synonymous")
  
  # Stop
  df_stop <- fit_data %>%
    filter(STOP == TRUE | STOP_readthrough == TRUE) %>%
    mutate(mut_type = "Stop")
  
  # Missense
  df_mis <- fit_data %>%
    filter(Nham_aa > 0 & indel == FALSE & STOP == FALSE) %>%
    mutate(mut_type = "Missense")
  
  # Combine all data
  df_plot <- bind_rows(df_wt, df_syn, df_stop, df_mis) %>%
    select(mut_type, fitness)
  
  # 绘图
  p <- ggplot(df_plot, aes(x = fitness, fill = mut_type, color = mut_type)) +
    geom_density(alpha = 0.4) +
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
    labs(
      x = "Fitness",
      y = "Density",
      title = assay_name,
      subtitle = block_name
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5)
    )
  
  return(p)
}

# =====================
# Helper function: combine three blocks for one assay
# =====================
combine_assay_blocks <- function(assay_name, block_paths, output_file) {
  
  plots <- list()
  
  for (i in seq_along(block_paths)) {
    load(block_paths[[i]])  # loads all_variants, wildtype, synonymous
    block_name <- paste0("Block ", i)
    
    plots[[i]] <- plot_fitness_density(
      fit_data = all_variants %>% filter(sigma < 0.5),
      wt_data = wildtype %>% filter(sigma < 0.5),
      syn_data = synonymous %>% filter(sigma < 0.5),
      assay_name = assay_name,
      block_name = block_name
    )
  }
  
  combined_plot <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = paste0("Fitness Density Distribution - ", assay_name),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 8))
    )
  
  ggsave(output_file, combined_plot, width = 12, height = 5, dpi = 300)
  cat(paste("Saved combined figure for", assay_name, "→", output_file, "\n"))
  
  return(combined_plot)
}

# =====================
# Example usage
# =====================
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Define paths (replace with actual file paths)
k13_paths <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData"
)
k19_paths <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData"
)
raf1_paths <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData"
)
abundance_paths <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData"
)

# Generate combined figures
combine_assay_blocks("K13", k13_paths, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/K13_combined_density.pdf")
combine_assay_blocks("K19", k19_paths, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/K19_combined_density.pdf")
combine_assay_blocks("RAF1", raf1_paths, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/RAF1_combined_density.pdf")
combine_assay_blocks("Abundance", abundance_paths, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/Abundance_combined_density.pdf")

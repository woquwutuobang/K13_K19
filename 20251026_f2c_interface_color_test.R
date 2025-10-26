# KRAS Single Pair Binding Comparison
# This script creates a scatterplot comparing mutation effects between two binding assays

library(data.table)
library(krasddpcams)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)

#' Create scatterplot for single pair ΔΔG comparison
#'
#' @param ddG_path1 File path to first ddG data
#' @param ddG_path2 File path to second ddG data
#' @param assay_name1 First assay name
#' @param assay_name2 Second assay name
#' @param keep_cols Columns to keep (default 23:27 for mean/std columns)
#' @param interface_sets List defining interface residue sets
#' @param color_map Color mapping for interface groups
#' @param point_size Point size (default 1.1)
#' @param alpha Point transparency (default 0.7)
#' @param base_size Base font size (default 10)
#' @param xlim X-axis limits (default c(-3.3, 3.3))
#' @param ylim Y-axis limits (default c(-3.3, 3.3))
#' @return A ggplot object

create_single_assay_comparison <- function(ddG_path1, 
                                           ddG_path2,
                                           assay_name1,
                                           assay_name2,
                                           keep_cols = 23:27,
                                           interface_sets = list(
                                             Binding_Interface_2 = c(107, 101, 102, 99, 136, 68, 95, 137, 94, 133, 90, 129, 87, 91, 88, 98),
                                             Binding_Interface_1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
                                             nucleotide_binding_pocket = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 
                                                                           116, 117, 119, 120, 145, 146, 147)
                                           ),
                                           color_map = c(
                                             "Other" = "grey",
                                             "Nucleotide Binding Pocket" = "#007A20",
                                             "Binding Interface 1" = "#F4270C",
                                             "Binding Interface 2" = "#1B38A6"
                                           ),
                                           point_size = 2,
                                           alpha = 0.5,
                                           base_size = 10,
                                           xlim = c(-1.5, 3.3),
                                           ylim = c(-1.5, 3.3)) {
  
  # Process single assay data
  process_assay_data <- function(path, assay_name) {
    # Read ddG data as data frame
    ddG <- as.data.frame(krasddpcams__read_ddG(ddG = path, assay_sele = assay_name))
    
    # Select required columns and remove NA values
    ddG_sel <- ddG %>% 
      select(mt, Pos_real, assay, `mean_kcal/mol`) %>% 
      filter(!is.na(`mean_kcal/mol`))
    
    # Convert to wide format
    ddG_wide <- ddG_sel %>% 
      pivot_wider(
        names_from = assay, 
        values_from = `mean_kcal/mol`,
        values_fn = mean  # Take mean if duplicate values exist
      )
    
    return(ddG_wide)
  }
  
  # Process both assay datasets
  assay_data1 <- process_assay_data(ddG_path1, assay_name1)
  assay_data2 <- process_assay_data(ddG_path2, assay_name2)
  
  # Merge two assay datasets
  merged <- merge(assay_data1, assay_data2,
                  by = c("mt", "Pos_real"),
                  all = TRUE,
                  suffixes = paste0("_", c(assay_name1, assay_name2)))
  
  # Add interface grouping information
  merged <- merged %>%
    mutate(interface_group = case_when(
      Pos_real %in% interface_sets$Binding_Interface_2 ~ "Binding Interface 2",
      Pos_real %in% interface_sets$Binding_Interface_1 ~ "Binding Interface 1",
      Pos_real %in% interface_sets$nucleotide_binding_pocket ~ "Nucleotide Binding Pocket",
      TRUE ~ "Other"
    ))
  
  # Reorder data for proper layering
  merged <- merged %>%
    mutate(plot_order = case_when(
      interface_group == "Other" ~ 1,
      interface_group == "Nucleotide Binding Pocket" ~ 2,
      interface_group == "Binding Interface 2" ~ 3,
      interface_group == "Binding Interface 1" ~ 4
    )) %>%
    arrange(plot_order)
  
  # Create scatter plot
  p <- ggplot(merged, aes(x = .data[[assay_name1]], y = .data[[assay_name2]], 
                          color = interface_group)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = color_map, name = "Interface Region") +
    labs(x = bquote(Delta*Delta*"Gb ("*.(assay_name1)*") (kcal/mol)"), 
         y = bquote(Delta*Delta*"Gb ("*.(assay_name2)*") (kcal/mol)")) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = base_size - 2),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  return(p)
}

# Example usage: Compare K13 and RAF1
comparison_plot <- create_single_assay_comparison(
  ddG_path1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  ddG_path2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_name1 = "K13",
  assay_name2 = "RAF1",
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3)
)

# Display plot
print(comparison_plot)

# Save plot
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K13_vs_RAF1_comparison.png", 
       plot = comparison_plot, 
       width = 4, 
       height = 4.5,
       dpi = 300)

# Save as PDF
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K13_vs_RAF1_comparison.pdf", 
       plot = comparison_plot, 
       width = 4, 
       height = 4.5, 
       device = cairo_pdf)

# You can also compare other pairs, for example:
# K13 vs K19
k13_k19_plot <- create_single_assay_comparison(
  ddG_path1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  ddG_path2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay_name1 = "K13",
  assay_name2 = "K19",
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3)
)

k13_k19_plot

# Save plot
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K13_vs_K19_comparison.png", 
       plot = k13_k19_plot, 
       width = 4, 
       height = 4.5,
       dpi = 300)

# Save as PDF
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K13_vs_K19_comparison.pdf", 
       plot = k13_k19_plot, 
       width = 4, 
       height = 4.5, 
       device = cairo_pdf)



# RAF1 vs K27
raf1_k27_plot <- create_single_assay_comparison(
  ddG_path1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  ddG_path2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_name1 = "K27",
  assay_name2 = "RAF1",
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3)
)


raf1_k27_plot

# Save plot
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K27_vs_RAF1_comparison.png", 
       plot = raf1_k27_plot, 
       width = 4, 
       height = 4.5,
       dpi = 300)

# Save as PDF
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/K27_vs_RAF1_comparison.pdf", 
       plot = raf1_k27_plot, 
       width = 4, 
       height = 4.5, 
       device = cairo_pdf)

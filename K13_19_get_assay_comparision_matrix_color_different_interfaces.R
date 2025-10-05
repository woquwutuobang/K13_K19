# KRAS Multi-Assay Binding Comparison
# This script creates a scatterplot matrix comparing mutation effects across multiple binding assays

library(data.table)
library(krasddpcams)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(patchwork)

#' Create scatterplot matrix for multi-assay ΔΔG comparisons
#'
#' @param ddG_paths Vector of file paths to ddG data
#' @param assay_names Vector of assay names
#' @param keep_cols Columns to keep (default 23:27 for mean/std columns)
#' @param interface_sets List defining interface residue sets
#' @param color_map Color mapping for interface groups
#' @param point_size Point size (default 1.1)
#' @param alpha Point transparency (default 0.7)
#' @param base_size Base font size (default 10)
#' @param xlim X-axis limits (default c(-3.3, 3.3))
#' @param ylim Y-axis limits (default c(-3.3, 3.3))
#' @return A patchwork combined plot

create_assay_comparison_matrix <- function(ddG_paths, 
                                           assay_names,
                                           keep_cols = 23:27,
                                           interface_sets = list(
                                             allosteric_interface = c(107, 101, 102, 99, 136, 68, 95, 137, 94, 133, 90, 129, 87, 91, 88, 98),
                                             effector_interface = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
                                             nucleotide_binding_pocket = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 
                                                                           116, 117, 119, 120, 145, 146, 147)
                                           ),
                                           color_map = c(
                                             "Other" = "grey",
                                             "Nucleotide Binding Pocket" = "#F1DD10",
                                             "K13 Binding Interface" = "#D62728",
                                             "RAF1 Binding Interface" = "#0C226F"
                                           ),
                                           point_size = 1.1,
                                           alpha = 0.7,
                                           base_size = 10,
                                           xlim = c(-3.3, 3.3),
                                           ylim = c(-3.3, 3.3)) {
  
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
  
  # Process all assay data
  assay_data <- map2(ddG_paths, assay_names, process_assay_data)
  
  # Create all pairwise combinations
  combinations <- combn(length(assay_names), 2, simplify = FALSE)
  
  # Generate plot list
  plot_list <- map(combinations, function(pair) {
    i <- pair[1]
    j <- pair[2]
    
    # Merge two assay datasets
    merged <- merge(assay_data[[i]], assay_data[[j]],
                    by = c("mt", "Pos_real"),
                    all = TRUE,
                    suffixes = paste0("_", assay_names[c(i, j)]))
    
    # Add interface grouping information
    merged <- merged %>%
      mutate(interface_group = case_when(
        Pos_real %in% interface_sets$allosteric_interface ~ "K13 Binding Interface",
        Pos_real %in% interface_sets$effector_interface ~ "RAF1 Binding Interface",
        Pos_real %in% interface_sets$nucleotide_binding_pocket ~ "Nucleotide Binding Pocket",
        TRUE ~ "Other"
      ))
    
    # Reorder data for proper layering
    merged <- merged %>%
      mutate(plot_order = case_when(
        interface_group == "Other" ~ 1,
        interface_group == "Nucleotide Binding Pocket" ~ 2,
        interface_group == "K13 Binding Interface" ~ 3,
        interface_group == "RAF1 Binding Interface" ~ 4
      )) %>%
      arrange(plot_order)
    
    # Get current assay names
    x_assay <- assay_names[i]
    y_assay <- assay_names[j]
    
    # Create scatter plot
    ggplot(merged, aes(x = .data[[x_assay]], y = .data[[y_assay]], 
                       color = interface_group)) +
      geom_point(size = point_size, alpha = alpha) +
      scale_color_manual(values = color_map) +
      labs(x = bquote(Delta*Delta*"G ("*.(x_assay)*")"), 
           y = bquote(Delta*Delta*"G ("*.(y_assay)*")")) +
      theme_classic(base_size = base_size) +
      theme(
        legend.position = "none",
        axis.text = element_text(size = base_size - 2),
        axis.title = element_text(size = base_size - 1)
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
      coord_cartesian(xlim = xlim, ylim = ylim)
  })
  
  # Create legend plot
  legend_plot <- ggplot(data.frame(x = 0, y = 0, 
                                   interface_group = factor(names(color_map), 
                                                            levels = c("Other", "Nucleotide Binding Pocket", 
                                                                       "K13 Binding Interface", "RAF1 Binding Interface")))) +
    geom_point(aes(x, y, color = interface_group)) +
    scale_color_manual(values = color_map, name = "Interface Region") +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = base_size - 2),
          legend.title = element_text(size = base_size - 1))
  
  # Combine plots with legend
  wrap_plots(plot_list) / guide_area() + 
    plot_layout(guides = "collect", heights = c(20, 1)) &
    theme(legend.position = "bottom")
}

# Define input parameters
ddG_paths <- c(
  "path/to/weights_Binding_K13.txt",
  "path/to/weights_Binding_RAF1.txt", 
  "path/to/weights_Binding_K19.txt",
  "path/to/weights_Binding_K27.txt",
  "path/to/weights_Binding_K55.txt",
  "path/to/weights_Binding_PI3.txt",
  "path/to/weights_Binding_RAL.txt",
  "path/to/weights_Binding_SOS.txt"
)

assay_names <- c("K13", "RAF1", "K19", "K27", "K55", "PI3KCG", "RALGDS", "SOS1")

# Generate comparison plot with fixed axis limits
comparison_plot <- create_assay_comparison_matrix(
  ddG_paths, 
  assay_names,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3)
)

# Display plot
print(comparison_plot)

# Save plot
ggsave("KRAS_multi_assay_comparison.png", 
       plot = comparison_plot, 
       width = 14.5, 
       height = 12,
       dpi = 100)

# Optional: Save as PDF
# ggsave("KRAS_multi_assay_comparison.pdf", 
#        plot = comparison_plot, 
#        width = 14.5, 
#        height = 12, 
#        device = cairo_pdf)
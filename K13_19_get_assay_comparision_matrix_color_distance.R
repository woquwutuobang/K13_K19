# KRAS Multi-Assay Binding Comparison Matrix
# This script creates a comprehensive matrix of pairwise scatter plots comparing mutation effects across KRAS binding partners

library(data.table)
library(krasddpcams)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(patchwork)
library(purrr)

# =====================
# Configuration Parameters
# =====================

# Define file paths for binding assay data
assay_paths <- list(
  K13    = "path/to/weights_Binding_K13.txt",
  RAF1   = "path/to/weights_Binding_RAF1.txt", 
  K19    = "path/to/weights_Binding_K19.txt",
  K27    = "path/to/weights_Binding_K27.txt",
  K55    = "path/to/weights_Binding_K55.txt",
  PI3KCG = "path/to/weights_Binding_PI3.txt",
  RALGDS = "path/to/weights_Binding_RAL.txt",
  SOS1   = "path/to/weights_Binding_SOS.txt"
)

# Output directory for results
output_directory <- "results/binding_comparison_matrix"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# =====================
# Data Loading Functions
# =====================

#' Load and process binding assay data
#'
#' @param file_path Path to the assay data file
#' @param assay_name Name of the binding assay
#' @return Processed data table in wide format

load_assay_data <- function(file_path, assay_name) {
  # Read raw ddG data
  raw_data <- krasddpcams__read_ddG(ddG = file_path, assay_sele = assay_name)
  
  # Select relevant columns and convert to wide format
  processed_data <- raw_data[, c(1:3, 23:27)] %>%
    spread(key = assay, value = `mean_kcal/mol`)
  
  return(processed_data)
}

# Load all assay datasets
assay_datasets <- map2(assay_paths, names(assay_paths), load_assay_data)
names(assay_datasets) <- names(assay_paths)

# Load structural annotation data
structural_annotations <- fread("path/to/structural_annotations.csv")
structural_annotations[, Pos_real := Pos]

# =====================
# Visualization Functions
# =====================

#' Create distance-colored scatter plot for assay comparison
#'
#' @param plot_data Data frame containing merged assay data
#' @param x_assay Name of the assay for x-axis
#' @param y_assay Name of the assay for y-axis  
#' @param color_column Column name for distance coloring
#' @param legend_title Title for the color legend
#' @param x_limits Limits for x-axis (default: c(-1.5, 3.3))
#' @param y_limits Limits for y-axis (default: c(-1.5, 3.3))
#' @return ggplot object with colored scatter plot

create_distance_scatter <- function(plot_data, x_assay, y_assay, color_column, 
                                    legend_title, x_limits = c(-1.5, 3.3), 
                                    y_limits = c(-1.5, 3.3)) {
  
  ggplot(plot_data, aes_string(x = x_assay, y = y_assay, color = color_column)) +
    geom_point(size = 1.3, alpha = 0.8) +
    scale_color_gradient(
      low = "#09B636",      # Green for close distances
      high = "#F4270C",     # Red for far distances
      name = legend_title
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(
      x = bquote(Delta*Delta*"G ("*.(x_assay)*")"),
      y = bquote(Delta*Delta*"G ("*.(y_assay)*")")
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6)
    )
}

# =====================
# Generate Comparison Matrix
# =====================

# Initialize list to store all comparison plots
comparison_plots <- list()

# Generate all pairwise combinations
assay_names <- names(assay_datasets)

for (i in 1:(length(assay_names) - 1)) {
  for (j in (i + 1):length(assay_names)) {
    
    current_assay1 <- assay_names[i]
    current_assay2 <- assay_names[j]
    
    # Merge data from both assays
    merged_data <- merge(
      assay_datasets[[i]], 
      assay_datasets[[j]],
      by = c("mt", "Pos_real"),
      all = TRUE,
      suffixes = c(paste0("_", current_assay1), paste0("_", current_assay2))
    )
    
    # Add structural annotations
    merged_data <- merge(merged_data, structural_annotations, by = "Pos_real", all = TRUE)
    
    # Convert distance columns to numeric
    distance_col1 <- paste0("scHAmin_ligand_", current_assay1)
    distance_col2 <- paste0("scHAmin_ligand_", current_assay2)
    merged_data[[distance_col1]] <- as.numeric(merged_data[[distance_col1]])
    merged_data[[distance_col2]] <- as.numeric(merged_data[[distance_col2]])
    
    # Filter complete cases
    clean_data <- merged_data %>%
      filter(
        !is.na(.data[[current_assay1]]), 
        !is.na(.data[[current_assay2]]),
        !is.na(.data[[distance_col1]]),
        !is.na(.data[[distance_col2]])
      )
    
    # Create plots colored by distance to each partner
    plot_distance1 <- create_distance_scatter(
      clean_data, current_assay1, current_assay2,
      distance_col1,
      paste("Distance to", current_assay1)
    )
    
    plot_distance2 <- create_distance_scatter(
      clean_data, current_assay1, current_assay2, 
      distance_col2,
      paste("Distance to", current_assay2)
    )
    
    # Combine both distance views
    combined_plot <- plot_distance1 + plot_distance2 + 
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    # Store in list
    plot_key <- paste0(current_assay1, "_vs_", current_assay2)
    comparison_plots[[plot_key]] <- combined_plot
  }
}

# =====================
# Create and Save Matrix Plot
# =====================

# Arrange all pairwise comparisons in a grid
matrix_plot <- wrap_plots(comparison_plots, ncol = 4)

# Save the comprehensive comparison matrix
ggsave(
  file.path(output_directory, "binding_assay_comparison_matrix.png"),
  plot = matrix_plot,
  width = 22, 
  height = 22,
  dpi = 100
)






# KRAS Pairwise Binding Comparison - Individual Plots
# This script generates individual pairwise scatter plots comparing mutation effects across KRAS binding partners

library(data.table)
library(krasddpcams)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(patchwork)
library(purrr)

# =====================
# Configuration Parameters
# =====================

# Define file paths for binding assay data
assay_paths <- list(
  K13    = "path/to/weights_Binding_K13.txt",
  RAF1   = "path/to/weights_Binding_RAF1.txt",
  K19    = "path/to/weights_Binding_K19.txt", 
  K27    = "path/to/weights_Binding_K27.txt",
  K55    = "path/to/weights_Binding_K55.txt",
  PI3KCG = "path/to/weights_Binding_PI3.txt",
  RALGDS = "path/to/weights_Binding_RAL.txt",
  SOS1   = "path/to/weights_Binding_SOS.txt"
)

# Output directory for individual plots
output_directory <- "results/individual_binding_comparisons"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# =====================
# Data Loading Functions
# =====================

#' Load and process binding assay data
#'
#' @param file_path Path to the assay data file
#' @param assay_name Name of the binding assay
#' @return Processed data table in wide format

load_assay_data <- function(file_path, assay_name) {
  # Read raw ddG data
  raw_data <- krasddpcams__read_ddG(ddG = file_path, assay_sele = assay_name)
  
  # Select relevant columns and convert to wide format
  processed_data <- raw_data[, c(1:3, 23:27)] %>%
    spread(key = assay, value = `mean_kcal/mol`)
  
  return(processed_data)
}

# Load all assay datasets
assay_datasets <- map2(assay_paths, names(assay_paths), load_assay_data)
names(assay_datasets) <- names(assay_paths)

# Load structural annotation data
structural_annotations <- fread("path/to/structural_annotations.csv")
structural_annotations[, Pos_real := Pos]

# =====================
# Visualization Functions
# =====================

#' Create distance-colored scatter plot for assay comparison
#'
#' @param plot_data Data frame containing merged assay data
#' @param x_assay Name of the assay for x-axis
#' @param y_assay Name of the assay for y-axis
#' @param color_column Column name for distance coloring
#' @param legend_title Title for the color legend
#' @param x_limits Limits for x-axis (default: c(-1.5, 3.3))
#' @param y_limits Limits for y-axis (default: c(-1.5, 3.3))
#' @return ggplot object with colored scatter plot

create_distance_scatter <- function(plot_data, x_assay, y_assay, color_column, 
                                    legend_title, x_limits = c(-1.5, 3.3), 
                                    y_limits = c(-1.5, 3.3)) {
  
  ggplot(plot_data, aes_string(x = x_assay, y = y_assay, color = color_column)) +
    geom_point(size = 1.3, alpha = 0.8) +
    scale_color_gradient(
      low = "#09B636",      # Green for close distances
      high = "#F4270C",     # Red for far distances
      name = legend_title
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(
      x = bquote(Delta*Delta*"G ("*.(x_assay)*")"),
      y = bquote(Delta*Delta*"G ("*.(y_assay)*")")
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
}

# =====================
# Generate Individual Pairwise Comparisons
# =====================

# Define consistent axis limits for all plots
x_axis_limits <- c(-1.5, 3.3)
y_axis_limits <- c(-1.5, 3.3)

assay_names <- names(assay_datasets)

# Generate and save all pairwise comparisons
for (i in 1:(length(assay_names) - 1)) {
  for (j in (i + 1):length(assay_names)) {
    
    current_assay1 <- assay_names[i]
    current_assay2 <- assay_names[j]
    
    # Merge data from both assays
    merged_data <- merge(
      assay_datasets[[i]], 
      assay_datasets[[j]],
      by = c("mt", "Pos_real"),
      all = TRUE,
      suffixes = c(paste0("_", current_assay1), paste0("_", current_assay2))
    )
    
    # Add structural annotations
    merged_data <- merge(merged_data, structural_annotations, by = "Pos_real", all = TRUE)
    
    # Convert distance columns to numeric
    distance_col1 <- paste0("scHAmin_ligand_", current_assay1)
    distance_col2 <- paste0("scHAmin_ligand_", current_assay2)
    merged_data[[distance_col1]] <- as.numeric(merged_data[[distance_col1]])
    merged_data[[distance_col2]] <- as.numeric(merged_data[[distance_col2]])
    
    # Filter complete cases
    clean_data <- merged_data %>%
      filter(
        !is.na(.data[[current_assay1]]), 
        !is.na(.data[[current_assay2]]),
        !is.na(.data[[distance_col1]]),
        !is.na(.data[[distance_col2]])
      )
    
    # Create plots colored by distance to each partner with consistent axis limits
    plot_distance1 <- create_distance_scatter(
      clean_data, current_assay1, current_assay2,
      distance_col1,
      paste("Distance to", current_assay1),
      x_limits = x_axis_limits,
      y_limits = y_axis_limits
    )
    
    plot_distance2 <- create_distance_scatter(
      clean_data, current_assay1, current_assay2, 
      distance_col2,
      paste("Distance to", current_assay2),
      x_limits = x_axis_limits,
      y_limits = y_axis_limits
    )
    
    # Combine both distance views
    combined_plot <- plot_distance1 + plot_distance2 + 
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    # Save individual plot as PDF
    output_filename <- paste0("comparison_", current_assay1, "_vs_", current_assay2, ".pdf")
    output_path <- file.path(output_directory, output_filename)
    
    ggsave(
      output_path,
      plot = combined_plot,
      width = 6,
      height = 4,
      device = cairo_pdf
    )
    
    # Print progress
    cat("Generated:", output_filename, "\n")
  }
}

# Print completion summary
total_comparisons <- length(assay_names) * (length(assay_names) - 1) / 2
cat("\n=== Analysis Complete ===\n")
cat("Total pairwise comparisons generated:", total_comparisons, "\n")
cat("Output directory:", output_directory, "\n")
cat("File format: PDF (vector graphics)\n")
cat("Axis limits: X(", paste(x_axis_limits, collapse = ", "), 
    "), Y(", paste(y_axis_limits, collapse = ", "), ")\n")
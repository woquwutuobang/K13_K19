# Library Fitness Comparison Analysis
# This script compares fitness data between two different libraries (e.g., synthetic vs nicking)

library(wlab.block)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

#' Compare fitness data between two libraries for single mutations
#'
#' @param lib1_block1 Library 1 block 1 data file path
#' @param lib1_block2 Library 1 block 2 data file path
#' @param lib1_block3 Library 1 block 3 data file path
#' @param lib2_block1 Library 2 block 1 data file path
#' @param lib2_block2 Library 2 block 2 data file path
#' @param lib2_block3 Library 2 block 3 data file path
#' @param wt_aa Wild-type amino acid sequence
#' @param output_file Output file path for saving the plot (optional)
#' @param x_lab X-axis label
#' @param y_lab Y-axis label
#' @param main_title Main plot title
#' @param point_alpha Point transparency (default: 0.3)
#' @param plot_width Plot width in inches (default: 16)
#' @param plot_height Plot height in inches (default: 5)
#' @return List containing merged data and plot object

compare_fitness_libraries_singlemut <- function(
    lib1_block1, lib1_block2, lib1_block3,
    lib2_block1, lib2_block2, lib2_block3,
    wt_aa,
    output_file = NULL,
    x_lab = "Library 1 fitness",
    y_lab = "Library 2 fitness", 
    main_title = "Comparison of fitness data between two libraries",
    point_alpha = 0.3,
    plot_width = 16,
    plot_height = 5
) {
  
  # Internal function: process single library data
  process_library_data <- function(block1, block2, block3, wt_aa, suffix) {
    nor_fit <- nor_fitness(block1 = block1, block2 = block2, block3 = block3)
    nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- pos_id(nor_fit_single, wt_aa)
    
    fitness_data <- nor_fit_single[, c(1, 40, 41, 46, 48, 50, 52)]
    colnames(fitness_data) <- c("block", 
                                paste0("fitness", suffix),
                                paste0("fitness_sigma", suffix),
                                "Pos", "wtcodon", "codon", "mt")
    return(fitness_data)
  }
  
  # Process data from both libraries
  cat("Processing library 1 data...\n")
  fitness_data_1 <- process_library_data(lib1_block1, lib1_block2, lib1_block3, wt_aa, "1")
  
  cat("Processing library 2 data...\n")
  fitness_data_2 <- process_library_data(lib2_block1, lib2_block2, lib2_block3, wt_aa, "2")
  
  # Merge data
  data <- merge(fitness_data_1, fitness_data_2, by = c("block", "Pos", "wtcodon", "codon", "mt"), all = FALSE)
  setDT(data)
  
  cat(paste("Total variants after merging:", nrow(data), "\n"))
  
  # Create correlation plot function
  create_cor_plot <- function(data, title, x_label, y_label) {
    # Remove missing values
    complete_cases <- complete.cases(data$fitness1, data$fitness2)
    data_complete <- data[complete_cases, ]
    
    if (nrow(data_complete) < 2) {
      warning(paste("Insufficient complete cases for plot:", title))
      return(ggplot() + 
               labs(title = title, x = x_label, y = y_label) +
               annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
               theme_minimal())
    }
    
    # Calculate Pearson correlation and p-value
    cor_test <- cor.test(data_complete$fitness1, data_complete$fitness2, 
                         method = "pearson", use = "complete.obs")
    r_value <- round(cor_test$estimate, 3)
    p_value <- round(cor_test$p.value, 4)
    
    # Create plot
    p <- ggplot(data_complete, aes(x = fitness1, y = fitness2)) +
      geom_point(color = "#48B3AF", alpha = point_alpha, size = 1.5) +
      geom_smooth(method = "lm", color = "#9A3F3F", se = TRUE, 
                  fill = "gray70", alpha = 0.3) +
      labs(
        title = title,
        x = x_label,
        y = y_label
      ) +
      annotate("text", 
               x = min(data_complete$fitness1, na.rm = TRUE) + 
                 0.1 * diff(range(data_complete$fitness1, na.rm = TRUE)),
               y = max(data_complete$fitness2, na.rm = TRUE) - 
                 0.1 * diff(range(data_complete$fitness2, na.rm = TRUE)),
               label = paste0("r = ", r_value, "\np = ", 
                              ifelse(p_value < 0.0001, "< 0.0001", p_value)),
               hjust = 0, vjust = 1, size = 3.5,
               color = "black", fontface = "bold") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    
    return(p)
  }
  
  # Create overall plot
  p_total <- create_cor_plot(data, "Overall", x_lab, y_lab)
  
  # Create block subplots
  block_plots <- list()
  blocks <- unique(data$block)
  
  for (i in seq_along(blocks)) {
    block_data <- data[block == blocks[i]]
    p_block <- create_cor_plot(block_data, paste("Region:", blocks[i]), x_lab, y_lab)
    block_plots[[i]] <- p_block
  }
  
  # Combine plots
  if (length(blocks) == 3) {
    combined_plot <- p_total + block_plots[[1]] + block_plots[[2]] + block_plots[[3]] +
      plot_annotation(
        title = main_title,
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      ) +
      plot_layout(ncol = 4)
  } else {
    # If number of blocks is not 3, use automatic layout
    combined_plot <- p_total + wrap_plots(block_plots) +
      plot_annotation(
        title = main_title,
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      ) +
      plot_layout(ncol = min(length(blocks) + 1, 4))
  }
  
  # Display plot
  print(combined_plot)
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(filename = output_file,
           plot = combined_plot,
           device = cairo_pdf,
           width = plot_width,
           height = plot_height,
           units = "in",
           dpi = 300)
    cat(paste("Plot saved to:", output_file, "\n"))
  }
  
  # Return data and plot object
  return(list(data = data, plot = combined_plot))
}

# =====================
# =====================

# Define wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

### Abundance comparison
result_abundance <- compare_fitness_libraries_singlemut(
  # Library 1 (nicking library)
  lib1_block1 = "path/to/abundance_block1_nicking.RData",
  lib1_block2 = "path/to/abundance_block2_nicking.RData", 
  lib1_block3 = "path/to/abundance_block3_nicking.RData",
  
  # Library 2 (synthetic library)
  lib2_block1 = "path/to/abundance_block1_synthetic.RData",
  lib2_block2 = "path/to/abundance_block2_synthetic.RData",
  lib2_block3 = "path/to/abundance_block3_synthetic.RData",
  
  wt_aa = wt_aa,
  output_file = "path/to/abundance_comparison.pdf",
  x_lab = "Abundance nicking library fitness",
  y_lab = "Abundance synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library",
  point_alpha = 0.3,
  plot_width = 16,
  plot_height = 5
)

### RAF1 binding comparison
result_raf1 <- compare_fitness_libraries_singlemut(
  # Library 1 (nicking library)
  lib1_block1 = "path/to/raf1_block1_nicking.RData",
  lib1_block2 = "path/to/raf1_block2_nicking.RData", 
  lib1_block3 = "path/to/raf1_block3_nicking.RData",
  
  # Library 2 (synthetic library)
  lib2_block1 = "path/to/raf1_block1_synthetic.RData",
  lib2_block2 = "path/to/raf1_block2_synthetic.RData",
  lib2_block3 = "path/to/raf1_block3_synthetic.RData",
  
  wt_aa = wt_aa,
  output_file = "path/to/raf1_comparison.pdf",
  x_lab = "RAF1 nicking library fitness",
  y_lab = "RAF1 synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library",
  point_alpha = 0.3,
  plot_width = 16,
  plot_height = 5
)
# Library Fitness Comparison Analysis
# Compare fitness data between two libraries (e.g., synthetic vs nicking)

library(wlab.block)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

#' Compare Fitness Data Between Two Libraries
#'
#' This function compares fitness data between two different libraries (e.g., synthetic vs nicking libraries)
#' by processing normalized fitness data from multiple blocks and creating correlation plots.
#' The function handles single mutation fitness data and generates comprehensive visualizations
#' showing correlations between libraries both overall and by genomic regions/blocks.
#'
#' @param lib1_block1 character. Path to the RData file containing fitness data for library 1, block 1.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param lib1_block2 character. Path to the RData file containing fitness data for library 1, block 2.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param lib1_block3 character. Path to the RData file containing fitness data for library 1, block 3.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param lib2_block1 character. Path to the RData file containing fitness data for library 2, block 1.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param lib2_block2 character. Path to the RData file containing fitness data for library 2, block 2.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param lib2_block3 character. Path to the RData file containing fitness data for library 2, block 3.
#'   Should contain fitness replicates data that can be processed by \code{nor_fitness()} function.
#' @param wt_aa character. Wild-type amino acid sequence of the protein being analyzed.
#'   Used for position identification and mutation mapping. Should be a single string containing
#'   the complete amino acid sequence in single-letter code.
#' @param output_file character, optional. Path where the output plot should be saved as PDF.
#'   If NULL (default), the plot will only be displayed without saving. Should include
#'   the full path and filename with .pdf extension.
#' @param x_lab character. Label for the x-axis of the correlation plots. Default: "Library 1 fitness".
#'   Should describe what library 1 represents (e.g., "Abundance nicking library fitness").
#' @param y_lab character. Label for the y-axis of the correlation plots. Default: "Library 2 fitness".
#'   Should describe what library 2 represents (e.g., "Abundance synthetic library fitness").
#' @param main_title character. Main title for the combined plot. Default: "Comparison of fitness data between two libraries".
#'   Should provide a clear description of what is being compared.
#' @param point_alpha numeric. Transparency level for scatter plot points. Range: 0 (completely transparent) to 1 (completely opaque).
#'   Default: 0.3. Lower values make overlapping points more visible in dense datasets.
#' @param plot_width numeric. Width of the output plot in inches. Default: 16.
#'   Used when saving the plot to file. Should be adjusted based on the number of subplots.
#' @param plot_height numeric. Height of the output plot in inches. Default: 5.
#'   Used when saving the plot to file. Should be adjusted based on the desired aspect ratio.
#'
#' @return A list containing:
#'   \item{data}{data.table. Merged fitness data from both libraries containing columns:
#'     \itemize{
#'       \item block: Genomic region/block identifier
#'       \item fitness1: Normalized fitness values from library 1
#'       \item fitness_sigma1: Standard error of fitness values from library 1
#'       \item fitness2: Normalized fitness values from library 2
#'       \item fitness_sigma2: Standard error of fitness values from library 2
#'       \item Pos: Amino acid position in the protein sequence
#'       \item wtcodon: Wild-type codon at the position
#'       \item codon: Mutant codon at the position
#'       \item mt: Mutation type (e.g., "A12V")
#'     }
#'   }
#'   \item{plot}{ggplot object. Combined correlation plot showing:
#'     \itemize{
#'       \item Overall correlation between libraries
#'       \item Block-specific correlations (one subplot per genomic region)
#'       \item Pearson correlation coefficients and p-values
#'       \item Linear regression lines with confidence intervals
#'     }
#'   }
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Loads and processes fitness data from both libraries using \code{nor_fitness()} and \code{nor_fitness_single_mut()}
#'   \item Merges data from both libraries based on position and mutation information
#'   \item Creates correlation plots showing fitness relationships
#'   \item Calculates Pearson correlation coefficients and statistical significance
#'   \item Generates both overall and block-specific visualizations
#'   \item Optionally saves the plot as a high-resolution PDF
#' }
#'
#' The function requires the \code{wlab.block} package and expects RData files containing
#' fitness replicate data that can be processed by the normalization functions.
#'
#' @examples
#' # Basic usage comparing abundance libraries
#' wt_sequence <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
#' 
#' result <- compare_fitness_libraries_singlemut(
#'   lib1_block1 = "path/to/nicking_block1.RData",
#'   lib1_block2 = "path/to/nicking_block2.RData", 
#'   lib1_block3 = "path/to/nicking_block3.RData",
#'   lib2_block1 = "path/to/synthetic_block1.RData",
#'   lib2_block2 = "path/to/synthetic_block2.RData",
#'   lib2_block3 = "path/to/synthetic_block3.RData",
#'   wt_aa = wt_sequence,
#'   output_file = "library_comparison.pdf",
#'   x_lab = "Nicking library fitness",
#'   y_lab = "Synthetic library fitness",
#'   main_title = "Fitness comparison: Nicking vs Synthetic libraries"
#' )
#' 
#' # Access the merged data
#' merged_data <- result$data
#' print(paste("Total variants compared:", nrow(merged_data)))
#' 
#' # Access the plot
#' comparison_plot <- result$plot
#' print(comparison_plot)
#' 
#' # Example with custom plot parameters
#' result_custom <- compare_fitness_libraries_singlemut(
#'   lib1_block1 = "path/to/lib1_block1.RData",
#'   lib1_block2 = "path/to/lib1_block2.RData",
#'   lib1_block3 = "path/to/lib1_block3.RData",
#'   lib2_block1 = "path/to/lib2_block1.RData",
#'   lib2_block2 = "path/to/lib2_block2.RData",
#'   lib2_block3 = "path/to/lib2_block3.RData",
#'   wt_aa = wt_sequence,
#'   point_alpha = 0.5,
#'   plot_width = 14,
#'   plot_height = 4,
#'   x_lab = "RAF1 binding fitness (Library A)",
#'   y_lab = "RAF1 binding fitness (Library B)",
#'   main_title = "RAF1 binding fitness correlation analysis"
#' )
#'
#' @seealso \code{\link{nor_fitness}}, \code{\link{nor_fitness_single_mut}}, \code{\link{pos_id}}
#' @import wlab.block data.table ggplot2 dplyr patchwork
#' @export

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
  
  # Process both libraries
  cat("Processing library 1 data...\n")
  fitness_data_1 <- process_library_data(lib1_block1, lib1_block2, lib1_block3, wt_aa, "1")
  
  cat("Processing library 2 data...\n")
  fitness_data_2 <- process_library_data(lib2_block1, lib2_block2, lib2_block3, wt_aa, "2")
  
  # Merge data
  data <- merge(fitness_data_1, fitness_data_2, by = c("block", "Pos", "wtcodon", "codon", "mt"), all = FALSE)
  setDT(data)
  
  cat(paste("Total variants after merging:", nrow(data), "\n"))
  
  # ---- Correlation plot ----
  create_cor_plot <- function(data, title, x_label, y_label) {
    complete_cases <- complete.cases(data$fitness1, data$fitness2)
    data_complete <- data[complete_cases, ]
    
    if (nrow(data_complete) < 2) {
      warning(paste("Insufficient complete cases for plot:", title))
      return(ggplot() + 
               labs(title = title, x = x_label, y = y_label) +
               annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
               theme_minimal(base_size = 8))
    }
    
    # Pearson correlation
    cor_test <- cor.test(data_complete$fitness1, data_complete$fitness2, 
                         method = "pearson", use = "complete.obs")
    r_value <- round(cor_test$estimate, 3)
    p_value <- round(cor_test$p.value, 4)
    
    # Plot
    p <- ggplot(data_complete, aes(x = fitness1, y = fitness2)) +
      geom_point(color = "#75C2F6", alpha = point_alpha, size = 1.5) +
      geom_smooth(method = "lm", color = "#FF6A56", se = TRUE, 
                  fill = "gray70", alpha = 0.3) +
      coord_cartesian(xlim = c(-1.5, 1.0), ylim = c(-1.5, 0.5)) +   # Fixed coordinate range
      labs(
        title = title,
        x = x_label,
        y = y_label
      ) +
      annotate("text", 
               x = -1.4, y = 0.45,
               label = paste0("r = ", r_value, "\np = ", 
                              ifelse(p_value < 0.0001, "< 0.0001", p_value)),
               hjust = 0, vjust = 1, size = 2.5,
               color = "black") +
      theme_classic(base_size = 8) + 
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),  
        axis.title = element_text(size = 8),
        legend.position = "none",
        plot.margin = margin(3, 3, 3, 3),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    
    return(p)
  }
  
  # Overall plot
  p_total <- create_cor_plot(data, "Overall", x_lab, y_lab)
  
  # Subplots by block
  block_plots <- list()
  blocks <- unique(data$block)
  
  for (i in seq_along(blocks)) {
    block_data <- data[block == blocks[i]]
    p_block <- create_cor_plot(block_data, paste("Region:", blocks[i]), x_lab, y_lab)
    block_plots[[i]] <- p_block
  }
  
  # Combine
  if (length(blocks) == 3) {
    combined_plot <- p_total + block_plots[[1]] + block_plots[[2]] + block_plots[[3]] +
      plot_annotation(
        title = main_title,
        theme = theme(plot.title = element_text(hjust = 0.5, size = 8))
      ) +
      plot_layout(ncol = 4)
  } else {
    combined_plot <- p_total + wrap_plots(block_plots) +
      plot_annotation(
        title = main_title,
        theme = theme(plot.title = element_text(hjust = 0.5, size = 8))
      ) +
      plot_layout(ncol = min(length(blocks) + 1, 4))
  }
  
  # Show
  print(combined_plot)
  
  # Save
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
  
  return(list(data = data, plot = combined_plot))
}



#Execute function
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"


###Abundance
result <- compare_fitness_libraries_singlemut(
  #nicking library data
  lib1_block1 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
  
  #synthetic library data
  lib2_block1 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251010/Comparison_of_fitness_data_Abundance.pdf",
  x_lab = "Abundance nicking library fitness",
  y_lab = "Abundance synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library",
  point_alpha = 0.3,
  plot_width = 12,
  plot_height = 3
)





###RAF1
result <- compare_fitness_libraries_singlemut(
  #nicking library data
  lib1_block1 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
  
  #synthetic library data
  lib2_block1 = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251010/Comparison_of_fitness_data_RAF1.pdf",
  x_lab = "RAF1 nicking library fitness",
  y_lab = "RAF1 synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library",
  point_alpha = 0.3,
  plot_width = 12,
  plot_height = 3
)





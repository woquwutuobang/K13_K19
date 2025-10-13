# Library Fitness Comparison Analysis
# Compare fitness data between two libraries (e.g., synthetic vs nicking)

library(wlab.block)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

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
      theme_classic(base_size = 8) +  # 全局字体 8pt
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),  # X轴旋转90度
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


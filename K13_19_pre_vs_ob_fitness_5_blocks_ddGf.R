# MoCHI Model Fitness Integration and Visualization
# This script integrates predicted and observed fitness data from 5 experimental blocks for model evaluation

library(data.table)
library(krasddpcams)

#' Integrate predicted and observed fitness data from 5 experimental blocks
#'
#' @param prediction Path to prediction data file
#' @param folding_ddG Path to folding free energy data file
#' @param block1_dimsum_df Path to block 1 fitness data
#' @param block2_dimsum_df Path to block 2 fitness data
#' @param block3_dimsum_df Path to block 3 fitness data
#' @param block4_dimsum_df Path to block 4 fitness data
#' @param block5_dimsum_df Path to block 5 fitness data
#' @param wt_aa_input Wild-type amino acid sequence
#' @return Integrated data table with normalized fitness values

krasddpcams__merge_ddGf_fitness_5blocks <- function(prediction = prediction, 
                                           folding_ddG = folding_ddG, 
                                           block1_dimsum_df = block1_dimsum_df, 
                                           block2_dimsum_df = block2_dimsum_df, 
                                           block3_dimsum_df = block3_dimsum_df,
                                           block4_dimsum_df = block4_dimsum_df,
                                           block5_dimsum_df = block5_dimsum_df,
                                           wt_aa_input = wt_aa_input) {
  
  # Read prediction and folding free energy data
  pre <- fread(prediction)
  folding_ddG <- fread(folding_ddG)
  pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
  
  # Load all 5 block data
  load(block1_dimsum_df)
  block1 <- as.data.table(all_variants)
  load(block2_dimsum_df)
  block2 <- as.data.table(all_variants)
  load(block3_dimsum_df)
  block3 <- as.data.table(all_variants)
  load(block4_dimsum_df)
  block4 <- as.data.table(all_variants)
  load(block5_dimsum_df)
  block5 <- as.data.table(all_variants)
  
  # Merge all block data
  data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3,
                           block4 = block4, block5 = block5, 
                           idcol = "block", fill = TRUE)
  
  # Calculate intermediate values for weighted fitness
  data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness / (data_before_nor$sigma)^2
  data_before_nor$one_over_fitness_sigmasquared <- 1 / (data_before_nor$sigma)^2
  
  # Calculate stop codon fitness for each block (weighted average)
  calculate_stop_fitness <- function(block_name) {
    dead_fitness <- data_before_nor[STOP == TRUE & block == block_name, ]
    stop_fitness <- sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(stop_fitness)
  }
  
  stop1_fitness <- calculate_stop_fitness("block1")
  stop2_fitness <- calculate_stop_fitness("block2")
  stop3_fitness <- calculate_stop_fitness("block3")
  stop4_fitness <- calculate_stop_fitness("block4")
  stop5_fitness <- calculate_stop_fitness("block5")
  
  # Calculate wild-type fitness for each block (weighted average)
  calculate_wt_fitness <- function(block_name) {
    wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
    wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(wt_fitness_block$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(wt_fitness)
  }
  
  wt1_fitness <- calculate_wt_fitness("block1")
  wt2_fitness <- calculate_wt_fitness("block2")
  wt3_fitness <- calculate_wt_fitness("block3")
  wt4_fitness <- calculate_wt_fitness("block4")
  wt5_fitness <- calculate_wt_fitness("block5")
  
  # Create scaling data for normalization
  scaling_data_fitness <- data.frame(
    block1 = c(stop1_fitness, wt1_fitness),
    block2 = c(stop2_fitness, wt2_fitness),
    block3 = c(stop3_fitness, wt3_fitness),
    block4 = c(stop4_fitness, wt4_fitness),
    block5 = c(stop5_fitness, wt5_fitness)
  )
  
  # Calculate scaling parameters relative to block1
  calculate_scaling_params <- function(target_block) {
    lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
    return(list(
      slope = lm_model$coefficients[[2]],
      intercept = lm_model$coefficients[[1]]
    ))
  }
  
  params_block2 <- calculate_scaling_params("block2")
  params_block3 <- calculate_scaling_params("block3")
  params_block4 <- calculate_scaling_params("block4")
  params_block5 <- calculate_scaling_params("block5")
  
  d2 <- params_block2$slope
  e2 <- params_block2$intercept
  d3 <- params_block3$slope
  e3 <- params_block3$intercept
  d4 <- params_block4$slope
  e4 <- params_block4$intercept
  d5 <- params_block5$slope
  e5 <- params_block5$intercept
  
  # Process prediction data
  pre_nor <- pre_pos
  
  # Extract predicted fitness (adjust indices based on actual data structure)
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # Extract additive traits (adjust indices based on actual data structure)
  extract_additive_trait0 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2 - 1])
  }
  extract_additive_trait1 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2])
  }
  
  pre_nor$additive_trait0 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait0)
  pre_nor$additive_trait0 <- as.numeric(pre_nor$additive_trait0)
  pre_nor$additive_trait1 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait1)
  pre_nor$additive_trait1 <- as.numeric(pre_nor$additive_trait1)
  pre_nor[, additive_trait := additive_trait0 + additive_trait1]
  
  # Apply normalization transformations (5 blocks)
  pre_nor[phenotype == 1, pre_nor_mean_fitness := mean]
  pre_nor[phenotype == 2, pre_nor_mean_fitness := mean * d2 + e2]
  pre_nor[phenotype == 3, pre_nor_mean_fitness := mean * d3 + e3]
  pre_nor[phenotype == 4, pre_nor_mean_fitness := mean * d4 + e4]
  pre_nor[phenotype == 5, pre_nor_mean_fitness := mean * d5 + e5]
  
  pre_nor[phenotype == 1, pre_nor_fitness_sigma := std]
  pre_nor[phenotype == 2, pre_nor_fitness_sigma := std * d2]
  pre_nor[phenotype == 3, pre_nor_fitness_sigma := std * d3]
  pre_nor[phenotype == 4, pre_nor_fitness_sigma := std * d4]
  pre_nor[phenotype == 5, pre_nor_fitness_sigma := std * d5]
  
  pre_nor[phenotype == 1, ob_nor_fitness := fitness]
  pre_nor[phenotype == 2, ob_nor_fitness := fitness * d2 + e2]
  pre_nor[phenotype == 3, ob_nor_fitness := fitness * d3 + e3]
  pre_nor[phenotype == 4, ob_nor_fitness := fitness * d4 + e4]
  pre_nor[phenotype == 5, ob_nor_fitness := fitness * d5 + e5]
  
  pre_nor[phenotype == 1, ob_nor_fitness_sigma := sigma]
  pre_nor[phenotype == 2, ob_nor_fitness_sigma := sigma * d2]
  pre_nor[phenotype == 3, ob_nor_fitness_sigma := sigma * d3]
  pre_nor[phenotype == 4, ob_nor_fitness_sigma := sigma * d4]
  pre_nor[phenotype == 5, ob_nor_fitness_sigma := sigma * d5]
  
  pre_nor[phenotype == 1, pre_nor_fitness := predicted_fitness]
  pre_nor[phenotype == 2, pre_nor_fitness := predicted_fitness * d2 + e2]
  pre_nor[phenotype == 3, pre_nor_fitness := predicted_fitness * d3 + e3]
  pre_nor[phenotype == 4, pre_nor_fitness := predicted_fitness * d4 + e4]
  pre_nor[phenotype == 5, pre_nor_fitness := predicted_fitness * d5 + e5]
  
  return(pre_nor)
}



krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock <- function(pre_nor = pre_nor, phenotypen = phenotypen, rotate_x_axis = TRUE) 
{
  pre_nor
  lm_mochi <- lm(pre_nor_fitness ~ ob_nor_fitness, pre_nor[phenotype == phenotypen, ])
  
  p <- ggplot2::ggplot() + 
    ggplot2::stat_binhex(data = pre_nor[phenotype == phenotypen, ], 
                         ggplot2::aes(x = ob_nor_fitness, y = pre_nor_fitness), 
                         bins = 50, size = 0, color = "black") + 
    ggplot2::scale_fill_gradient(low = "white", 
                                 high = "black", 
                                 trans = "log10", 
                                 guide = ggplot2::guide_colorbar(barwidth = 0.5, 
                                                                 barheight = 1.5)) + 
    ggplot2::geom_hline(yintercept = 0) + 
    ggplot2::geom_vline(xintercept = 0) + 
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    ggplot2::annotate("text", 
                      x = -0.8, y = 0.3, 
                      label = paste0("R² = ", round(summary(lm_mochi)$r.squared, 2)), 
                      size = 8 * 0.35) + 
    ggplot2::theme_classic() + 
    ggplot2::xlab("Observed fitness") + 
    ggplot2::ylab("Predicted fitness") + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 8), 
      axis.text = ggplot2::element_text(size = 8), 
      legend.text = ggplot2::element_text(size = 8), 
      legend.key.size = ggplot2::unit(1, "cm"), 
      plot.title = ggplot2::element_text(size = 8)
    ) + 
    ggplot2::coord_fixed()
  
  # Determines whether to rotate the X axis according to the parameters
  if (rotate_x_axis) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  
  return(p)
}


# =====================
# Model Evaluation Plots
# =====================

# Define wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Integrate abundance fitness data
abundance_fitness_data <- krasddpcams__merge_ddGf_fitness_5blocks(
  prediction="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/predictions/predicted_phenotypes_all.txt",
  folding_ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt",
  block1_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  block2_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
  block3_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block4_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
  block5_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  wt_aa_input = wt_aa
)

# Generate evaluation plots for each block
krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = abundance_fitness_data, phenotypen = 1, rotate_x_axis = TRUE)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251019/abundance_block1_evaluation_fitness_pre_vs_ob.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = abundance_fitness_data, phenotypen = 2, rotate_x_axis = TRUE)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251019/abundance_block2_1_evaluation_fitness_pre_vs_ob.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = abundance_fitness_data, phenotypen = 3, rotate_x_axis = TRUE)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251019/abundance_block2_2_evaluation_fitness_pre_vs_ob.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = abundance_fitness_data, phenotypen = 4, rotate_x_axis = TRUE)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251019/abundance_block3_1_evaluation_fitness_pre_vs_ob.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = abundance_fitness_data, phenotypen = 5, rotate_x_axis = TRUE)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251019/abundance_block3_2_evaluation_fitness_pre_vs_ob.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

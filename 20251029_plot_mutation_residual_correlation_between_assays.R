library(ggplot2)
library(data.table)
library(ggpubr)
library(dplyr)

# ===============================
# 主函数：两种 assay 的残差相关性分析与绘图
# ===============================
plot_mutation_residual_correlation_between_assays <- function(
    input1, assay1,
    input2, assay2,
    anno_file,
    output_pdf,
    interface_sites1 = NULL,
    interface_sites2 = NULL
) {
  
  # ========== Step 1: 计算两个 assay 的残差 ==========
  cat("=== 计算", assay1, "的残差 ===\n")
  result1 <- calculate_distance_residuals(input1, assay1, anno_file)
  
  cat("=== 计算", assay2, "的残差 ===\n")
  result2 <- calculate_distance_residuals(input2, assay2, anno_file)
  
  # 检查结果是否为空
  if (is.null(result1) || is.null(result2)) {
    cat("❌ 一个或多个assay的残差计算失败\n")
    return(NULL)
  }
  
  # ========== Step 2: 处理第一个 assay 数据 ==========
  ddG1 <- fread(input1)
  ddG1 <- ddG1[, Pos_real := Pos + 1]
  ddG1 <- ddG1[, c(1, 20, 23)]
  colnames(ddG1)[2] <- paste0(colnames(ddG1)[2], "_", assay1)
  
  # 重命名result1的列
  setnames(result1, "original_y", paste0("mean_kcal/mol_", assay1))
  
  # 合并第一个assay数据
  merge1 <- merge(ddG1, result1, by = c("Pos_real", paste0("mean_kcal/mol_", assay1)), all = FALSE)
  
  # ========== Step 3: 处理第二个 assay 数据 ==========
  ddG2 <- fread(input2)
  ddG2 <- ddG2[, Pos_real := Pos + 1]
  ddG2 <- ddG2[, c(1, 20, 23)]
  colnames(ddG2)[2] <- paste0(colnames(ddG2)[2], "_", assay2)
  
  # 重命名result2的列
  setnames(result2, "original_y", paste0("mean_kcal/mol_", assay2))
  
  # 合并第二个assay数据
  merge2 <- merge(ddG2, result2, by = c("Pos_real", paste0("mean_kcal/mol_", assay2)), all = FALSE)
  
  # ========== Step 4: 重命名残差列并合并 ==========
  setnames(merge1, "residual", paste0("residual_", assay1))
  setnames(merge2, "residual", paste0("residual_", assay2))
  
  # 使用明确的列名进行合并
  data_plot <- merge(merge1, merge2, 
                     by = c("Pos_real", "id", "Pos"), 
                     all = FALSE)
  
  cat("合并后数据维度:", dim(data_plot), "\n")
  
  # ========== Step 5: 检查残差列是否存在 ==========
  residual_col1 <- paste0("residual_", assay1)
  residual_col2 <- paste0("residual_", assay2)
  
  if (!all(c(residual_col1, residual_col2) %in% names(data_plot))) {
    cat("❌ 错误: 残差列不存在！可用列名为:\n")
    print(names(data_plot))
    stop("请检查列名设置")
  }
  
  # ========== Step 6: 过滤结合界面位点 ==========
  if (!is.null(interface_sites1) && !is.null(interface_sites2)) {
    all_interface_sites <- unique(c(interface_sites1, interface_sites2))
    data_plot_filtered <- data_plot[!Pos_real %in% all_interface_sites]
    
    cat("过滤统计:\n")
    cat("过滤前数据点:", nrow(data_plot), "\n")
    cat("过滤后数据点:", nrow(data_plot_filtered), "\n")
    cat("排除的结合界面位点数量:", length(all_interface_sites), "\n")
  } else {
    data_plot_filtered <- data_plot
    cat("未提供结合界面位点，使用所有数据点:", nrow(data_plot_filtered), "\n")
  }
  
  # ========== Step 7: 准备相关性分析数据 ==========
  df <- data_plot_filtered %>%
    select(residual1 = all_of(residual_col1),
           residual2 = all_of(residual_col2)) %>%
    filter(!is.na(residual1) & !is.na(residual2))
  
  cat("去除NA后数据点:", nrow(df), "\n")
  
  if (nrow(df) == 0) {
    cat("❌ 没有有效数据点计算相关性\n")
    return(NULL)
  }
  
  # ========== Step 8: 计算相关系数 ==========
  cor_value <- cor(df$residual1, df$residual2, method = "pearson")
  cat("Pearson correlation =", round(cor_value, 3), "\n")
  
  # ========== Step 9: 绘图 ==========
  p <- ggplot(df, aes(x = residual1, y = residual2)) +
    geom_point(alpha = 0.6, color = "#75C2F6", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", fill = "#FFB0A5", alpha = 0.3) +
    stat_cor(method = "pearson", 
             label.x = min(df$residual1, na.rm = TRUE), 
             label.y = max(df$residual2, na.rm = TRUE) - 0.5, 
             size = 2.5,
             label.sep = "\n") +
    labs(
      x = paste0("Residual ", assay1, " (kcal/mol)"),
      y = paste0("Residual ", assay2, " (kcal/mol)"),
      title = paste0("Correlation of Residuals between ", assay1, " and ", assay2,
                     if (!is.null(interface_sites1) && !is.null(interface_sites2)) 
                       "\n(Excluding Binding Interface Sites)" else "")
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.border = element_rect(color = "black", linewidth = 0.8)
    ) +
    coord_cartesian(xlim = c(-2, 2.1), ylim = c(-1.7, 2.3))
  
  # ========== Step 10: 保存图表 ==========
  ggsave(output_pdf, p, device = cairo_pdf, width = 4, height = 4)
  
  # ========== Step 11: 返回结果 ==========
  return(list(
    plot = p,
    correlation = cor_value,
    data = df,
    merged_data = data_plot_filtered
  ))
}

# ===============================
# 使用示例
# ===============================

#### K13 vs RAF1
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

# 示例1: K13 vs RAF1
result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5a_mutation_residues_Correlation_K13_vs_RAF1.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$plot)
  cat("K13 vs RAF1 相关系数:", round(result$correlation, 3), "\n")
}




#### K27 vs RAF1
# 定义结合界面位点
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)


result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5a_mutation_residues_Correlation_K27_vs_RAF1.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$plot)
  cat("K27 vs RAF1 相关系数:", round(result$correlation, 3), "\n")
}



#### K13 vs K19
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay2 = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5a_mutation_residues_Correlation_K13_vs_K19.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = K19_Binding_interface_site
)

if (!is.null(result)) {
  print(result$plot)
  cat("K13 vs K19 相关系数:", round(result$correlation, 3), "\n")
}

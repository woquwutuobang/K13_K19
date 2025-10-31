library(data.table)
library(ggplot2)
library(dplyr)
library(krasddpcams)

# 给定的变构位点
K13_allosteric_site <- c(15, 16, 17, 35, 145, 10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)
K27_allosteric_site <- c(16, 17, 28, 57, 60, 117, 119, 10, 63, 76, 84, 90, 144, 151, 155, 156)
K19_allosteric_site <- c(15, 16, 17, 145,8, 10, 19, 21, 55, 77, 78, 79, 82, 93, 151, 159, 163)
# 计算位点级别中位数的函数
calculate_site_median <- function(input, assay_name, wt_aa) {
  # Read and process ddG data
  ddG <- fread(input)
  ddG[, `:=`(Pos_real = Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon = substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon = substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt = paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]]))
  )
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate median ddG for each position
  site_data <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), .(
    median_ddG = median(`mean_kcal/mol`, na.rm = TRUE),
    n_mutations = sum(!is.na(`mean_kcal/mol`))
  ), by = .(Pos_real)]
  
  # Add assay identifier
  site_data[, assay := assay_name]
  
  return(site_data)
}

# 绘制位点级别相关性散点图的函数
plot_site_level_correlation <- function(input1, input2, assay1, assay2, wt_aa, 
                                        point_size = 3, alpha = 0.7, base_size = 10,
                                        xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3)) {
  
  # 获取两个assay的位点级别数据
  site_data1 <- calculate_site_median(input1, assay1, wt_aa)
  site_data2 <- calculate_site_median(input2, assay2, wt_aa)
  
  # 选择需要的列并重命名
  data1_clean <- site_data1[, .(Pos_real, median_ddG, n_mutations)]
  setnames(data1_clean, "median_ddG", paste0("median_", assay1))
  setnames(data1_clean, "n_mutations", paste0("n_mutations_", assay1))
  
  data2_clean <- site_data2[, .(Pos_real, median_ddG, n_mutations)]
  setnames(data2_clean, "median_ddG", paste0("median_", assay2))
  setnames(data2_clean, "n_mutations", paste0("n_mutations_", assay2))
  
  # 合并两个assay的数据
  merged_data <- merge(data1_clean, data2_clean, by = "Pos_real", all = TRUE)
  
  # 根据给定的变构位点定义颜色分组
  merged_data[, color_group := "Other"]
  
  # 获取对应assay的变构位点
  allosteric_sites1 <- get(paste0(assay1, "_allosteric_site"))
  allosteric_sites2 <- get(paste0(assay2, "_allosteric_site"))
  
  # 设置颜色分组
  merged_data[Pos_real %in% allosteric_sites1 & !Pos_real %in% allosteric_sites2, 
              color_group := paste0(assay1, " allosteric")]
  
  merged_data[Pos_real %in% allosteric_sites2 & !Pos_real %in% allosteric_sites1, 
              color_group := paste0(assay2, " allosteric")]
  
  merged_data[Pos_real %in% allosteric_sites1 & Pos_real %in% allosteric_sites2, 
              color_group := "Overlap allosteric"]
  
  # 动态设置颜色映射
  color_levels <- c("Other", 
                    paste0(assay1, " allosteric"), 
                    paste0(assay2, " allosteric"), 
                    "Overlap allosteric")
  
  # 动态颜色映射
  color_map <- setNames(
    c("grey", "#FF0066", "#A31300", "#C68EFD"),
    color_levels
  )
  
  # 设置因子水平以确保正确的绘图顺序
  merged_data[, color_group := factor(color_group, levels = color_levels)]
  
  # 按颜色分组排序数据，确保灰色(Other)在最底层
  merged_data <- merged_data[order(as.numeric(color_group))]
  
  # 创建散点图
  p <- ggplot(merged_data, aes(x = .data[[paste0("median_", assay1)]], 
                               y = .data[[paste0("median_", assay2)]], 
                               color = color_group)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = color_map, name = "Allosteric Sites",
                       breaks = color_levels[color_levels != "Other"]) +  # 图例中不显示Other
    labs(x = bquote("Median "*Delta*Delta*"G ("*.(assay1)*") (kcal/mol)"), 
         y = bquote("Median "*Delta*Delta*"G ("*.(assay2)*") (kcal/mol)")) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = base_size - 2),
      axis.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1),
      legend.key.size = unit(0.5, "cm")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  # 打印统计信息
  cat("=== Site Level Correlation Summary ===\n")
  cat("Total sites:", nrow(merged_data), "\n")
  cat(paste0(assay1, " allosteric sites: ", sum(merged_data$color_group == paste0(assay1, " allosteric")), "\n"))
  cat(paste0(assay2, " allosteric sites: ", sum(merged_data$color_group == paste0(assay2, " allosteric")), "\n"))
  cat("Overlap allosteric sites:", sum(merged_data$color_group == "Overlap allosteric"), "\n")
  cat("Other sites:", sum(merged_data$color_group == "Other"), "\n")
  
  # 计算相关系数
  cor_value <- cor(merged_data[[paste0("median_", assay1)]], 
                   merged_data[[paste0("median_", assay2)]], 
                   use = "complete.obs")
  cat("Correlation coefficient:", round(cor_value, 3), "\n\n")
  
  return(list(plot = p, data = merged_data))
}

# 使用示例
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 比较K13和RAF1
result_k13_raf1 <- plot_site_level_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K13",
  assay2 = "RAF1",
  wt_aa = wt_aa
)

# 显示图形
print(result_k13_raf1$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251031/K13_vs_RAF1_site_level_correlation.pdf", 
       plot = result_k13_raf1$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

# 比较K13和K19
result_k13_k19 <- plot_site_level_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K13",
  assay2 = "K19",
  wt_aa = wt_aa
)

# 显示图形
print(result_k13_k19$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251031/K13_vs_K19_site_level_correlation.pdf", 
       plot = result_k13_k19$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

# 比较K27和RAF1
result_k27_raf1 <- plot_site_level_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K27",
  assay2 = "RAF1",
  wt_aa = wt_aa
)

# 显示图形
print(result_k27_raf1$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251031/K27_vs_RAF1_site_level_correlation.pdf", 
       plot = result_k27_raf1$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

# 比较K13和K27
result_k13_k27 <- plot_site_level_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K13",
  assay2 = "K27",
  wt_aa = wt_aa
)

# 显示图形
print(result_k13_k27$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251027/K13_vs_K27_site_level_correlation.pdf", 
       plot = result_k13_k27$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

# 比较K19和RAF1
result_k19_raf1 <- plot_site_level_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K19",
  assay2 = "RAF1",
  wt_aa = wt_aa
)

# 显示图形
print(result_k19_raf1$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251027/K19_vs_RAF1_site_level_correlation.pdf", 
       plot = result_k19_raf1$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)
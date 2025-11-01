library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)
library(dplyr)

# ===============================
# 参数设置
# ===============================
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 定义变构位点信息
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146,10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)
RALGDS_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 117, 146,10, 20, 58, 59, 85)
PI3KCG_allosteric_site <- c(16, 17, 18, 28, 32, 34, 35, 57, 60, 117, 146,20, 55, 58, 59, 68, 85)
SOS1_allosteric_site <- c(16, 17, 18, 28, 32, 35, 57, 117, 146,10, 40, 54, 55, 58, 63, 68, 69, 85, 144, 148)
K55_allosteric_site <- c(15, 16, 17, 28, 32, 35, 57, 60, 117,10, 20, 58, 59, 68, 69, 71, 72, 78, 79, 85)
K27_allosteric_site <- c(16, 17, 28, 57, 60, 117, 119,10, 63, 76, 84, 90, 144, 151, 155, 156)
K13_allosteric_site <- c(15, 16, 17, 35, 145,10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
K19_allosteric_site <- c(15, 16, 17, 145,8, 10, 19, 21, 55, 77, 78, 79, 82, 93, 151, 159, 163)

# 定义结合界面位点信息
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)
PI3KCG_binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
GTP_Binding_pocket_site <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# 创建变构位点列表
allosteric_sites <- list(
  RAF1 = RAF1_allosteric_site,
  RALGDS = RALGDS_allosteric_site,
  PI3KCG = PI3KCG_allosteric_site,
  SOS1 = SOS1_allosteric_site,
  K55 = K55_allosteric_site,
  K27 = K27_allosteric_site,
  K13 = K13_allosteric_site,
  K19 = K19_allosteric_site
)

# 创建结合界面位点列表
binding_sites <- list(
  RAF1 = RAF1_Binding_interface_site,
  RALGDS = RALGDS_Binding_interface_site,
  PI3KCG = PI3KCG_binding_interface_site,
  SOS1 = SOS1_Binding_interface_site,
  K55 = K55_Binding_interface_site,
  K27 = K27_Binding_interface_site,
  K13 = K13_Binding_interface_site,
  K19 = K19_Binding_interface_site
)

# 文件路径
input_RAF1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt"
input_RALGDS <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt"
input_PI3KCG <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt"
input_SOS1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt"
input_K55 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt"
input_K27 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt"
input_K13 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt"
input_K19 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt"

# ===============================
# 数据预处理函数
# ===============================
read_ddG_data <- function(input, wt_aa) {
  ddG <- fread(input)
  ddG[, `:=`(Pos_real = Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon = substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon = substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt = paste0(wt_codon, Pos_real, mt_codon))]
  
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:(nchar(wt_aa) + 1), each = 20),
    mt_codon = rep(aa_list, times = nchar(wt_aa))
  )
  
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  return(ddG)
}

# 修改为使用median计算位点级别的ddG
calculate_median_ddG <- function(ddG) {
  # 计算每个位点的median ddG
  output <- ddG[Pos_real > 1, .(
    median_ddG = median(`mean_kcal/mol`, na.rm = TRUE),
    mad_ddG = mad(`mean_kcal/mol`, na.rm = TRUE),
    n_mutations = sum(!is.na(`mean_kcal/mol`))
  ), by = "Pos_real"]
  
  # 获取野生型氨基酸
  wt_codon_info <- ddG[Pos_real > 1, .(codon = first(wt_codon)), by = "Pos_real"]
  
  # 合并信息
  final_output <- merge(wt_codon_info, output, by = "Pos_real")
  final_output[, Pos := Pos_real]
  
  return(final_output)
}

# 基于median ddG和预定义的变构位点进行分类
classify_allosteric_sites <- function(median_ddG_data, assay_sele) {
  # 获取该assay的变构位点列表
  assay_allosteric_sites <- allosteric_sites[[assay_sele]]
  
  # 标记变构位点
  median_ddG_data[, is_allosteric := Pos_real %in% assay_allosteric_sites]
  
  return(median_ddG_data)
}

# ===============================
# 计算单个assay对的OR值（基于预定义的变构位点，只保留GTP pocket位点）
# ===============================
calculate_pair_or_pvalue_allosteric_gtp_only <- function(assay1, assay2) {
  cat("处理", assay1, "vs", assay2, "基于预定义变构位点（只保留GTP pocket）...\n")
  
  # 分别读取两个assay的数据
  input1 <- get(paste0("input_", assay1))
  input2 <- get(paste0("input_", assay2))
  
  # 处理assay1数据
  ddG_data1 <- read_ddG_data(input1, wt_aa)
  median_ddG1 <- calculate_median_ddG(ddG_data1)
  result1 <- classify_allosteric_sites(median_ddG1, assay1)
  
  # 处理assay2数据
  ddG_data2 <- read_ddG_data(input2, wt_aa)
  median_ddG2 <- calculate_median_ddG(ddG_data2)
  result2 <- classify_allosteric_sites(median_ddG2, assay2)
  
  # 合并这两个assay的数据（基于位点）
  data1 <- result1[, c("Pos_real", "is_allosteric")]
  setnames(data1, "is_allosteric", paste0("is_allosteric_", assay1))
  
  data2 <- result2[, c("Pos_real", "is_allosteric")]
  setnames(data2, "is_allosteric", paste0("is_allosteric_", assay2))
  
  pair_data <- merge(data1, data2, by = "Pos_real", all = FALSE)
  
  # 只保留GTP pocket位点
  filtered_data <- pair_data[Pos_real %in% GTP_Binding_pocket_site]
  
  # 检查是否有足够的数据
  if (nrow(filtered_data) == 0) {
    warning(paste("在GTP pocket中，", assay1, "vs", assay2, "没有有效位点"))
    return(NULL)
  }
  
  # 计算OR值
  is_effect1 <- filtered_data[[paste0("is_allosteric_", assay1)]] == TRUE
  is_effect2 <- filtered_data[[paste0("is_allosteric_", assay2)]] == TRUE
  
  contingency_table <- table(is_effect1, is_effect2)
  
  # 检查列联表是否有效
  if (any(dim(contingency_table) < 2) || any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) {
    warning(paste("在GTP pocket中，", assay1, "vs", assay2, "的列联表无效"))
    return(NULL)
  }
  
  tryCatch({
    fisher_test <- fisher.test(contingency_table)
    
    or_value <- fisher_test$estimate
    p_value <- fisher_test$p.value
    
    # 获取分类信息
    group1 <- ifelse(assay1 %in% c("K13", "K19"), "BI2", "BI1")
    group2 <- ifelse(assay2 %in% c("K13", "K19"), "BI2", "BI1")
    comparison_type <- paste0(group1, "_vs_", group2)
    
    return(data.table(
      assay1 = assay1,
      assay2 = assay2,
      group1 = group1,
      group2 = group2,
      comparison_type = comparison_type,
      OR = or_value,
      p_value = p_value,
      n_sites = nrow(filtered_data),
      n_effect1 = sum(is_effect1),
      n_effect2 = sum(is_effect2),
      n_both = sum(is_effect1 & is_effect2),
      sites_analyzed = length(GTP_Binding_pocket_site),
      status = "valid"
    ))
  }, error = function(e) {
    warning(paste("在GTP pocket中，", assay1, "vs", assay2, "的fisher检验出错:", e$message))
    return(NULL)
  })
}

# ===============================
# 主分析流程
# ===============================
cat("=== 开始分析所有assay组合（基于预定义变构位点，只分析GTP pocket） ===\n")

# 定义所有assay
assay_names <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# 计算所有组合的OR值
or_results <- data.table()
skipped_combinations <- data.table()

combinations <- combn(assay_names, 2, simplify = FALSE)

for (comb in combinations) {
  assay1 <- comb[1]
  assay2 <- comb[2]
  
  result <- calculate_pair_or_pvalue_allosteric_gtp_only(assay1, assay2)
  if (!is.null(result)) {
    or_results <- rbind(or_results, result)
  } else {
    skipped_combinations <- rbind(skipped_combinations, 
                                  data.table(assay1 = assay1, assay2 = assay2, reason = "无效的列联表"))
  }
}

# 添加显著性标记
if (nrow(or_results) > 0) {
  or_results[, significance := fcase(
    p_value < 0.001, "***",
    p_value < 0.01, "**",
    p_value < 0.05, "*",
    default = "NS"
  )]
}

# 打印结果
cat("\n=== 基于预定义变构位点的OR值结果（只分析GTP pocket） ===\n")
if (nrow(or_results) > 0) {
  print(or_results[order(comparison_type, OR)])
} else {
  cat("没有有效的OR结果\n")
}

# 打印跳过的组合
if (nrow(skipped_combinations) > 0) {
  cat("\n=== 跳过的组合 ===\n")
  print(skipped_combinations)
}

# ===============================
# 绘制柱状图 - 完整调试版本
# ===============================
cat("\n=== 绘制柱状图 - 调试信息 ===\n")

# 首先检查数据
cat("1. or_results 数据框行数:", nrow(or_results), "\n")
cat("2. or_results 列名:", names(or_results), "\n")

if (nrow(or_results) > 0) {
  # 详细检查OR值
  cat("3. OR值详情:\n")
  print(or_results[, .(assay1, assay2, OR, p_value, significance)])
  
  cat("4. OR值统计:\n")
  cat("   - 最小值:", min(or_results$OR, na.rm = TRUE), "\n")
  cat("   - 最大值:", max(or_results$OR, na.rm = TRUE), "\n")
  cat("   - 中位数:", median(or_results$OR, na.rm = TRUE), "\n")
  cat("   - NA数量:", sum(is.na(or_results$OR)), "\n")
  cat("   - Inf数量:", sum(is.infinite(or_results$OR)), "\n")
  
  # 检查是否有有效的OR值
  valid_or <- or_results[!is.na(OR) & is.finite(OR) & OR > 0]
  cat("5. 有效OR值数量:", nrow(valid_or), "\n")
  
  if (nrow(valid_or) == 0) {
    cat("!!! 警告: 没有有效的OR值可以绘图 !!!\n")
    # 创建一个显示错误信息的图
    p_error <- ggplot() +
      geom_text(aes(x = 1, y = 1, label = "No valid OR values to plot"), size = 8) +
      theme_void() +
      labs(title = "错误: 没有有效的OR值")
    print(p_error)
  } else {
    # 使用有效数据
    plot_data <- valid_or
    
    # 设置颜色
    comparison_colors <- c(
      "BI1_vs_BI1" = "#F4270C",
      "BI1_vs_BI2" = "#F4AD0C", 
      "BI2_vs_BI2" = "#1B38A6"
    )
    
    # 创建组合标签
    plot_data[, combination_label := paste0(assay1, " vs ", assay2)]
    
    # 排序
    plot_data[, comparison_type := factor(comparison_type, 
                                          levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2"))]
    plot_data[, within_group_order := frank(OR, ties.method = "first"), by = comparison_type]
    plot_data[, final_order := as.numeric(comparison_type) * 100 + within_group_order]
    plot_data <- plot_data[order(final_order)]
    plot_data[, combination_label := factor(combination_label, levels = unique(combination_label[order(final_order)]))]
    
    # 计算y轴范围
    max_or <- max(plot_data$OR, na.rm = TRUE)
    min_or <- min(plot_data$OR, na.rm = TRUE)
    
    cat("6. 绘图数据范围: OR从", min_or, "到", max_or, "\n")
    
    # 确保y轴最小值至少为0
    y_min <- 0
    y_max <- max(1, max_or * 1.4)  # 确保至少能看到柱子
    
    # 绘制柱状图 - 强制显示
    p_or <- ggplot(plot_data, aes(x = combination_label, y = OR, fill = comparison_type)) +
      geom_col(width = 0.7, color = "black", alpha = 0.8) +
      # 添加数值标签
      geom_text(aes(label = sprintf("OR=%.2f", OR)), 
                vjust = -0.3, size = 2, color = "black") +
      # 添加显著性标记
      geom_text(aes(label = significance), 
                vjust = -1.2, size = 3, color = "#F1DD10") +
      scale_fill_manual(values = comparison_colors,
                        name = "Comparison Type",
                        labels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")) +
      labs(x = "Assay Combination",
           y = "Odds Ratio (OR)",
           title = "Odds Ratios for Allosteric Site Co-occurrence\n per site (GTP pocket only)") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(),
        panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(hjust = 0.5,  size = 10),
        legend.position = "bottom"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      coord_cartesian(ylim = c(y_min, y_max))
    
    # 强制显示图形
    cat("7. 正在显示图形...\n")
    print(p_or)
    cat("8. 图形显示完成\n")
    
    # 保存图形
    cat("9. 正在保存图形...\n")
    ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf5/20251101/OR_comparison_predefined_allosteric_sites_per_site_GTP_pocket_only.pdf", 
           p_or, device = cairo_pdf, width = 14, height = 10)
    
    # 同时保存PNG格式
    ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf5/20251101/OR_comparison_predefined_allosteric_sites_per_site_GTP_pocket_only.png", 
           p_or, width = 14, height = 10, dpi = 300, bg = "white")
    
    cat("10. 图形保存完成\n")
  }
  
} else {
  cat("!!! 错误: or_results 数据框为空 !!!\n")
}

cat("=== 绘图过程结束 ===\n")

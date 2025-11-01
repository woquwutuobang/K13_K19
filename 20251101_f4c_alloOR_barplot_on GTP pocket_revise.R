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

# 定义位点信息
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)
PI3KCG_binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
GTP_Binding_pocket_site <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

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

anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"

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

calculate_weighted_mean_ddG <- function(ddG) {
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  wt_codon_info <- ddG[Pos_real > 1, .(codon = first(wt_codon)), by = "Pos_real"]
  final_output <- merge(wt_codon_info, weighted_mean_ddG, by = "Pos_real")
  
  return(final_output)
}

classify_site_mutation_types <- function(ddG, weighted_mean_ddG, anno, assay_sele, reg_threshold = NULL) {
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  if (is.null(reg_threshold)) {
    reg_threshold <- data_plot[binding_type == "binding site",
                               sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                               .SDcols = c("mean", "sigma")]
    cat("Calculated regulatory threshold for", assay_sele, ":", reg_threshold, "\n")
  }
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  return(data_plot_mutation)
}

# ===============================
# 计算单个assay对的OR值（分别处理每个assay对，只保留GTP pocket）
# ===============================
calculate_pair_or_pvalue_gtp <- function(assay1, assay2) {
  cat("处理", assay1, "vs", assay2, "在GTP pocket上的OR...\n")
  
  # 分别读取两个assay的数据
  input1 <- get(paste0("input_", assay1))
  input2 <- get(paste0("input_", assay2))
  
  # 处理assay1数据
  ddG_data1 <- read_ddG_data(input1, wt_aa)
  ddG_weighted1 <- calculate_weighted_mean_ddG(ddG_data1)
  anno <- fread(anno_file)
  result1 <- classify_site_mutation_types(ddG_data1, ddG_weighted1, anno, assay1)
  
  # 处理assay2数据
  ddG_data2 <- read_ddG_data(input2, wt_aa)
  ddG_weighted2 <- calculate_weighted_mean_ddG(ddG_data2)
  result2 <- classify_site_mutation_types(ddG_data2, ddG_weighted2, anno, assay2)
  
  # 合并这两个assay的数据
  data1 <- result1[, c("Pos_real", "wt_codon", "mt_codon", "allosteric_mutation")]
  setnames(data1, "allosteric_mutation", paste0("allosteric_mutation_", assay1))
  
  data2 <- result2[, c("Pos_real", "wt_codon", "mt_codon", "allosteric_mutation")]
  setnames(data2, "allosteric_mutation", paste0("allosteric_mutation_", assay2))
  
  pair_data <- merge(data1, data2, by = c("Pos_real", "wt_codon", "mt_codon"), all = FALSE)
  
  # 只保留GTP pocket位点
  filtered_data <- pair_data[Pos_real %in% GTP_Binding_pocket_site]
  
  # 计算OR值
  is_effect1 <- filtered_data[[paste0("allosteric_mutation_", assay1)]] == TRUE
  is_effect2 <- filtered_data[[paste0("allosteric_mutation_", assay2)]] == TRUE
  
  contingency_table <- table(is_effect1, is_effect2)
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
    n_mutations = nrow(filtered_data),
    n_effect1 = sum(is_effect1),
    n_effect2 = sum(is_effect2),
    n_both = sum(is_effect1 & is_effect2)
  ))
}

# ===============================
# 主分析流程
# ===============================
cat("=== 开始分析所有assay组合在GTP pocket上的OR值 ===\n")

# 定义所有assay
assay_names <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# 计算所有组合的OR值
or_results <- data.table()
combinations <- combn(assay_names, 2, simplify = FALSE)

for (comb in combinations) {
  assay1 <- comb[1]
  assay2 <- comb[2]
  
  result <- calculate_pair_or_pvalue_gtp(assay1, assay2)
  or_results <- rbind(or_results, result)
}

# 添加显著性标记
or_results[, significance := fcase(
  p_value < 0.001, "***",
  p_value < 0.01, "**",
  p_value < 0.05, "*",
  default = "NS"
)]

# 打印结果
cat("\n=== GTP pocket OR值结果 ===\n")
print(or_results[order(comparison_type, OR)])

# ===============================
# 绘制柱状图（带OR值标注）
# ===============================
cat("\n=== 绘制GTP pocket OR柱状图 ===\n")

# 设置颜色
comparison_colors <- c(
  "BI1_vs_BI1" = "#F4270C",
  "BI1_vs_BI2" = "#F4AD0C", 
  "BI2_vs_BI2" = "#1B38A6"
)

# 创建组合标签
or_results[, combination_label := paste0(assay1, " vs ", assay2)]

# 按照比较类型分区，然后在每个分区内按OR值排序
or_results[, comparison_type := factor(comparison_type, 
                                       levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2"))]

# 在每个比较类型内按OR值排序
or_results[, within_group_order := frank(OR, ties.method = "first"), by = comparison_type]
or_results[, final_order := as.numeric(comparison_type) * 100 + within_group_order]

# 重新排序数据
or_results <- or_results[order(final_order)]
or_results[, combination_label := factor(combination_label, levels = unique(combination_label[order(final_order)]))]

# 计算最大OR值用于调整标注位置
max_or <- max(or_results$OR, na.rm = TRUE)

# 绘制柱状图
p_or <- ggplot(or_results, aes(x = combination_label, y = OR, fill = comparison_type)) +
  geom_col(width = 0.7) +
  # 添加OR值标注
  geom_text(aes(label = sprintf("OR=%.2f", OR)), 
            vjust = -0.5, size = 2, color = "black") +
  # 添加显著性标记
  geom_text(aes(label = significance), 
            vjust = -1.5, size = 3, color = "#F1DD10") +
  scale_fill_manual(values = comparison_colors,
                    name = "Comparison Type",
                    labels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")) +
  labs(x = "Assay Combination",
       y = "Odds Ratio (OR)",
       title = "Odds Ratios for Allosteric Mutation Co-occurrence in GTP Pocket\n per mutation") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, max_or * 1.4))  # 增加y轴上限以容纳标注

# 显示图形
print(p_or)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/OR_comparison_GTP_pocket_separate_analysis_per_mutation_new.pdf", 
       p_or, device = cairo_pdf, width = 12, height = 8)





# 保存OR结果到文件
fwrite(or_results, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/OR_results_GTP_pocket_separate_analysis.csv")

cat("分析完成！\n")
cat("图形已保存\n")
cat("结果数据已保存\n")

# 输出统计摘要
cat("\n=== 统计摘要 ===\n")
cat("GTP pocket位点数:", length(GTP_Binding_pocket_site), "\n")
cat("BI1 assays: RAF1, RALGDS, PI3KCG, SOS1, K55, K27\n")
cat("BI2 assays: K13, K19\n")
cat("比较类型分布:\n")
print(table(or_results$comparison_type))
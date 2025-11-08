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

#anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"

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
# 计算单个assay对的OR值（基于预定义的变构位点，去除被比较assay的界面）
# ===============================
calculate_pair_or_pvalue_allosteric <- function(assay1, assay2) {
  cat("处理", assay1, "vs", assay2, "基于预定义变构位点...\n")
  
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
  
  # 移除被比较的两个assay的结合界面位点
  sites_to_remove <- unique(c(binding_sites[[assay1]], binding_sites[[assay2]]))
  filtered_data <- pair_data[!Pos_real %in% sites_to_remove]
  
  # 计算OR值
  is_effect1 <- filtered_data[[paste0("is_allosteric_", assay1)]] == TRUE
  is_effect2 <- filtered_data[[paste0("is_allosteric_", assay2)]] == TRUE
  
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
    n_sites = nrow(filtered_data),
    n_effect1 = sum(is_effect1),
    n_effect2 = sum(is_effect2),
    n_both = sum(is_effect1 & is_effect2),
    sites_removed = length(sites_to_remove)
  ))
}

# ===============================
# 主分析流程
# ===============================
cat("=== 开始分析所有assay组合（基于预定义变构位点，去除结合界面） ===\n")

# 定义所有assay
assay_names <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# 计算所有组合的OR值
or_results <- data.table()
combinations <- combn(assay_names, 2, simplify = FALSE)

for (comb in combinations) {
  assay1 <- comb[1]
  assay2 <- comb[2]
  
  result <- calculate_pair_or_pvalue_allosteric(assay1, assay2)
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
cat("\n=== 基于预定义变构位点的OR值结果（去除结合界面） ===\n")
print(or_results[order(comparison_type, OR)])


####===========================================================
## plot2

# ============================================================
# OR Results Visualization with Multi-layer Sorting (No OR-based sorting)
# ============================================================

# ===============================
# 分类定义
# ===============================
binders <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

non_darpin <- c("RAF1", "RALGDS", "PI3KCG", "SOS1")
darpin     <- c("K55", "K27", "K13", "K19")

active     <- c("RAF1", "RALGDS", "PI3KCG", "K55")
inactive   <- c("SOS1","K27", "K13", "K19")

# 定义BI1和BI2
BI1 <- c("RAF1", "RALGDS", "PI3KCG", "SOS1","K55", "K27")
BI2 <- c("K13", "K19")

# ===============================
# 数据准备
# ===============================
or_results$Comparison <- paste0(or_results$assay1, " vs ", or_results$assay2)

# 分类信息
or_results$Type1 <- ifelse(or_results$assay1 %in% darpin, "DARPin", "non-DARPin")
or_results$Type2 <- ifelse(or_results$assay2 %in% darpin, "DARPin", "non-DARPin")

or_results$Act1 <- ifelse(or_results$assay1 %in% active, "Active", "Inactive")
or_results$Act2 <- ifelse(or_results$assay2 %in% active, "Active", "Inactive")

# ===============================
# 层级标签生成
# ===============================
or_results$DARPinClass <- with(
  or_results,
  ifelse(Type1 == "non-DARPin" & Type2 == "non-DARPin", "nonDARPin_vs_nonDARPin",
         ifelse(Type1 == "DARPin" & Type2 == "DARPin", "DARPin_vs_DARPin", "nonDARPin_vs_DARPin"))
)

or_results$ActivityClass <- with(
  or_results,
  ifelse(Act1 == "Active" & Act2 == "Active", "Active_vs_Active",
         ifelse(Act1 == "Inactive" & Act2 == "Inactive", "Inactive_vs_Inactive", "Active_vs_Inactive"))
)

# 生成comparison_type分类
or_results$comparison_type <- apply(or_results, 1, function(x) {
  assay1 <- x["assay1"]
  assay2 <- x["assay2"]
  
  if (assay1 %in% BI1 && assay2 %in% BI1) {
    return("BI1_vs_BI1")
  } else if (assay1 %in% BI2 && assay2 %in% BI2) {
    return("BI2_vs_BI2")
  } else {
    return("BI1_vs_BI2")
  }
})

# ===============================
# 三层层级排序
# ===============================
or_results <- or_results[order(
  # 第一层：comparison_type (BI1_vs_BI1, BI1_vs_BI2, BI2_vs_BI2)
  factor(or_results$comparison_type, levels = c(
    "BI1_vs_BI1",
    "BI1_vs_BI2", 
    "BI2_vs_BI2"
  )),
  # 第二层：DARPin分类
  factor(or_results$DARPinClass, levels = c(
    "nonDARPin_vs_nonDARPin",
    "nonDARPin_vs_DARPin",
    "DARPin_vs_DARPin"
  )),
  # 第三层：活性分类
  factor(or_results$ActivityClass, levels = c(
    "Active_vs_Active",
    "Active_vs_Inactive", 
    "Inactive_vs_Inactive"
  )),
  # 最后按比较名称排序以确保一致性
  factor(or_results$Comparison, levels = unique(or_results$Comparison))
), ]

# 设置因子顺序
or_results$Comparison <- factor(or_results$Comparison, levels = or_results$Comparison)

# 设置comparison_type因子顺序
or_results$comparison_type <- factor(
  or_results$comparison_type, 
  levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2")
)

# ===============================
# 绘图
# ===============================
p_or <- ggplot(or_results, aes(x = Comparison, y = OR, fill = comparison_type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("OR=%.2f", OR)), vjust = -0.5, size = 2, color = "black") +
  geom_text(aes(label = significance), vjust = -1.5, size = 2.5, color = "#F1DD10") +
  scale_fill_manual(values = c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C",
    "BI2_vs_BI2" = "#1B38A6"
  )) +
  labs(
    title = "Odds Ratios for Allosteric site Co-occurrence (delete 2BIs)\nper site",
    x = "Binder Pair Comparison",
    y = "Odds Ratio (OR)",
    fill = "Comparison Type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  ylim(0, max(or_results$OR, na.rm = TRUE) * 1.2)

print(p_or)

# ===============================
# 保存图表
# ===============================
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251106/Binder_OR_AllostericSites_multilayer_sorted per site logic sort all.pdf",
  plot = p_or,
  device = cairo_pdf,
  width = 12,
  height = 8
)


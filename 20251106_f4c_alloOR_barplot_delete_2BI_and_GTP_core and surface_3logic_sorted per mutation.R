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
core_surface_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv"

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
# 计算单个assay对的OR值（分别处理每个assay对）
# ===============================
calculate_pair_or_pvalue <- function(assay1, assay2, region_type) {
  cat("处理", assay1, "vs", assay2, "-", region_type, "...\n")
  
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
  
  # 读取core/surface信息并合并
  core_surface <- fread(core_surface_file)
  core_surface <- core_surface[, Pos_real := Pos]
  core_surface_threshold <- 0.25
  core_surface <- core_surface %>%
    mutate(type = case_when(
      RASA <= core_surface_threshold ~ "core",
      RASA > core_surface_threshold ~ "surface"
    )) %>%
    filter(!is.na(type))
  
  pair_data <- merge(pair_data, core_surface[, .(Pos_real, type)], by = "Pos_real", all.x = TRUE)
  
  # 移除两个assay各自的结合界面位点和GTP pocket
  sites_to_remove <- unique(c(binding_sites[[assay1]], binding_sites[[assay2]], GTP_Binding_pocket_site))
  
  # 根据区域类型筛选数据
  if (region_type == "core") {
    filtered_data <- pair_data[!Pos_real %in% sites_to_remove & type == "core"]
  } else if (region_type == "surface") {
    filtered_data <- pair_data[!Pos_real %in% sites_to_remove & type == "surface"]
  }
  
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
    region_type = region_type,
    OR = or_value,
    p_value = p_value,
    n_mutations = nrow(filtered_data),
    n_effect1 = sum(is_effect1),
    n_effect2 = sum(is_effect2),
    n_both = sum(is_effect1 & is_effect2),
    sites_removed = length(sites_to_remove)
  ))
}

# ===============================
# 主分析流程
# ===============================
cat("=== 开始分析所有assay组合 ===\n")

# 定义所有assay
assay_names <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# 计算所有组合的OR值（分区域）
or_results_core <- data.table()
or_results_surface <- data.table()

combinations <- combn(assay_names, 2, simplify = FALSE)

for (comb in combinations) {
  assay1 <- comb[1]
  assay2 <- comb[2]
  
  # 计算core区域
  result_core <- calculate_pair_or_pvalue(assay1, assay2, "core")
  or_results_core <- rbind(or_results_core, result_core)
  
  # 计算surface区域
  result_surface <- calculate_pair_or_pvalue(assay1, assay2, "surface")
  or_results_surface <- rbind(or_results_surface, result_surface)
}

# 为每个结果集添加显著性标记
add_significance <- function(df) {
  df[, significance := fcase(
    p_value < 0.001, "***",
    p_value < 0.01, "**",
    p_value < 0.05, "*",
    default = "NS"
  )]
  return(df)
}

or_results_core <- add_significance(or_results_core)
or_results_surface <- add_significance(or_results_surface)

####===========================================================
## plot2_core

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
or_results_core$Comparison <- paste0(or_results_core$assay1, " vs ", or_results_core$assay2)

# 分类信息
or_results_core$Type1 <- ifelse(or_results_core$assay1 %in% darpin, "DARPin", "non-DARPin")
or_results_core$Type2 <- ifelse(or_results_core$assay2 %in% darpin, "DARPin", "non-DARPin")

or_results_core$Act1 <- ifelse(or_results_core$assay1 %in% active, "Active", "Inactive")
or_results_core$Act2 <- ifelse(or_results_core$assay2 %in% active, "Active", "Inactive")

# ===============================
# 层级标签生成
# ===============================
or_results_core$DARPinClass <- with(
  or_results_core,
  ifelse(Type1 == "non-DARPin" & Type2 == "non-DARPin", "nonDARPin_vs_nonDARPin",
         ifelse(Type1 == "DARPin" & Type2 == "DARPin", "DARPin_vs_DARPin", "nonDARPin_vs_DARPin"))
)

or_results_core$ActivityClass <- with(
  or_results_core,
  ifelse(Act1 == "Active" & Act2 == "Active", "Active_vs_Active",
         ifelse(Act1 == "Inactive" & Act2 == "Inactive", "Inactive_vs_Inactive", "Active_vs_Inactive"))
)

# 生成comparison_type分类
or_results_core$comparison_type <- apply(or_results_core, 1, function(x) {
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
or_results_core <- or_results_core[order(
  # 第一层：comparison_type (BI1_vs_BI1, BI1_vs_BI2, BI2_vs_BI2)
  factor(or_results_core$comparison_type, levels = c(
    "BI1_vs_BI1",
    "BI1_vs_BI2", 
    "BI2_vs_BI2"
  )),
  # 第二层：DARPin分类
  factor(or_results_core$DARPinClass, levels = c(
    "nonDARPin_vs_nonDARPin",
    "nonDARPin_vs_DARPin",
    "DARPin_vs_DARPin"
  )),
  # 第三层：活性分类
  factor(or_results_core$ActivityClass, levels = c(
    "Active_vs_Active",
    "Active_vs_Inactive", 
    "Inactive_vs_Inactive"
  )),
  # 最后按比较名称排序以确保一致性
  factor(or_results_core$Comparison, levels = unique(or_results_core$Comparison))
), ]

# 设置因子顺序
or_results_core$Comparison <- factor(or_results_core$Comparison, levels = or_results_core$Comparison)

# 设置comparison_type因子顺序
or_results_core$comparison_type <- factor(
  or_results_core$comparison_type, 
  levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2")
)

# ===============================
# 绘图
# ===============================
p_or <- ggplot(or_results_core, aes(x = Comparison, y = OR, fill = comparison_type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("OR=%.2f", OR)), vjust = -0.5, size = 2, color = "black") +
  geom_text(aes(label = significance), vjust = -1.5, size = 2.5, color = "#F1DD10") +
  scale_fill_manual(values = c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C",
    "BI2_vs_BI2" = "#1B38A6"
  )) +
  labs(
    title = "Odds Ratios for Allosteric Mutation Co-occurrence (delete 2BIs and GTP pocket)\nper mutation-core region",
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
  ylim(0, max(or_results_core$OR, na.rm = TRUE) * 1.2)

print(p_or)

# ===============================
# 保存图表
# ===============================
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251106/Binder_OR_AllostericSites_multilayer_sorted per mutation logic sort_delete 2BIs and GTP pocket core.pdf",
  plot = p_or,
  device = cairo_pdf,
  width = 12,
  height = 8
)

####===========================================================
## plot2_surface

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
or_results_surface$Comparison <- paste0(or_results_surface$assay1, " vs ", or_results_surface$assay2)

# 分类信息
or_results_surface$Type1 <- ifelse(or_results_surface$assay1 %in% darpin, "DARPin", "non-DARPin")
or_results_surface$Type2 <- ifelse(or_results_surface$assay2 %in% darpin, "DARPin", "non-DARPin")

or_results_surface$Act1 <- ifelse(or_results_surface$assay1 %in% active, "Active", "Inactive")
or_results_surface$Act2 <- ifelse(or_results_surface$assay2 %in% active, "Active", "Inactive")

# ===============================
# 层级标签生成
# ===============================
or_results_surface$DARPinClass <- with(
  or_results_surface,
  ifelse(Type1 == "non-DARPin" & Type2 == "non-DARPin", "nonDARPin_vs_nonDARPin",
         ifelse(Type1 == "DARPin" & Type2 == "DARPin", "DARPin_vs_DARPin", "nonDARPin_vs_DARPin"))
)

or_results_surface$ActivityClass <- with(
  or_results_surface,
  ifelse(Act1 == "Active" & Act2 == "Active", "Active_vs_Active",
         ifelse(Act1 == "Inactive" & Act2 == "Inactive", "Inactive_vs_Inactive", "Active_vs_Inactive"))
)

# 生成comparison_type分类
or_results_surface$comparison_type <- apply(or_results_surface, 1, function(x) {
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
or_results_surface <- or_results_surface[order(
  # 第一层：comparison_type (BI1_vs_BI1, BI1_vs_BI2, BI2_vs_BI2)
  factor(or_results_surface$comparison_type, levels = c(
    "BI1_vs_BI1",
    "BI1_vs_BI2", 
    "BI2_vs_BI2"
  )),
  # 第二层：DARPin分类
  factor(or_results_surface$DARPinClass, levels = c(
    "nonDARPin_vs_nonDARPin",
    "nonDARPin_vs_DARPin",
    "DARPin_vs_DARPin"
  )),
  # 第三层：活性分类
  factor(or_results_surface$ActivityClass, levels = c(
    "Active_vs_Active",
    "Active_vs_Inactive", 
    "Inactive_vs_Inactive"
  )),
  # 最后按比较名称排序以确保一致性
  factor(or_results_surface$Comparison, levels = unique(or_results_surface$Comparison))
), ]

# 设置因子顺序
or_results_surface$Comparison <- factor(or_results_surface$Comparison, levels = or_results_surface$Comparison)

# 设置comparison_type因子顺序
or_results_surface$comparison_type <- factor(
  or_results_surface$comparison_type, 
  levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2")
)

# ===============================
# 绘图
# ===============================
p_or <- ggplot(or_results_surface, aes(x = Comparison, y = OR, fill = comparison_type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("OR=%.2f", OR)), vjust = -0.5, size = 2, color = "black") +
  geom_text(aes(label = significance), vjust = -1.5, size = 2.5, color = "#F1DD10") +
  scale_fill_manual(values = c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C",
    "BI2_vs_BI2" = "#1B38A6"
  )) +
  labs(
    title = "Odds Ratios for Allosteric Mutation Co-occurrence (delete 2BIs and GTP pocket)\nper mutation-surface region",
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
  ylim(0, max(or_results_surface$OR, na.rm = TRUE) * 1.2)

print(p_or)

# ===============================
# 保存图表
# ===============================
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251106/Binder_OR_AllostericSites_multilayer_sorted per mutation logic sort_delete 2BIs and GTP pocket surface.pdf",
  plot = p_or,
  device = cairo_pdf,
  width = 12,
  height = 8
)

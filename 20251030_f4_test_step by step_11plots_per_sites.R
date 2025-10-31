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
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
GTP_Binding_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# 新增allosteric位点
K13_allosteric_site <- c(15, 16, 17, 35, 145, 10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)

# 文件路径
input1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt"
input2 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
core_surface_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv"
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251030"

# 坐标轴范围
xlim_all = c(-0.7, 2.3)
ylim_all = c(-0.7, 1.8)
xlim_effect = c(-0.7, 2.3)
ylim_effect = c(-0.7, 1.8)

# ===============================
# 数据预处理函数 - 修改为使用median
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

classify_site_mutation_types <- function(ddG, median_ddG, anno, assay_sele, reg_threshold = NULL) {
  data_plot <- merge(median_ddG, anno, by = "Pos", all = TRUE)
  
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  if (is.null(reg_threshold)) {
    reg_threshold <- data_plot[binding_type == "binding site",
                               median(abs(median_ddG), na.rm = TRUE)]
    cat("Calculated regulatory threshold for", assay_sele, ":", reg_threshold, "\n")
  }
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # 添加allosteric site分类
  if (assay_sele == "K13") {
    data_plot[Pos %in% K13_allosteric_site, site_type := "Allosteric site"]
  } else if (assay_sele == "RAF1") {
    data_plot[Pos %in% RAF1_allosteric_site, site_type := "Allosteric site"]
  }
  
  # 合并突变数据与位点类型信息
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # 识别变构突变 - 使用median-based方法
  data_plot_mutation[, allosteric_mutation := FALSE]
  if (assay_sele == "K13") {
    data_plot_mutation[Pos %in% K13_allosteric_site, allosteric_mutation := TRUE]
  } else if (assay_sele == "RAF1") {
    data_plot_mutation[Pos %in% RAF1_allosteric_site, allosteric_mutation := TRUE]
  }
  
  # 分类突变类型
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Orthosteric site huge differences"]
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "Orthosteric site small differences"]
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "GTP binding allosteric mutation"]
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "GTP binding other mutation"]
  data_plot_mutation[Pos %in% data_plot[site_type == "Allosteric site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Allosteric mutation"]
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site", "Allosteric site") &
                       allosteric_mutation == FALSE, mutation_type := "Other mutation"]
  
  data_plot_mutation <- within(data_plot_mutation,
                               mutation_type <- factor(mutation_type,
                                                       levels = c("Orthosteric site huge differences",
                                                                  "Orthosteric site small differences",
                                                                  "GTP binding allosteric mutation",
                                                                  "GTP binding other mutation",
                                                                  "Allosteric mutation",
                                                                  "Other mutation")))
  
  return(data_plot_mutation)
}

# ===============================
# 数据预处理
# ===============================
cat("=== 数据预处理 ===\n")

# 读取RAF1数据
cat("处理 RAF1 数据...\n")
ddG_data1 <- read_ddG_data(input1, wt_aa)
ddG_median1 <- calculate_median_ddG(ddG_data1)

# 读取K13数据
cat("处理 K13 数据...\n")
ddG_data2 <- read_ddG_data(input2, wt_aa)
ddG_median2 <- calculate_median_ddG(ddG_data2)

# 读取注释文件
anno <- fread(anno_file)

# 分类突变类型
result1 <- classify_site_mutation_types(ddG_data1, ddG_median1, anno, "RAF1")
result2 <- classify_site_mutation_types(ddG_data2, ddG_median2, anno, "K13")

# 提取需要的列并合并
result1 <- result1[, c("Pos_real", "wt_codon", "mt_codon", "mean_kcal/mol", "std_kcal/mol", 
                       "allosteric_mutation", "mutation_type", "site_type")]
result2 <- result2[, c("Pos_real", "wt_codon", "mt_codon", "mean_kcal/mol", "std_kcal/mol", 
                       "allosteric_mutation", "mutation_type", "site_type")]

data_plot_mutation1 <- merge(result1, result2, by = c("Pos_real", "wt_codon", "mt_codon"), 
                             suffixes = c("_RAF1", "_K13"), all = FALSE)

# 读取core/surface信息
core_surface <- fread(core_surface_file)
core_surface <- core_surface[, Pos_real := Pos]
core_surface_threshold <- 0.25
core_surface <- core_surface %>%
  mutate(type = case_when(
    RASA <= core_surface_threshold ~ "core",
    RASA > core_surface_threshold ~ "surface"
  )) %>%
  filter(!is.na(type))

data_plot_mutation <- merge(data_plot_mutation1, core_surface, by = "Pos_real", all = TRUE)

# ===============================
# 创建位点级别的数据（使用median）
# ===============================
cat("\n=== 创建位点级别数据（基于median） ===\n")

# 计算每个位点的median值
site_level_data_RAF1 <- ddG_data1[Pos_real > 1, .(
  median_ddG_RAF1 = median(`mean_kcal/mol`, na.rm = TRUE),
  n_mutations_RAF1 = sum(!is.na(`mean_kcal/mol`))
), by = "Pos_real"]

site_level_data_K13 <- ddG_data2[Pos_real > 1, .(
  median_ddG_K13 = median(`mean_kcal/mol`, na.rm = TRUE),
  n_mutations_K13 = sum(!is.na(`mean_kcal/mol`))
), by = "Pos_real"]

# 合并位点级别数据
site_level_data <- merge(site_level_data_RAF1, site_level_data_K13, by = "Pos_real", all = TRUE)

# 添加位点类型信息
site_level_data[, is_RAF1_interface := Pos_real %in% RAF1_Binding_interface_site]
site_level_data[, is_K13_interface := Pos_real %in% K13_Binding_interface_site]
site_level_data[, is_gtp_pocket := Pos_real %in% GTP_Binding_pocket]
site_level_data[, is_RAF1_allosteric := Pos_real %in% RAF1_allosteric_site]
site_level_data[, is_K13_allosteric := Pos_real %in% K13_allosteric_site]

# 添加core/surface信息
site_level_data <- merge(site_level_data, core_surface[, .(Pos_real, type)], by = "Pos_real", all.x = TRUE)

# ===============================
# 图1: 所有位点
# ===============================
cat("\n=== 绘制图1: 所有位点 ===\n")

p1 <- ggplot(site_level_data, aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(color = "grey70", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "All sites") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p1
ggsave(file.path(output_dir, "p1_all_sites_median.pdf"), p1, 
       device = cairo_pdf, width = 4, height = 4)
cat("图1保存完成\n")

# ===============================
# 图2: RAF1结合界面位点
# ===============================
cat("\n=== 绘制图2: RAF1结合界面位点 ===\n")

p2 <- ggplot(site_level_data[is_RAF1_interface == TRUE], 
             aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(color = "#F4270C", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "RAF1 Binding Interface") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
p2
ggsave(file.path(output_dir, "p2_RAF1_binding_interface_median.pdf"), p2, 
       device = cairo_pdf, width = 4, height = 4)
cat("图2保存完成\n")

# ===============================
# 图3: K13结合界面位点
# ===============================
cat("\n=== 绘制图3: K13结合界面位点 ===\n")

p3 <- ggplot(site_level_data[is_K13_interface == TRUE], 
             aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(color = "#1B38A6", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "K13 Binding Interface") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
p3
ggsave(file.path(output_dir, "p3_K13_binding_interface_median.pdf"), p3, 
       device = cairo_pdf, width = 4, height = 4)
cat("图3保存完成\n")

# ===============================
# 准备变构效应分析数据 - 位点级别
# ===============================
cat("\n=== 准备变构效应分析数据 - 位点级别 ===\n")

# 移除结合界面位点的位点
remaining_sites <- site_level_data[!Pos_real %in% unique(c(RAF1_Binding_interface_site, K13_Binding_interface_site))]

# 基于位点的功能效应判断
remaining_sites[, allosteric_type_by_effect := "Other"]
remaining_sites[is_RAF1_allosteric == TRUE & is_K13_allosteric == FALSE, 
                allosteric_type_by_effect := "RAF1_allosteric_only"]
remaining_sites[is_RAF1_allosteric == FALSE & is_K13_allosteric == TRUE, 
                allosteric_type_by_effect := "K13_allosteric_only"]
remaining_sites[is_RAF1_allosteric == TRUE & is_K13_allosteric == TRUE, 
                allosteric_type_by_effect := "Both_allosteric"]

# 创建颜色映射
allosteric_colors <- c("#F4AD0C", "#FFB0A5", "#C68EFD", "grey70")
names(allosteric_colors) <- c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other")

# ===============================
# 图4: 所有位点（排除界面后）
# ===============================
cat("\n=== 绘制图4: 所有位点（排除界面后） ===\n")

calculate_or_pvalue_site <- function(data) {
  is_effect1 <- data$is_RAF1_allosteric == TRUE
  is_effect2 <- data$is_K13_allosteric == TRUE
  
  contingency_table <- table(is_effect1, is_effect2)
  fisher_test <- fisher.test(contingency_table)
  
  or_value <- round(fisher_test$estimate, 2)
  p_value <- fisher_test$p.value
  
  if (p_value < 0.001) {
    p_label <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_label <- "p < 0.01" 
  } else if (p_value < 0.05) {
    p_label <- "p < 0.05"
  } else {
    p_label <- paste0("p = ", round(p_value, 3))
  }
  
  cat("Contingency Table:\n")
  print(contingency_table)
  cat("Total sites:", nrow(data), "\n")
  cat("RAF1 allosteric sites:", sum(is_effect1), "\n")
  cat("K13 allosteric sites:", sum(is_effect2), "\n")
  cat("Both allosteric sites:", sum(is_effect1 & is_effect2), "\n")
  cat("OR =", or_value, ",", p_label, "\n\n")
  
  return(list(or = or_value, p_label = p_label, contingency = contingency_table))
}

or_result_all <- calculate_or_pvalue_site(remaining_sites)
or_label_all <- paste0("OR = ", or_result_all$or, "\n", or_result_all$p_label)

p4 <- ggplot(remaining_sites, aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type_by_effect), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "All sites outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_all, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p4
ggsave(file.path(output_dir, "p4_all_sites_outside_interfaces_median.pdf"), p4, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图4保存完成\n")

# ===============================
# 图5: Core位点（排除界面后）
# ===============================
cat("\n=== 绘制图5: Core位点（排除界面后） ===\n")

core_sites <- remaining_sites[type == "core"]
or_result_core <- calculate_or_pvalue_site(core_sites)
or_label_core <- paste0("OR = ", or_result_core$or, "\n", or_result_core$p_label)

p5 <- ggplot(core_sites, aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type_by_effect), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Core sites outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_core, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p5
ggsave(file.path(output_dir, "p5_core_sites_outside_interfaces_median.pdf"), p5, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图5保存完成\n")

# ===============================
# 图6: Surface位点（排除界面后）
# ===============================
cat("\n=== 绘制图6: Surface位点（排除界面后） ===\n")

surface_sites <- remaining_sites[type == "surface"]
or_result_surface <- calculate_or_pvalue_site(surface_sites)
or_label_surface <- paste0("OR = ", or_result_surface$or, "\n", or_result_surface$p_label)

p6 <- ggplot(surface_sites, aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type_by_effect), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Surface sites outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_surface, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p6
ggsave(file.path(output_dir, "p6_surface_sites_outside_interfaces_median.pdf"), p6, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图6保存完成\n")

# ===============================
# 图7: GTP位点着色
# ===============================
cat("\n=== 绘制图7: GTP位点着色 ===\n")

p7 <- ggplot(site_level_data, aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(data = site_level_data[is_gtp_pocket == TRUE], 
             color = "#09B636", size = 2.5, alpha = 0.8) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "GTP pocket sites") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
p7
ggsave(file.path(output_dir, "p7_gtp_pocket_sites_median.pdf"), p7, 
       device = cairo_pdf, width = 4, height = 4)
cat("图7保存完成\n")

# ===============================
# 图8: 排除所有界面后的变构位点OR值计算
# ===============================
cat("\n=== 绘制图8: 排除所有界面后的变构位点OR值计算 ===\n")

excluded_sites <- unique(c(RAF1_Binding_interface_site, K13_Binding_interface_site, GTP_Binding_pocket))
remaining_sites_clean <- site_level_data[!Pos_real %in% excluded_sites]

# 分类变构位点类型
remaining_sites_clean[, allosteric_type := "Other"]
remaining_sites_clean[is_RAF1_allosteric == TRUE & is_K13_allosteric == FALSE, 
                      allosteric_type := "RAF1_allosteric_only"]
remaining_sites_clean[is_RAF1_allosteric == FALSE & is_K13_allosteric == TRUE, 
                      allosteric_type := "K13_allosteric_only"]
remaining_sites_clean[is_RAF1_allosteric == TRUE & is_K13_allosteric == TRUE, 
                      allosteric_type := "Both_allosteric"]

# 计算OR值
or_result_clean <- calculate_or_pvalue_site(remaining_sites_clean)
or_label_clean <- paste0("OR = ", or_result_clean$or, "\n", or_result_clean$p_label)

p8 <- ggplot(remaining_sites_clean, 
             aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Allosteric sites (exclude all interfaces)") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_clean, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p8
ggsave(file.path(output_dir, "p8_allosteric_exclude_all_interfaces_median.pdf"), p8, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图8保存完成\n")

# ===============================
# 图9: Core位点变构（排除所有界面）
# ===============================
cat("\n=== 绘制图9: Core位点变构（排除所有界面） ===\n")

core_sites_clean <- remaining_sites_clean[type == "core"]
or_result_core_clean <- calculate_or_pvalue_site(core_sites_clean)
or_label_core_clean <- paste0("OR = ", or_result_core_clean$or, "\n", or_result_core_clean$p_label)

p9 <- ggplot(core_sites_clean, 
             aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Core allosteric sites (exclude all interfaces)") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_core_clean, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p9
ggsave(file.path(output_dir, "p9_core_allosteric_exclude_all_interfaces_median.pdf"), p9, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图9保存完成\n")

# ===============================
# 图10: Surface位点变构（排除所有界面）
# ===============================
cat("\n=== 绘制图10: Surface位点变构（排除所有界面） ===\n")

surface_sites_clean <- remaining_sites_clean[type == "surface"]

# 安全地计算OR值，如果无法计算则只画图
safe_calculate_or_pvalue <- function(data) {
  is_effect1 <- data$is_RAF1_allosteric == TRUE
  is_effect2 <- data$is_K13_allosteric == TRUE
  
  contingency_table <- table(is_effect1, is_effect2)
  
  # 检查contingency table是否有足够的维度
  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
    cat("Contingency table does not have enough dimensions for Fisher test\n")
    cat("Contingency Table:\n")
    print(contingency_table)
    cat("Total sites:", nrow(data), "\n")
    cat("RAF1 allosteric sites:", sum(is_effect1), "\n")
    cat("K13 allosteric sites:", sum(is_effect2), "\n")
    cat("Both allosteric sites:", sum(is_effect1 & is_effect2), "\n")
    cat("Skipping OR calculation, only plotting\n\n")
    return(NULL)
  }
  
  fisher_test <- fisher.test(contingency_table)
  
  or_value <- round(fisher_test$estimate, 2)
  p_value <- fisher_test$p.value
  
  if (p_value < 0.001) {
    p_label <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_label <- "p < 0.01" 
  } else if (p_value < 0.05) {
    p_label <- "p < 0.05"
  } else {
    p_label <- paste0("p = ", round(p_value, 3))
  }
  
  cat("Contingency Table:\n")
  print(contingency_table)
  cat("Total sites:", nrow(data), "\n")
  cat("RAF1 allosteric sites:", sum(is_effect1), "\n")
  cat("K13 allosteric sites:", sum(is_effect2), "\n")
  cat("Both allosteric sites:", sum(is_effect1 & is_effect2), "\n")
  cat("OR =", or_value, ",", p_label, "\n\n")
  
  return(list(or = or_value, p_label = p_label, contingency = contingency_table))
}

or_result_surface_clean <- safe_calculate_or_pvalue(surface_sites_clean)

# 如果有OR结果则添加标签，否则只画图
if (!is.null(or_result_surface_clean)) {
  or_label_surface_clean <- paste0("OR = ", or_result_surface_clean$or, "\n", or_result_surface_clean$p_label)
  plot_with_label <- geom_annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, 
                                   label = or_label_surface_clean, 
                                   hjust = 0, vjust = 1, size = 3.3, color = "black")
} else {
  plot_with_label <- NULL
}

p10 <- ggplot(surface_sites_clean, 
              aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Surface allosteric sites (exclude all interfaces)") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

# 如果有OR标签则添加
if (!is.null(plot_with_label)) {
  p10 <- p10 + plot_with_label
}

p10
ggsave(file.path(output_dir, "p10_surface_allosteric_exclude_all_interfaces_median.pdf"), p10, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图10保存完成\n")
# ===============================
# 图11: GTP口袋位点变构
# ===============================
cat("\n=== 绘制图11: GTP口袋位点变构 ===\n")

gtp_sites <- site_level_data[is_gtp_pocket == TRUE]

# 分类GTP口袋的变构位点类型
gtp_sites[, allosteric_type := "Other"]
gtp_sites[is_RAF1_allosteric == TRUE & is_K13_allosteric == FALSE, 
          allosteric_type := "RAF1_allosteric_only"]
gtp_sites[is_RAF1_allosteric == FALSE & is_K13_allosteric == TRUE, 
          allosteric_type := "K13_allosteric_only"]
gtp_sites[is_RAF1_allosteric == TRUE & is_K13_allosteric == TRUE, 
          allosteric_type := "Both_allosteric"]

or_result_gtp <- calculate_or_pvalue_site(gtp_sites)
or_label_gtp <- paste0("OR = ", or_result_gtp$or, "\n", or_result_gtp$p_label)

p11 <- ggplot(gtp_sites, 
              aes(x = median_ddG_RAF1, y = median_ddG_K13)) +
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Site Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 allosteric only", "K13 allosteric only", "Both allosteric", "Other")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "GTP pocket allosteric sites") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_gtp, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
p11
ggsave(file.path(output_dir, "p11_gtp_pocket_allosteric_median.pdf"), p11, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图11保存完成\n")

# ===============================
# 统计摘要
# ===============================
cat("\n=== 统计摘要 ===\n")
cat("总位点数:", nrow(site_level_data), "\n")
cat("RAF1结合界面位点数:", sum(site_level_data$is_RAF1_interface), "\n")
cat("K13结合界面位点数:", sum(site_level_data$is_K13_interface), "\n")
cat("GTP口袋位点数:", sum(site_level_data$is_gtp_pocket), "\n")
cat("RAF1变构位点数:", sum(site_level_data$is_RAF1_allosteric), "\n")
cat("K13变构位点数:", sum(site_level_data$is_K13_allosteric), "\n")
cat("共享变构位点数:", sum(site_level_data$is_RAF1_allosteric & site_level_data$is_K13_allosteric), "\n")
cat("RAF1特异性变构位点数:", sum(site_level_data$is_RAF1_allosteric & !site_level_data$is_K13_allosteric), "\n")
cat("K13特异性变构位点数:", sum(site_level_data$is_K13_allosteric & !site_level_data$is_RAF1_allosteric), "\n")
cat("排除所有界面后的位点数:", nrow(remaining_sites_clean), "\n")
cat("  - Core位点数:", nrow(core_sites_clean), "\n")
cat("  - Surface位点数:", nrow(surface_sites_clean), "\n")

cat("\n=== 所有11张位点级别图绘制完成！ ===\n")
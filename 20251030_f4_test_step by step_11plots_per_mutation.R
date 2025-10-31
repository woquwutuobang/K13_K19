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

# 文件路径
input1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt"
input2 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
core_surface_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv"
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251030"

# 坐标轴范围
xlim_all = c(-1.3, 3)
ylim_all = c(-1.7, 2.6)
xlim_effect = c(-1.5, 3)
ylim_effect = c(-1.7, 2.6)

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
  
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Orthosteric site huge differences"]
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "Orthosteric site small differences"]
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "GTP binding allosteric mutation"]
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "GTP binding other mutation"]
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == TRUE, mutation_type := "Allosteric mutation"]
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
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
ddG_weighted1 <- calculate_weighted_mean_ddG(ddG_data1)

# 读取K13数据
cat("处理 K13 数据...\n")
ddG_data2 <- read_ddG_data(input2, wt_aa)
ddG_weighted2 <- calculate_weighted_mean_ddG(ddG_data2)

# 读取注释文件
anno <- fread(anno_file)

# 分类突变类型
result1 <- classify_site_mutation_types(ddG_data1, ddG_weighted1, anno, "RAF1")
result2 <- classify_site_mutation_types(ddG_data2, ddG_weighted2, anno, "K13")

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
# 图1: 所有突变
# ===============================
cat("\n=== 绘制图1: 所有突变 ===\n")

p1 <- ggplot(data_plot_mutation, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(color = "grey70", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
       title = "All mutations") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p1
ggsave(file.path(output_dir, "p1_all_mutations.pdf"), p1, 
       device = cairo_pdf, width = 4, height = 4)
cat("图1保存完成\n")

# ===============================
# 图2: RAF1结合界面位点的突变
# ===============================
cat("\n=== 绘制图2: RAF1结合界面位点的突变 ===\n")

p2 <- ggplot(data_plot_mutation[Pos_real %in% RAF1_Binding_interface_site], 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(color = "#F4270C", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
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
ggsave(file.path(output_dir, "p2_RAF1_binding_interface.pdf"), p2, 
       device = cairo_pdf, width = 4, height = 4)
cat("图2保存完成\n")

# ===============================
# 图3: K13结合界面位点的突变
# ===============================
cat("\n=== 绘制图3: K13结合界面位点的突变 ===\n")

p3 <- ggplot(data_plot_mutation[Pos_real %in% K13_Binding_interface_site], 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(color = "#1B38A6", size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
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
ggsave(file.path(output_dir, "p3_K13_binding_interface.pdf"), p3, 
       device = cairo_pdf, width = 4, height = 4)
cat("图3保存完成\n")

# ===============================
# 准备变构效应分析数据
# ===============================
cat("\n=== 准备变构效应分析数据 ===\n")

# 移除结合界面位点的突变
remaining_mutations <- data_plot_mutation[!Pos_real %in% unique(c(RAF1_Binding_interface_site, K13_Binding_interface_site))]

# 添加GTP结合口袋标记
remaining_mutations[, is_gtp := Pos_real %in% GTP_Binding_pocket]

# 基于突变的功能效应判断
remaining_mutations[, allosteric_type_by_effect := "Other"]
remaining_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == FALSE, 
                    allosteric_type_by_effect := "RAF1_allosteric_only"]
remaining_mutations[allosteric_mutation_RAF1 == FALSE & allosteric_mutation_K13 == TRUE, 
                    allosteric_type_by_effect := "K13_allosteric_only"]
remaining_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == TRUE, 
                    allosteric_type_by_effect := "Both_allosteric"]

# 创建颜色映射
allosteric_colors <- c("#F4AD0C", "#FFB0A5", "#C68EFD", "grey70")
names(allosteric_colors) <- c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other")

# ===============================
# 图4: 所有突变（排除界面后）
# ===============================
cat("\n=== 绘制图4: 所有突变（排除界面后） ===\n")

calculate_or_pvalue <- function(data) {
  is_effect1 <- data$allosteric_mutation_RAF1 == TRUE
  is_effect2 <- data$allosteric_mutation_K13 == TRUE
  
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
  cat("Total mutations:", nrow(data), "\n")
  cat("RAF1 allosteric effect:", sum(is_effect1), "\n")
  cat("K13 allosteric effect:", sum(is_effect2), "\n")
  cat("Both effects:", sum(is_effect1 & is_effect2), "\n")
  cat("OR =", or_value, ",", p_label, "\n\n")
  
  return(paste0("OR = ", or_value, "\n", p_label))
}

or_label_all <- calculate_or_pvalue(remaining_mutations)

# 将数据分为灰色点和其他颜色点
grey_points <- remaining_mutations[allosteric_type_by_effect == "Other"]
colored_points <- remaining_mutations[allosteric_type_by_effect != "Other"]

# 进一步将彩色点分为GTP和非GTP
colored_non_gtp <- colored_points[is_gtp == FALSE]
colored_gtp <- colored_points[is_gtp == TRUE]

p4 <- ggplot(remaining_mutations, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  # 第一层：灰色非GTP点（圆形）
  geom_point(data = grey_points[is_gtp == FALSE], 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  # 第二层：灰色GTP点（三角形）
  geom_point(data = grey_points[is_gtp == TRUE], 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  # 第三层：彩色非GTP点（圆形）
  geom_point(data = colored_non_gtp, 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  # 第四层：彩色GTP点（三角形）
  geom_point(data = colored_gtp, 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),
                     labels = c("RAF1 effect only", "K13 effect only", "Both effects", "Other")) +
  scale_shape_manual(name = "GTP Binding Pocket",
                     values = c("Non-GTP pocket" = 16, "GTP pocket" = 16),
                     labels = c("GTP pocket" = "Yes", "Non-GTP pocket" = "No")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
       title = "All mutations outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_all, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )

p4
ggsave(file.path(output_dir, "p4_all_mutations_outside_interfaces 2.pdf"), p4, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图4保存完成\n")

# ===============================
# 图5: Core位点的突变（排除界面后）
# ===============================
cat("\n=== 绘制图5: Core位点的突变（排除界面后） ===\n")

core_mutations <- remaining_mutations[type == "core"]
or_label_core <- calculate_or_pvalue(core_mutations)

grey_points_core <- core_mutations[allosteric_type_by_effect == "Other"]
colored_points_core <- core_mutations[allosteric_type_by_effect != "Other"]
colored_non_gtp_core <- colored_points_core[is_gtp == FALSE]
colored_gtp_core <- colored_points_core[is_gtp == TRUE]

p5 <- ggplot(core_mutations, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(data = grey_points_core[is_gtp == FALSE], 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = grey_points_core[is_gtp == TRUE], 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = colored_non_gtp_core, 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = colored_gtp_core, 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type") +
  scale_shape_manual(name = "GTP Binding Pocket",
                     values = c("Non-GTP pocket" = 16, "GTP pocket" = 16),
                     labels = c("GTP pocket" = "Yes", "Non-GTP pocket" = "No")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
       title = "Core mutations outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_core, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
p5
ggsave(file.path(output_dir, "p5_core_mutations_outside_interfaces 2.pdf"), p5, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图5保存完成\n")

# ===============================
# 图6: Surface位点的突变（排除界面后）
# ===============================
cat("\n=== 绘制图6: Surface位点的突变（排除界面后） ===\n")

surface_mutations <- remaining_mutations[type == "surface"]
or_label_surface <- calculate_or_pvalue(surface_mutations)

grey_points_surface <- surface_mutations[allosteric_type_by_effect == "Other"]
colored_points_surface <- surface_mutations[allosteric_type_by_effect != "Other"]
colored_non_gtp_surface <- colored_points_surface[is_gtp == FALSE]
colored_gtp_surface <- colored_points_surface[is_gtp == TRUE]

p6 <- ggplot(surface_mutations, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(data = grey_points_surface[is_gtp == FALSE], 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = grey_points_surface[is_gtp == TRUE], 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = colored_non_gtp_surface, 
             aes(color = allosteric_type_by_effect, shape = "Non-GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  geom_point(data = colored_gtp_surface, 
             aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
             size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type") +
  scale_shape_manual(name = "GTP Binding Pocket",
                     values = c("Non-GTP pocket" = 16, "GTP pocket" = 16),
                     labels = c("GTP pocket" = "Yes", "Non-GTP pocket" = "No")) +
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
       title = "Surface mutations outside interfaces") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_surface, 
           hjust = 0, vjust = 1, size = 3.3, color = "black") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
p6
ggsave(file.path(output_dir, "p6_surface_mutations_outside_interfaces 2.pdf"), p6, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图6保存完成\n")

# ===============================
# 图7: GTP位点突变着色
# ===============================
cat("\n=== 绘制图7: GTP位点突变着色 ===\n")

data_plot_mutation[, is_gtp_pocket := Pos_real %in% GTP_Binding_pocket]

p7 <- ggplot(data_plot_mutation, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +
  geom_point(data = data_plot_mutation[is_gtp_pocket == TRUE], 
             color = "#09B636", size = 2.5, alpha = 0.8) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",
       title = "GTP pocket mutations") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
p7
ggsave(file.path(output_dir, "p7_gtp_pocket_mutations.pdf"), p7, 
       device = cairo_pdf, width = 4, height = 4)
cat("图7保存完成\n")

# ===============================
# 修正的OR计算函数 - 返回列表
# ===============================
calculate_or_pvalue <- function(data) {
  is_effect1 <- data$allosteric_mutation_RAF1 == TRUE
  is_effect2 <- data$allosteric_mutation_K13 == TRUE  # 改为K13
  
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
  cat("Total mutations:", nrow(data), "\n")
  cat("RAF1 allosteric effect:", sum(is_effect1), "\n")
  cat("K13 allosteric effect:", sum(is_effect2), "\n")  # 改为K13
  cat("Both effects:", sum(is_effect1 & is_effect2), "\n")
  cat("OR =", or_value, ",", p_label, "\n\n")
  
  # 确保返回列表
  return(list(or = or_value, p_label = p_label, contingency = contingency_table))
}

# ===============================
# 图8: 排除所有界面后的变构突变分析 - 修正版本
# ===============================
cat("\n=== 绘制图8: 排除所有界面后的变构突变分析 ===\n")

excluded_sites <- unique(c(RAF1_Binding_interface_site, K13_Binding_interface_site, GTP_Binding_pocket))  # 改为K13
remaining_mutations_clean <- data_plot_mutation[!Pos_real %in% excluded_sites]

# 分类变构突变类型 - 改为K13
remaining_mutations_clean[, allosteric_type := "Other"]
remaining_mutations_clean[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == FALSE, 
                          allosteric_type := "RAF1_allosteric_only"]
remaining_mutations_clean[allosteric_mutation_RAF1 == FALSE & allosteric_mutation_K13 == TRUE, 
                          allosteric_type := "K13_allosteric_only"]
remaining_mutations_clean[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == TRUE, 
                          allosteric_type := "Both_allosteric"]

# 计算OR值
or_result_clean <- calculate_or_pvalue(remaining_mutations_clean)

# 检查or_result_clean的类型
cat("or_result_clean类型:", class(or_result_clean), "\n")
cat("or_result_clean内容:\n")
print(or_result_clean)

# 安全地提取OR值和p_label
if (is.list(or_result_clean)) {
  or_label_clean <- paste0("OR = ", or_result_clean$or, "\n", or_result_clean$p_label)
} else {
  # 如果是字符串，直接使用
  or_label_clean <- or_result_clean
}


p8 <- ggplot(remaining_mutations_clean, 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +  # 改为K13
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),  # 改为K13
                     labels = c("RAF1 effect only", "K13 effect only", "Both effects", "Other")) +  # 改为K13
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",  # 改为K13
       title = "Allosteric mutations (exclude all interfaces)") +
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
ggsave(file.path(output_dir, "p8_allosteric_exclude_all_interfaces.pdf"), p8, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图8保存完成\n")

# ===============================
# 图9: Core位点变构突变（排除所有界面）- 修正版本
# ===============================
cat("\n=== 绘制图9: Core位点变构突变（排除所有界面） ===\n")

core_mutations_clean <- remaining_mutations_clean[type == "core"]
or_result_core_clean <- calculate_or_pvalue(core_mutations_clean)

# 安全地提取OR值和p_label
if (is.list(or_result_core_clean)) {
  or_label_core_clean <- paste0("OR = ", or_result_core_clean$or, "\n", or_result_core_clean$p_label)
} else {
  or_label_core_clean <- or_result_core_clean
}

p9 <- ggplot(core_mutations_clean, 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +  # 改为K13
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),  # 改为K13
                     labels = c("RAF1 effect only", "K13 effect only", "Both effects", "Other")) +  # 改为K13
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",  # 改为K13
       title = "Core allosteric mutations (exclude all interfaces)") +
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
ggsave(file.path(output_dir, "p9_core_allosteric_exclude_all_interfaces.pdf"), p9, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图9保存完成\n")

# ===============================
# 图10: Surface位点变构突变（排除所有界面）- 修正版本
# ===============================
cat("\n=== 绘制图10: Surface位点变构突变（排除所有界面） ===\n")

surface_mutations_clean <- remaining_mutations_clean[type == "surface"]
or_result_surface_clean <- calculate_or_pvalue(surface_mutations_clean)

# 安全地提取OR值和p_label
if (is.list(or_result_surface_clean)) {
  or_label_surface_clean <- paste0("OR = ", or_result_surface_clean$or, "\n", or_result_surface_clean$p_label)
} else {
  or_label_surface_clean <- or_result_surface_clean
}

p10 <- ggplot(surface_mutations_clean, 
              aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +  # 改为K13
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),  # 改为K13
                     labels = c("RAF1 effect only", "K13 effect only", "Both effects", "Other")) +  # 改为K13
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",  # 改为K13
       title = "Surface allosteric mutations (exclude all interfaces)") +
  annotate("text", x = xlim_effect[1] + 0.5, y = ylim_effect[2] - 0.1, label = or_label_surface_clean, 
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
p10
ggsave(file.path(output_dir, "p10_surface_allosteric_exclude_all_interfaces.pdf"), p10, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图10保存完成\n")

# ===============================
# 图11: GTP口袋位点变构突变 - 修正版本
# ===============================
cat("\n=== 绘制图11: GTP口袋位点变构突变 ===\n")

gtp_mutations <- data_plot_mutation[Pos_real %in% GTP_Binding_pocket]

# 分类GTP口袋的变构突变类型 - 改为K13
gtp_mutations[, allosteric_type := "Other"]
gtp_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == FALSE, 
              allosteric_type := "RAF1_allosteric_only"]
gtp_mutations[allosteric_mutation_RAF1 == FALSE & allosteric_mutation_K13 == TRUE, 
              allosteric_type := "K13_allosteric_only"]
gtp_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K13 == TRUE, 
              allosteric_type := "Both_allosteric"]

or_result_gtp <- calculate_or_pvalue(gtp_mutations)

# 安全地提取OR值和p_label
if (is.list(or_result_gtp)) {
  or_label_gtp <- paste0("OR = ", or_result_gtp$or, "\n", or_result_gtp$p_label)
} else {
  or_label_gtp <- or_result_gtp
}

p11 <- ggplot(gtp_mutations, 
              aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K13`)) +  # 改为K13
  geom_point(aes(color = allosteric_type), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = allosteric_colors,
                     name = "Allosteric Effect Type",
                     breaks = c("RAF1_allosteric_only", "K13_allosteric_only", "Both_allosteric", "Other"),  # 改为K13
                     labels = c("RAF1 effect only", "K13 effect only", "Both effects", "Other")) +  # 改为K13
  coord_cartesian(xlim = xlim_effect, ylim = ylim_effect) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K13 (kcal/mol)",  # 改为K13
       title = "GTP pocket allosteric mutations") +
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
ggsave(file.path(output_dir, "p11_gtp_pocket_allosteric.pdf"), p11, 
       device = cairo_pdf, width = 4, height = 4.5)
cat("图11保存完成\n")

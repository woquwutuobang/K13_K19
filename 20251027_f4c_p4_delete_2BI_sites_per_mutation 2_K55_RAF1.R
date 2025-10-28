library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)
library(dplyr)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 函数1：读取和预处理ddG数据
read_ddG_data <- function(input, wt_aa) {
  ddG <- fread(input)
  ddG[, `:=`(Pos_real = Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon = substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon = substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt = paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:(nchar(wt_aa) + 1), each = 20),
    mt_codon = rep(aa_list, times = nchar(wt_aa))
  )
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  return(ddG)
}

# 函数2：计算加权平均ddG
calculate_weighted_mean_ddG <- function(ddG) {
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # 获取野生型氨基酸
  wt_codon_info <- ddG[Pos_real > 1, .(codon = first(wt_codon)), by = "Pos_real"]
  
  # 合并信息
  final_output <- merge(wt_codon_info, weighted_mean_ddG, by = "Pos_real")
  
  return(final_output)
}




# 函数3：判断site_type和mutation_type
classify_site_mutation_types <- function(ddG, weighted_mean_ddG, anno, assay_sele, reg_threshold = NULL) {
  
  # 合并数据
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  # 定义结合位点类型
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # 计算调节阈值（如果没有提供）
  if (is.null(reg_threshold)) {
    reg_threshold <- data_plot[binding_type == "binding site",
                               sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                               .SDcols = c("mean", "sigma")]
    
    # 打印reg_threshold
    cat("Calculated regulatory threshold for", assay_sele, ":", reg_threshold, "\n")
  } else {
    cat("Using provided regulatory threshold for", assay_sele, ":", reg_threshold, "\n")
  }
  
  # 分类位点类型
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # 合并突变数据与位点类型信息
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # 识别变构突变
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # 分类突变类型
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
  
  # 设置突变类型的因子水平
  data_plot_mutation <- within(data_plot_mutation,
                               mutation_type <- factor(mutation_type,
                                                       levels = c("Orthosteric site huge differences",
                                                                  "Orthosteric site small differences",
                                                                  "GTP binding allosteric mutation",
                                                                  "GTP binding other mutation",
                                                                  "Allosteric mutation",
                                                                  "Other mutation")))
  
  # 返回结果列表
  #  result <- list(
  #    site_data = data_plot,
  #    mutation_data = data_plot_mutation,
  #    regulatory_threshold = reg_threshold
  #  )
  
  return(data_plot_mutation)
}



ddG_RAF1<-read_ddG_data("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",wt_aa)
ddG_K55<-read_ddG_data("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",wt_aa)


ddG_RAF1_weighted<-calculate_weighted_mean_ddG(ddG_RAF1)
ddG_K55_weighted<-calculate_weighted_mean_ddG(ddG_K55)


anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv")



result_RAF1 <- classify_site_mutation_types(
  ddG = ddG_RAF1,                    
  weighted_mean_ddG = ddG_RAF1_weighted, 
  anno = anno,                
  assay_sele = "RAF1"                    
)



result_K55 <- classify_site_mutation_types(
  ddG = ddG_K55,                    
  weighted_mean_ddG = ddG_K55_weighted, 
  anno = anno,                
  assay_sele = "K55"                    
)

#names(result_RAF1)
#names(result_K55)
result_RAF1<-result_RAF1[,c(2:4,23:29)]
result_K55<-result_K55[,c(2:4,23:29)]
#mutation_data <- result_RAF1$mutation_data    
#reg_threshold <- result_RAF1$regulatory_threshold 
#table(mutation_data$mutation_type)

data_plot_mutation1 <- merge(result_RAF1, result_K55, by = c("Pos_real","wt_codon","mt_codon"), suffixes = c("_RAF1", "_K55"),all = F)


core_surface <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv")
core_surface <- core_surface[, Pos_real := Pos]
core_surface_threshold <- 0.25
#surface_threshold_sasa <- 0.26
core_surface <- core_surface %>%
  mutate(type = case_when(
    RASA <= core_surface_threshold ~ "core",
    RASA > core_surface_threshold ~ "surface"
  )) %>%
  filter(!is.na(type))


data_plot_mutation<-merge(data_plot_mutation1, core_surface, by = "Pos_real",all=T)

names(data_plot_mutation)

K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
GTP_Binding_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)



# 直接使用data_plot_mutation数据框绘制突变级别的图
### p1 - 所有突变
p1 <- ggplot(data_plot_mutation, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K55`)) +
  geom_point(color = "grey70", size = 2.5, alpha = 0.6) +
  #geom_errorbar(aes(ymin = `mean_kcal/mol_K55` - `std_kcal/mol_K55`, 
  #                  ymax = `mean_kcal/mol_K55` + `std_kcal/mol_K55`), 
  #              width = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  #geom_errorbarh(aes(xmin = `mean_kcal/mol_RAF1` - `std_kcal/mol_RAF1`, 
  #                   xmax = `mean_kcal/mol_RAF1` + `std_kcal/mol_RAF1`), 
  #               height = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(-1.3, 3), ylim = c(-1.7, 2.6)) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K55 (kcal/mol)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p1)

### p2 - K55结合界面位点的突变
p2 <- ggplot(data_plot_mutation[Pos_real %in% K55_Binding_interface_site], 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K55`)) +
  geom_point(color = "#F4270C", size = 2.5, alpha = 0.6) +
  #geom_errorbar(aes(ymin = `mean_kcal/mol_K55` - `std_kcal/mol_K55`, 
  #                  ymax = `mean_kcal/mol_K55` + `std_kcal/mol_K55`), 
  #              width = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  #geom_errorbarh(aes(xmin = `mean_kcal/mol_RAF1` - `std_kcal/mol_RAF1`, 
  #                   xmax = `mean_kcal/mol_RAF1` + `std_kcal/mol_RAF1`), 
  #               height = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(-1.3, 3), ylim = c(-1.7, 2.6)) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K55 (kcal/mol)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p2)

### p3 - RAF1结合界面位点的突变
p3 <- ggplot(data_plot_mutation[Pos_real %in% RAF1_Binding_interface_site], 
             aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K55`)) +
  geom_point(color = "#1B38A6", size = 2.5, alpha = 0.6) +
  #geom_errorbar(aes(ymin = `mean_kcal/mol_K55` - `std_kcal/mol_K55`, 
  #                  ymax = `mean_kcal/mol_K55` + `std_kcal/mol_K55`), 
  #              width = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  #geom_errorbarh(aes(xmin = `mean_kcal/mol_RAF1` - `std_kcal/mol_RAF1`, 
  #                   xmax = `mean_kcal/mol_RAF1` + `std_kcal/mol_RAF1`), 
  #               height = 0, color = "grey30", size = 0.1, alpha = 0.3) +
  coord_cartesian(xlim = c(-1.3, 3), ylim = c(-1.7, 2.6)) +
  labs(x = "ΔΔGb RAF1 (kcal/mol)",
       y = "ΔΔGb K55 (kcal/mol)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p3)

# 组合图形
combined_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3)

print(combined_plot)

# 保存组合图
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251028/combined_mutation_scatter_plots_p1_p2_p3 RAF1 VS K55 2.pdf", 
       combined_plot, device = cairo_pdf, width = 15, height = 5)




### p4/p5/p6

### 移除两个结合界面的位点的突变
remaining_mutations <- data_plot_mutation[!Pos_real %in% c(K55_Binding_interface_site, RAF1_Binding_interface_site)]

# 添加GTP结合口袋标记
remaining_mutations[, is_gtp := Pos_real %in% GTP_Binding_pocket]

# 基于突变的功能效应（allosteric_mutation列）来判断
remaining_mutations[, allosteric_type_by_effect := "Other"]
remaining_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K55 == FALSE, 
                    allosteric_type_by_effect := "RAF1_allosteric_only"]  # 只对RAF1有变构效应
remaining_mutations[allosteric_mutation_RAF1 == FALSE & allosteric_mutation_K55 == TRUE, 
                    allosteric_type_by_effect := "K55_allosteric_only"]  # 只对K55有变构效应
remaining_mutations[allosteric_mutation_RAF1 == TRUE & allosteric_mutation_K55 == TRUE, 
                    allosteric_type_by_effect := "Both_allosteric"]  # 对两者都有变构效应

# 创建颜色映射
allosteric_colors <- c(
  "Other" = "grey70",
  "K55_allosteric_only" = "#FFB0A5",  
  "RAF1_allosteric_only" = "#F4AD0C",  
  "Both_allosteric" = "#C68EFD"        
)

# 函数：基于突变的功能效应计算OR值和创建标签文本
calculate_mutation_or_label_by_effect <- function(mutation_data) {
  # 计算每个突变是否对K55/RAF1有变构效应
  is_k55_effect <- mutation_data$allosteric_mutation_K55 == TRUE
  is_raf1_effect <- mutation_data$allosteric_mutation_RAF1 == TRUE
  
  contingency_table <- table(is_k55_effect, is_raf1_effect)
  fisher_test <- fisher.test(contingency_table)
  
  or_value <- round(fisher_test$estimate, 2)
  if (fisher_test$p.value < 0.001) {
    p_label <- "p < 0.001"
  } else if (fisher_test$p.value < 0.01) {
    p_label <- "p < 0.01" 
  } else if (fisher_test$p.value < 0.05) {
    p_label <- "p < 0.05"
  } else {
    p_label <- paste0("p = ", round(fisher_test$p.value, 3))
  }
  
  # 打印详细的列联表信息
  cat("Contingency Table (Functional Effects):\n")
  print(contingency_table)
  cat("Total mutations:", nrow(mutation_data), "\n")
  cat("Mutations with K55 allosteric effect:", sum(is_k55_effect), "\n")
  cat("Mutations with RAF1 allosteric effect:", sum(is_raf1_effect), "\n")
  cat("Mutations with both effects:", sum(is_k55_effect & is_raf1_effect), "\n")
  cat("Proportion with K55 effect:", round(mean(is_k55_effect) * 100, 1), "%\n")
  cat("Proportion with RAF1 effect:", round(mean(is_raf1_effect) * 100, 1), "%\n")
  cat("Proportion with both effects:", round(mean(is_k55_effect & is_raf1_effect) * 100, 1), "%\n")
  
  return(paste0("OR = ", or_value, "\n", p_label))
}

# 函数：创建突变级别的散点图（基于功能效应，灰色点在底层）
create_mutation_scatter_plot_by_effect <- function(data, title, or_label) {
  
  # 将数据分为灰色点和其他颜色点 - 使用allosteric_type_by_effect
  grey_points <- data[allosteric_type_by_effect == "Other"]
  colored_points <- data[allosteric_type_by_effect != "Other"]
  
  # 进一步将彩色点分为GTP和非GTP
  colored_non_gtp <- colored_points[is_gtp == FALSE]
  colored_gtp <- colored_points[is_gtp == TRUE]
  
  ggplot(data, aes(x = `mean_kcal/mol_RAF1`, y = `mean_kcal/mol_K55`)) +
    # 第一层：灰色非GTP点（圆形）- 最底层
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
    # 第四层：彩色GTP点（三角形）- 最顶层
    geom_point(data = colored_gtp, 
               aes(color = allosteric_type_by_effect, shape = "GTP pocket"), 
               size = 2.5, alpha = 0.6) +
    # 颜色标度
    scale_color_manual(values = allosteric_colors,
                       name = "Allosteric Effect Type",
                       breaks = c("K55_allosteric_only", "RAF1_allosteric_only", "Both_allosteric", "Other"),
                       labels = c("K55 effect only", "RAF1 effect only", "Both effects", "Other")) +
    # 形状标度
    scale_shape_manual(name = "GTP Binding Pocket",
                       values = c("Non-GTP pocket" = 16, "GTP pocket" = 17),
                       labels = c("GTP pocket" = "Yes", "Non-GTP pocket" = "No")) +
    coord_cartesian(xlim = c(-1.5, 3), ylim = c(-1.7, 2.6)) +
    labs(x = "ΔΔGb RAF1 (kcal/mol)",
         y = "ΔΔGb K55 (kcal/mol)",
         title = title) +
    # 添加OR和p值标签
    annotate("text", x = -1, y = 2.5, label = or_label, 
             hjust = 0, vjust = 1, size = 3.3, color = "black") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}

# 计算整体的OR值（用于p4）- 基于功能效应
cat("\n=== CALCULATING OR BASED ON FUNCTIONAL EFFECTS ===\n")
or_label_all <- calculate_mutation_or_label_by_effect(remaining_mutations)

# p4: 所有突变
p4 <- create_mutation_scatter_plot_by_effect(remaining_mutations, "All mutations outside interfaces", or_label_all)

# p5: Core位点的突变
core_mutations <- remaining_mutations[type == "core"]
or_label_core <- calculate_mutation_or_label_by_effect(core_mutations)
p5 <- create_mutation_scatter_plot_by_effect(core_mutations, "Core mutations outside interfaces", or_label_core)

# p6: Surface位点的突变
surface_mutations <- remaining_mutations[type == "surface"]
or_label_surface <- calculate_mutation_or_label_by_effect(surface_mutations)
p6 <- create_mutation_scatter_plot_by_effect(surface_mutations, "Surface mutations outside interfaces", or_label_surface)

# 打印各个OR值
cat("\n=== SUMMARY ===\n")
cat("All mutations OR:", or_label_all, "\n")
cat("Core mutations OR:", or_label_core, "\n")
cat("Surface mutations OR:", or_label_surface, "\n")

# 打印各个组的突变数量和功能效应分布
cat("\n=== MUTATION COUNTS AND EFFECT DISTRIBUTION ===\n")
cat("All mutations count:", nrow(remaining_mutations), "\n")
print(table(remaining_mutations$allosteric_type_by_effect))

if(exists("core_mutations")) {
  cat("\nCore mutations count:", nrow(core_mutations), "\n")
  print(table(core_mutations$allosteric_type_by_effect))
}

if(exists("surface_mutations")) {
  cat("\nSurface mutations count:", nrow(surface_mutations), "\n")
  print(table(surface_mutations$allosteric_type_by_effect))
}

print(p4)
print(p5)
print(p6)

# 组合图形
combined_plot <- p4 + p5 + p6 + 
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

# 保存组合图
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251028/core_surface_mutation_scatter_plots_by_effect_p4_p5_p6 RAF1 VS K55.pdf", 
       combined_plot, device = cairo_pdf, width = 15, height = 5.5)

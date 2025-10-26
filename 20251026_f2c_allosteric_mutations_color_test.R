library(data.table)
library(ggplot2)
library(dplyr)
library(krasddpcams)

# p值计算函数
krasddpcams__pvalue <- function(x, std) {
  2 * pnorm(-abs(x/std))
}

# 判断mutation_type的函数
classify_mutation_type <- function(input, assay_sele, anno, wt_aa) {
  # Read and process ddG data
  ddG <- fread(input)
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  anno_data <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  # Define binding site types
  data_plot[, binding_type := "allosteric site"]
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # Merge mutation data with site type information
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(`mean_kcal/mol`) - reg_threshold, `std_kcal/mol`),
    method = "BH") < 0.05 & (abs(`mean_kcal/mol`) - reg_threshold) > 0]
  
  # Classify mutation types
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
  
  # Set factor levels for mutation types
  data_plot_mutation <- within(data_plot_mutation,
                               mutation_type <- factor(mutation_type,
                                                       levels = c("Orthosteric site huge differences",
                                                                  "Orthosteric site small differences",
                                                                  "GTP binding allosteric mutation",
                                                                  "GTP binding other mutation",
                                                                  "Allosteric mutation",
                                                                  "Other mutation")))
  
  # 添加assay标识
  data_plot_mutation[, assay := assay_sele]
  
  return(data_plot_mutation)
}




###=================================================
# 绘制两个assay相关性散点图的函数
plot_assay_correlation <- function(input1, input2, assay1, assay2, anno, wt_aa, 
                                   point_size = 2, alpha = 0.7, base_size = 10,
                                   xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3)) {
  
  # 获取两个assay的mutation type数据
  data_assay1 <- classify_mutation_type(input1, assay1, anno, wt_aa)
  data_assay2 <- classify_mutation_type(input2, assay2, anno, wt_aa)
  
  # 选择需要的列并重命名
  data1_clean <- data_assay1[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type)]
  setnames(data1_clean, "mean_kcal/mol", paste0("mean_", assay1))
  setnames(data1_clean, "mutation_type", paste0("mutation_type_", assay1))
  
  data2_clean <- data_assay2[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type)]
  setnames(data2_clean, "mean_kcal/mol", paste0("mean_", assay2))
  setnames(data2_clean, "mutation_type", paste0("mutation_type_", assay2))
  
  # 合并两个assay的数据
  merged_data <- merge(data1_clean, data2_clean, by = c("mt", "Pos_real"), all = TRUE)
  
  # 定义颜色分组 - 使用动态的assay名称
  merged_data[, color_group := "Other"]
  
  # 设置颜色分组
  merged_data[get(paste0("mutation_type_", assay1)) %in% 
                c("Allosteric mutation", "GTP binding allosteric mutation") &
                !(get(paste0("mutation_type_", assay2)) %in% 
                    c("Allosteric mutation", "GTP binding allosteric mutation")), 
              color_group := paste0(assay1, " allosteric")]
  
  merged_data[get(paste0("mutation_type_", assay2)) %in% 
                c("Allosteric mutation", "GTP binding allosteric mutation") &
                !(get(paste0("mutation_type_", assay1)) %in% 
                    c("Allosteric mutation", "GTP binding allosteric mutation")), 
              color_group := paste0(assay2, " allosteric")]
  
  merged_data[get(paste0("mutation_type_", assay1)) %in% 
                c("Allosteric mutation", "GTP binding allosteric mutation") &
                get(paste0("mutation_type_", assay2)) %in% 
                c("Allosteric mutation", "GTP binding allosteric mutation"), 
              color_group := "Overlap allosteric"]
  
  # 动态设置颜色映射
  color_levels <- c("Other", 
                    paste0(assay1, " allosteric"), 
                    paste0(assay2, " allosteric"), 
                    "Overlap allosteric")
  
  # 动态颜色映射 - 使用固定的颜色但动态的名称
  color_map <- setNames(
    c("grey", "#F4AD0C", "#F1DD10", "#C68EFD"),
    color_levels
  )
  
  # 设置因子水平以确保正确的绘图顺序
  merged_data[, color_group := factor(color_group, levels = color_levels)]
  
  # 按颜色分组排序数据，确保灰色(Other)在最底层
  merged_data <- merged_data[order(as.numeric(color_group))]
  
  # 创建散点图
  p <- ggplot(merged_data, aes(x = .data[[paste0("mean_", assay1)]], 
                               y = .data[[paste0("mean_", assay2)]], 
                               color = color_group)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = color_map, name = "Allosteric Mutations",
                       breaks = color_levels[color_levels != "Other"]) +  # 图例中不显示Other
    labs(x = bquote(Delta*Delta*"G ("*.(assay1)*") (kcal/mol)"), 
         y = bquote(Delta*Delta*"G ("*.(assay2)*") (kcal/mol)")) +
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
  
  # 打印调试信息
  print(paste("Total points:", nrow(merged_data)))
  print(table(merged_data$color_group))
  
  return(list(plot = p, data = merged_data))
}

# 使用示例
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"



########=====
# 比较K13和RAF1
result <- plot_assay_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K13",
  assay2 = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K13_vs_RAF1_allosteric_correlation.pdf", 
       plot = result$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K13_vs_RAF1_allosteric_correlation.png", 
#       plot = result$plot, 
#      width = 4, 
#       height = 4.3)





########=====
# 比较K13和K19
result <- plot_assay_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K13",
  assay2 = "K19",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K13_vs_K19_allosteric_correlation.pdf", 
       plot = result$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K13_vs_RAF1_allosteric_correlation.png", 
#       plot = result$plot, 
#      width = 4, 
#       height = 4.3)




########=====
# 比较K27和RAF1
result <- plot_assay_correlation(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K27",
  assay2 = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K27_vs_RAF1_allosteric_correlation.pdf", 
       plot = result$plot, 
       width = 4, 
       height = 4.3, 
       device = cairo_pdf)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/K13_vs_RAF1_allosteric_correlation.png", 
#       plot = result$plot, 
#      width = 4, 
#       height = 4.3)

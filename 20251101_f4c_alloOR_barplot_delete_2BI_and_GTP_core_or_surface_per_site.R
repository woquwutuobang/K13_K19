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

# 定义变构位点
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146,10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)
RALGDS_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 117, 146,10, 20, 58, 59, 85)
PI3KCG_allosteric_site <- c(16, 17, 18, 28, 32, 34, 35, 57, 60, 117, 146,20, 55, 58, 59, 68, 85)
SOS1_allosteric_site <- c(16, 17, 18, 28, 32, 35, 57, 117, 146,10, 40, 54, 55, 58, 63, 68, 69, 85, 144, 148)
K55_allosteric_site <- c(15, 16, 17, 28, 32, 35, 57, 60, 117,10, 20, 58, 59, 68, 69, 71, 72, 78, 79, 85)
K27_allosteric_site <- c(16, 17, 28, 57, 60, 117, 119,10, 63, 76, 84, 90, 144, 151, 155, 156)
K13_allosteric_site <- c(15, 16, 17, 35, 145,10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
K19_allosteric_site <- c(15, 16, 17, 145,8, 10, 19, 21, 55, 77, 78, 79, 82, 93, 151, 159, 163)

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

# ===============================
# 辅助函数
# ===============================
create_empty_result <- function(assay1, assay2, region_type, status) {
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
    OR = NA_real_,
    p_value = NA_real_,
    n_sites = 0,
    n_allo1 = 0,
    n_allo2 = 0,
    n_both = 0,
    sites_removed = NA_integer_,
    status = status
  ))
}

# ===============================
# OR计算函数
# ===============================
calculate_pair_or_pvalue_site_level <- function(assay1, assay2, region_type) {
  cat("处理", assay1, "vs", assay2, "-", region_type, "...\n")
  
  tryCatch({
    # 分别读取两个assay的数据
    input1 <- get(paste0("input_", assay1))
    input2 <- get(paste0("input_", assay2))
    
    # 处理assay1数据
    ddG_data1 <- read_ddG_data(input1, wt_aa)
    
    # 处理assay2数据
    ddG_data2 <- read_ddG_data(input2, wt_aa)
    
    # 合并这两个assay的数据（基于位点）
    data1 <- unique(ddG_data1[Pos_real > 1, .(Pos_real, wt_codon)])
    data2 <- unique(ddG_data2[Pos_real > 1, .(Pos_real, wt_codon)])
    
    site_data <- merge(data1, data2, by = c("Pos_real", "wt_codon"), all = FALSE)
    
    # 检查合并后是否有数据
    if (nrow(site_data) == 0) {
      warning(paste("No overlapping sites between", assay1, "and", assay2))
      return(create_empty_result(assay1, assay2, region_type, "No overlapping sites"))
    }
    
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
    
    site_data <- merge(site_data, core_surface[, .(Pos_real, type)], by = "Pos_real", all.x = TRUE)
    
    # 添加变构位点信息
    site_data[, paste0("is_allosteric_", assay1) := Pos_real %in% allosteric_sites[[assay1]]]
    site_data[, paste0("is_allosteric_", assay2) := Pos_real %in% allosteric_sites[[assay2]]]
    
    # 移除两个assay各自的结合界面位点和GTP pocket
    sites_to_remove <- unique(c(binding_sites[[assay1]], binding_sites[[assay2]], GTP_Binding_pocket_site))
    
    # 根据区域类型筛选数据
    if (region_type == "core") {
      filtered_data <- site_data[!Pos_real %in% sites_to_remove & type == "core"]
    } else if (region_type == "surface") {
      filtered_data <- site_data[!Pos_real %in% sites_to_remove & type == "surface"]
    }
    
    # 检查过滤后是否有数据
    if (nrow(filtered_data) == 0) {
      warning(paste("No data remaining after filtering for", assay1, "vs", assay2, "-", region_type))
      return(create_empty_result(assay1, assay2, region_type, "No data after filtering"))
    }
    
    # 计算OR值（基于位点）
    is_allo1 <- filtered_data[[paste0("is_allosteric_", assay1)]]
    is_allo2 <- filtered_data[[paste0("is_allosteric_", assay2)]]
    
    # 检查是否有TRUE值
    if (sum(is_allo1, na.rm = TRUE) == 0 || sum(is_allo2, na.rm = TRUE) == 0) {
      warning(paste("No allosteric sites found for", assay1, "or", assay2, "-", region_type))
      return(create_empty_result(assay1, assay2, region_type, "No allosteric sites"))
    }
    
    # 创建列联表
    contingency_table <- table(is_allo1, is_allo2)
    
    # 检查列联表维度
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
      warning(paste("Contingency table too small for", assay1, "vs", assay2, "-", region_type))
      return(create_empty_result(assay1, assay2, region_type, "Invalid contingency table"))
    }
    
    # 检查是否所有值都在同一行或列（无法计算OR）
    if (any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) {
      warning(paste("All values in one row/column for", assay1, "vs", assay2, "-", region_type))
      return(create_empty_result(assay1, assay2, region_type, "All values in one category"))
    }
    
    # 执行Fisher检验
    fisher_test <- fisher.test(contingency_table)
    
    or_value <- fisher_test$estimate
    p_value <- fisher_test$p.value
    
    # 检查是否为Inf值，如果是则丢弃
    if (is.infinite(or_value)) {
      warning(paste("OR is infinite for", assay1, "vs", assay2, "-", region_type, "- discarding"))
      return(create_empty_result(assay1, assay2, region_type, "OR is infinite"))
    }
    
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
      n_sites = nrow(filtered_data),
      n_allo1 = sum(is_allo1, na.rm = TRUE),
      n_allo2 = sum(is_allo2, na.rm = TRUE),
      n_both = sum(is_allo1 & is_allo2, na.rm = TRUE),
      sites_removed = length(sites_to_remove),
      status = "Success"
    ))
    
  }, error = function(e) {
    warning(paste("Error calculating OR for", assay1, "vs", assay2, "-", region_type, ":", e$message))
    return(create_empty_result(assay1, assay2, region_type, paste("Error:", e$message)))
  })
}

# ===============================
# 绘图函数
# ===============================
plot_or_barplot <- function(or_data, title_suffix) {
  # 只保留成功的分析
  plot_data <- or_data[status == "Success"]
  
  if (nrow(plot_data) == 0) {
    cat("No successful analyses to plot for", title_suffix, "\n")
    return(NULL)
  }
  
  cat("绘制", title_suffix, "，数据行数:", nrow(plot_data), "\n")
  print(plot_data[, .(assay1, assay2, OR, p_value)])
  
  # 设置颜色
  comparison_colors <- c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C", 
    "BI2_vs_BI2" = "#1B38A6"
  )
  
  # 创建组合标签
  plot_data[, combination_label := paste0(assay1, " vs ", assay2)]
  
  # 按照比较类型分区，然后在每个分区内按OR值排序
  plot_data[, comparison_type := factor(comparison_type, 
                                        levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2"))]
  
  # 在每个比较类型内按OR值排序
  plot_data[, within_group_order := frank(OR, ties.method = "first"), by = comparison_type]
  plot_data[, final_order := as.numeric(comparison_type) * 100 + within_group_order]
  
  # 重新排序数据
  plot_data <- plot_data[order(final_order)]
  plot_data[, combination_label := factor(combination_label, levels = unique(combination_label[order(final_order)]))]
  
  # 计算标注位置
  max_or <- max(plot_data$OR, na.rm = TRUE)
  min_or <- min(plot_data$OR, na.rm = TRUE)
  
  # 绘制柱状图
  p <- ggplot(plot_data, aes(x = combination_label, y = OR, fill = comparison_type)) +
    geom_col(width = 0.7) +
    # 添加OR值标注
    geom_text(aes(label = sprintf("OR=%.2f", OR), 
                  y = OR + max_or * 0.05), 
              vjust = 0, size = 2.5, color = "black") +
    # 添加显著性标记
    geom_text(aes(label = significance, 
                  y = OR + max_or * 0.12), 
              vjust = 0, size = 3, color = "#F1DD10") +
    scale_fill_manual(values = comparison_colors,
                      name = "Comparison Type",
                      labels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")) +
    labs(x = "Assay Combination",
         y = "Odds Ratio (OR)",
         title = paste("Odds Ratios for Allosteric Site Co-occurrence(delete 2BIs and GTP pocket)\n per site -", title_suffix)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom"
    ) +
    coord_cartesian(ylim = c(0, max_or * 1.25))
  
  return(p)
}

# ===============================
# 主分析流程
# ===============================
cat("=== 开始完整位点级别OR分析 ===\n")

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
  result_core <- calculate_pair_or_pvalue_site_level(assay1, assay2, "core")
  or_results_core <- rbind(or_results_core, result_core)
  
  # 计算surface区域
  result_surface <- calculate_pair_or_pvalue_site_level(assay1, assay2, "surface")
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


# 只对成功的分析添加显著性标记
or_results_core_success <- or_results_core[status == "Success"]
or_results_surface_success <- or_results_surface[status == "Success"]

if (nrow(or_results_core_success) > 0) {
  or_results_core_success <- add_significance(or_results_core_success)
  # 为失败的数据添加空的significance列
  or_results_core_failure <- or_results_core[status != "Success"]
  if (nrow(or_results_core_failure) > 0) {
    or_results_core_failure[, significance := NA_character_]
  }
  # 合并
  or_results_core <- rbind(or_results_core_success, or_results_core_failure, fill = TRUE)
} else {
  or_results_core[, significance := NA_character_]
}

if (nrow(or_results_surface_success) > 0) {
  or_results_surface_success <- add_significance(or_results_surface_success)
  # 为失败的数据添加空的significance列
  or_results_surface_failure <- or_results_surface[status != "Success"]
  if (nrow(or_results_surface_failure) > 0) {
    or_results_surface_failure[, significance := NA_character_]
  }
  # 合并
  or_results_surface <- rbind(or_results_surface_success, or_results_surface_failure, fill = TRUE)
} else {
  or_results_surface[, significance := NA_character_]
}
# 合并回原表
or_results_core <- rbind(or_results_core_success, or_results_core[status != "Success"])
or_results_surface <- rbind(or_results_surface_success, or_results_surface[status != "Success"])

# 检查结果
cat("Core区域成功分析数:", nrow(or_results_core[status == "Success"]), "/", nrow(or_results_core), "\n")
cat("Surface区域成功分析数:", nrow(or_results_surface[status == "Success"]), "/", nrow(or_results_surface), "\n")

# 打印失败的分析
failed_core <- or_results_core[status != "Success"]
failed_surface <- or_results_surface[status != "Success"]

if (nrow(failed_core) > 0) {
  cat("\nCore区域失败的分析:\n")
  print(failed_core[, .(assay1, assay2, status)])
}

if (nrow(failed_surface) > 0) {
  cat("\nSurface区域失败的分析:\n")
  print(failed_surface[, .(assay1, assay2, status)])
}

# 绘制并保存柱状图
cat("\n=== 绘制并保存柱状图 ===\n")

# 绘制core区域的图
p_core <- plot_or_barplot(or_results_core, "Core Regions")
if (!is.null(p_core)) {
  print(p_core)
  ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf5/20251101/OR_comparison_core_regions_delete_2BIs_and_GTP_pocket_per_site.pdf", 
         p_core, device = cairo_pdf, width = 12, height = 8)
  cat("Core区域图已保存\n")
} else {
  cat("No core region results to plot\n")
}

# 绘制surface区域的图
p_surface <- plot_or_barplot(or_results_surface, "Surface Regions")
if (!is.null(p_surface)) {
  print(p_surface)
  ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf5/20251101/OR_comparison_surface_regions_delete_2BIs_and_GTP_pocket_per_site.pdf", 
         p_surface, device = cairo_pdf, width = 12, height = 8)
  cat("Surface区域图已保存\n")
} else {
  cat("No surface region results to plot\n")
}






# 保存结果到文件
fwrite(or_results_core, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/OR_results_core_regions_site_level.csv")
fwrite(or_results_surface, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/OR_results_surface_regions_site_level.csv")

cat("分析完成！\n")
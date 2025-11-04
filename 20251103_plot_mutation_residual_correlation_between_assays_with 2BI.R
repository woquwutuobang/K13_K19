library(ggplot2)
library(data.table)
library(ggpubr)
library(dplyr)

# ===============================
# 函数：距离效应拟合和残差计算
# ===============================
calculate_distance_residuals <- function(input, assay_sele, anno_file) {
  # 1. 读取数据
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]  # 处理Pos
  data <- data[, c(1, 20:23)]
  colnames(data)[2:4] <- paste0(colnames(data)[2:4], "_", assay_sele)  # 修改列名后缀
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  # 2. 设置参数 - 根据assay自适应选择距离列
  x_cols <- paste0("scHAmin_ligand_", assay_sele)
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  cat("当前assay:", assay_sele, "\n")
  cat("使用的距离列:", x_cols, "\n")
  cat("使用的能量列:", y_col, "\n")
  
  # 3. 存储结果的列表
  result_list <- list()
  
  # 4. 对距离列进行拟合和残差计算
  for (i in seq_along(x_cols)) {
    cat("正在处理:", x_cols[i], "\n")
    
    # 提取数据
    xvector <- anno_final[[x_cols[i]]]
    yvector <- abs(anno_final[[y_col]])  # 取能量的绝对值
    
    df <- data.frame(
      Pos_real = anno_final$Pos_real,
      Pos = anno_final$Pos,
      x = xvector, 
      y = yvector,
      original_y = anno_final[[y_col]]  # 保留原始能量值（有正负号）
    )
    df <- df[complete.cases(df), ]  # 去除缺失值
    
    cat("有效数据点:", nrow(df), "\n")
    
    # --- 初始参数优化 ---
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]; b = params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess <- c(a = 1, b = -0.1)
    opt_params <- tryCatch(
      optim(initial_guess, residual_sum_of_squares, x = df$x, y = df$y)$par,
      error = function(e) {
        cat("初始参数优化失败，使用默认参数\n")
        c(a = 1, b = -0.1)
      }
    )
    
    cat("优化后的初始参数: a =", round(opt_params[1], 3), "b =", round(opt_params[2], 3), "\n")
    
    # --- 拟合指数模型 ---
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), data = df, start = list(a = opt_params[1], b = opt_params[2])),
      error = function(e) {
        cat("模型拟合失败:", e$message, "\n")
        NULL
      }
    )
    
    # --- 计算拟合值和残差 ---
    if (!is.null(fit_model)) {
      # 计算每个数据点的拟合值
      df$fitted_value <- predict(fit_model, newdata = data.frame(x = df$x))
      
      # 计算残差：原始能量值 - 拟合值
      df$residual <- df$original_y - df$fitted_value
      
      # 获取拟合参数
      fit_summary <- summary(fit_model)
      coefs <- fit_summary$coefficients
      a_val <- round(coefs["a", "Estimate"], 3)
      b_val <- round(coefs["b", "Estimate"], 3)
      p_val_b <- coefs["b", "Pr(>|t|)"]
      
      p_text <- if (is.na(p_val_b)) {
        "p = NA"
      } else if (p_val_b < 0.001) {
        "p < 0.001"
      } else if (p_val_b < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b, 3))
      }
      
      cat("拟合成功: a =", a_val, "b =", b_val, p_text, "\n")
      cat("残差统计: 均值 =", round(mean(df$residual, na.rm = TRUE), 3), 
          "标准差 =", round(sd(df$residual, na.rm = TRUE), 3), "\n")
      
    } else {
      # 如果拟合失败，设置默认值
      df$fitted_value <- NA
      df$residual <- NA
      cat("未进行拟合，残差设为NA\n")
    }
    
    # 添加距离类型信息
    df$distance_type <- x_cols[i]
    df$assay <- assay_sele
    
    # 存储结果
    result_list[[i]] <- df
    
    cat("----------------------------------------\n")
  }
  
  # 5. 合并所有结果
  if (length(result_list) > 0) {
    all_results <- rbindlist(result_list, fill = TRUE)
    
    # 查看结果摘要
    cat("\n=== 结果摘要 ===\n")
    cat("总数据点:", nrow(all_results), "\n")
    cat("有效拟合点:", sum(!is.na(all_results$residual)), "\n")
    
    # 按距离类型查看残差统计
    residual_summary <- all_results[!is.na(residual), .(
      n_points = .N,
      mean_residual = round(mean(residual), 3),
      sd_residual = round(sd(residual), 3),
      min_residual = round(min(residual), 3),
      max_residual = round(max(residual), 3)
    ), by = distance_type]
    
    print(residual_summary)
    
    return(all_results)
    
  } else {
    cat("没有成功生成结果\n")
    return(NULL)
  }
}

# ===============================
# 主函数：两种 assay 的残差相关性分析与绘图（不排除结合界面）
# ===============================
plot_mutation_residual_correlation_between_assays <- function(
    input1, assay1,
    input2, assay2,
    anno_file,
    output_pdf,
    interface_sites1 = NULL,
    interface_sites2 = NULL
) {
  
  # ========== Step 1: 计算两个 assay 的残差 ==========
  cat("=== 计算", assay1, "的残差 ===\n")
  result1 <- calculate_distance_residuals(input1, assay1, anno_file)
  
  cat("=== 计算", assay2, "的残差 ===\n")
  result2 <- calculate_distance_residuals(input2, assay2, anno_file)
  
  # 检查结果是否为空
  if (is.null(result1) || is.null(result2)) {
    cat("❌ 一个或多个assay的残差计算失败\n")
    return(NULL)
  }
  
  # ========== Step 2: 处理第一个 assay 数据 ==========
  ddG1 <- fread(input1)
  ddG1 <- ddG1[, Pos_real := Pos + 1]
  ddG1 <- ddG1[, c(1, 20, 23)]
  colnames(ddG1)[2] <- paste0(colnames(ddG1)[2], "_", assay1)
  
  # 重命名result1的列
  setnames(result1, "original_y", paste0("mean_kcal/mol_", assay1))
  
  # 合并第一个assay数据
  merge1 <- merge(ddG1, result1, by = c("Pos_real", paste0("mean_kcal/mol_", assay1)), all = FALSE)
  
  # ========== Step 3: 处理第二个 assay 数据 ==========
  ddG2 <- fread(input2)
  ddG2 <- ddG2[, Pos_real := Pos + 1]
  ddG2 <- ddG2[, c(1, 20, 23)]
  colnames(ddG2)[2] <- paste0(colnames(ddG2)[2], "_", assay2)
  
  # 重命名result2的列
  setnames(result2, "original_y", paste0("mean_kcal/mol_", assay2))
  
  # 合并第二个assay数据
  merge2 <- merge(ddG2, result2, by = c("Pos_real", paste0("mean_kcal/mol_", assay2)), all = FALSE)
  
  # ========== Step 4: 重命名残差列并合并 ==========
  setnames(merge1, "residual", paste0("residual_", assay1))
  setnames(merge2, "residual", paste0("residual_", assay2))
  
  # 使用明确的列名进行合并
  data_plot <- merge(merge1, merge2, 
                     by = c("Pos_real", "id", "Pos"), 
                     all = FALSE)
  
  cat("合并后数据维度:", dim(data_plot), "\n")
  
  # ========== Step 5: 检查残差列是否存在 ==========
  residual_col1 <- paste0("residual_", assay1)
  residual_col2 <- paste0("residual_", assay2)
  
  if (!all(c(residual_col1, residual_col2) %in% names(data_plot))) {
    cat("❌ 错误: 残差列不存在！可用列名为:\n")
    print(names(data_plot))
    stop("请检查列名设置")
  }
  
  # ========== Step 6: 不排除结合界面位点，使用所有数据 ==========
  data_plot_filtered <- data_plot  # 直接使用所有数据，不进行过滤
  
  cat("数据使用统计:\n")
  cat("总数据点:", nrow(data_plot_filtered), "\n")
  if (!is.null(interface_sites1) && !is.null(interface_sites2)) {
    all_interface_sites <- unique(c(interface_sites1, interface_sites2))
    interface_count <- sum(data_plot_filtered$Pos_real %in% all_interface_sites)
    cat("包含的结合界面位点数量:", interface_count, "\n")
    cat("结合界面位点:", paste(all_interface_sites, collapse = ", "), "\n")
  } else {
    cat("未提供结合界面位点信息\n")
  }
  
  # ========== Step 7: 准备相关性分析数据 ==========
  df <- data_plot_filtered %>%
    select(residual1 = all_of(residual_col1),
           residual2 = all_of(residual_col2),
           Pos_real) %>%  # 保留Pos_real用于标记
    filter(!is.na(residual1) & !is.na(residual2))
  
  cat("去除NA后数据点:", nrow(df), "\n")
  
  if (nrow(df) == 0) {
    cat("❌ 没有有效数据点计算相关性\n")
    return(NULL)
  }
  
  # ========== Step 8: 计算相关系数 ==========
  cor_value <- cor(df$residual1, df$residual2, method = "pearson")
  cat("Pearson correlation =", round(cor_value, 3), "\n")
  
  # ========== Step 9: 为位点添加颜色信息 ==========
  # 定义颜色 - X轴assay用红色，Y轴assay用蓝色
  x_axis_color <- "#F4270C"  # 红色
  y_axis_color <- "#1B38A6"  # 蓝色
  both_color <- "#C68EFD"    # 紫色（重叠位点）
  other_color <- "#75C2F6"   # 浅蓝色（其他位点）
  
  # 为每个位点添加颜色信息 - 使用固定的类别名称
  df$point_color <- "Other"
  
  if (!is.null(interface_sites1) && !is.null(interface_sites2)) {
    # 属于X轴assay界面位点的点
    df$point_color[df$Pos_real %in% interface_sites1] <- "X_Interface"
    # 属于Y轴assay界面位点的点
    df$point_color[df$Pos_real %in% interface_sites2] <- "Y_Interface"
    # 同时属于两个界面位点的点
    df$point_color[df$Pos_real %in% intersect(interface_sites1, interface_sites2)] <- "Both_Interfaces"
  }
  
  # 设置颜色映射
  color_mapping <- c(
    "Other" = other_color,
    "Both_Interfaces" = both_color,
    "X_Interface" = x_axis_color,
    "Y_Interface" = y_axis_color
  )
  
  # 创建图例标签
  legend_labels <- c(
    "Other" = "Other Sites",
    "Both_Interfaces" = "Both Interfaces",
    "X_Interface" = paste0(assay1, " Interface"),
    "Y_Interface" = paste0(assay2, " Interface")
  )
  
  # ========== Step 10: 绘图 ==========
  # 先绘制Other点（在最底层）
  p <- ggplot(df, aes(x = residual1, y = residual2, color = point_color)) +
    # 先绘制Other点
    geom_point(data = subset(df, point_color == "Other"), 
               alpha = 0.7, size = 2) +
    # 再绘制其他颜色的点
    geom_point(data = subset(df, point_color != "Other"), 
               alpha = 0.7, size = 2) +
    scale_color_manual(
      name = "Binding Interface",
      values = color_mapping,
      labels = legend_labels
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", fill = "#FFB0A5", alpha = 0.3) +
    stat_cor(method = "pearson", 
             label.x = min(df$residual1, na.rm = TRUE), 
             label.y = max(df$residual2, na.rm = TRUE) - 0.5, 
             size = 2.5,
             color = "black",
             label.sep = "\n") +
    labs(
      x = paste0("Residual ", assay1, " (kcal/mol)"),
      y = paste0("Residual ", assay2, " (kcal/mol)"),
      title = paste0("Correlation of Residuals between ", assay1, " and ", assay2,
                     "\n(All Sites Included)")
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6)
    ) +
    coord_cartesian(xlim = c(-2, 2.8), ylim = c(-1.7, 2.8)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  # ========== Step 11: 保存图表 ==========
  ggsave(output_pdf, p, device = cairo_pdf, width = 5, height = 5)
  
  # ========== Step 12: 返回结果 ==========
  return(list(
    plot = p,
    correlation = cor_value,
    data = df,
    merged_data = data_plot_filtered,
    color_info = list(
      x_axis_color = x_axis_color,
      y_axis_color = y_axis_color,
      both_color = both_color,
      other_color = other_color
    )
  ))
}
# ===============================
# 使用示例
# ===============================
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay1 = "RALGDS",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_RALGDS_vs_RAF1_with_2BIs.pdf",
  interface_sites1 = RALGDS_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}



#### PI3KCG vs RAF1
# 定义结合界面位点
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay1 = "PI3KCG",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_PI3KCG_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = PI3KCG_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}



#### SOS1 vs RAF1
# 定义结合界面位点
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay1 = "SOS1",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_SOS1_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = SOS1_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}



#### K55 vs RAF1
# 定义结合界面位点
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay1 = "K55",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K55_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = K55_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}


#### K27 vs RAF1
# 定义结合界面位点
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K27_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K13 vs RAF1
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}



#### K19 vs RAF1
# 定义结合界面位点
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_RAF1 with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = RAF1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





################===========================================================
#### PI3KCG vs RALGDS
# 定义结合界面位点
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay1 = "PI3KCG",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_PI3KCG_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = PI3KCG_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### SOS1 vs RALGDS
# 定义结合界面位点
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay1 = "SOS1",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_SOS1_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = SOS1_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### K55 vs RALGDS
# 定义结合界面位点
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay1 = "K55",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K55_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = K55_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





#### K27 vs RALGDS
# 定义结合界面位点
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K27_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





#### K13 vs RALGDS
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





#### K19 vs RALGDS
# 定义结合界面位点
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
RALGDS_Binding_interface_site <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_RALGDS with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = RALGDS_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### SOS1 vs PI3KCG
# 定义结合界面位点
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay1 = "SOS1",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay2 = "PI3KCG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_SOS1_vs_PI3KCG with 2BIs.pdf",
  interface_sites1 = SOS1_Binding_interface_site,
  interface_sites2 = PI3KCG_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K55 vs PI3KCG
# 定义结合界面位点
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay1 = "K55",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay2 = "PI3KCG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K55_vs_PI3KCG with 2BIs.pdf",
  interface_sites1 = K55_Binding_interface_site,
  interface_sites2 = PI3KCG_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





#### K27 vs PI3KCG
# 定义结合界面位点
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay2 = "PI3KCG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K27_vs_PI3KCG with 2BIs.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = PI3KCG_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### K13 vs PI3KCG
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay2 = "PI3KCG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_PI3KCG with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = PI3KCG_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K19 vs PI3KCG
# 定义结合界面位点
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
PI3KCG_Binding_interface_site <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay2 = "PI3KCG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_PI3KCG with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = PI3KCG_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### K55 vs SOS1
# 定义结合界面位点
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay1 = "K55",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay2 = "SOS1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K55_vs_SOS1 with 2BIs.pdf",
  interface_sites1 = K55_Binding_interface_site,
  interface_sites2 = SOS1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K27 vs SOS1
# 定义结合界面位点
K27_Binding_interface_site <- c( 21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay2 = "SOS1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K27_vs_SOS1 with 2BIs.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = SOS1_Binding_interface_site
)


if (!is.null(result)) {
  print(result$p)}





#### K13 vs SOS1
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay2 = "SOS1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_SOS1 with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = SOS1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### K19 vs SOS1
# 定义结合界面位点
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
SOS1_Binding_interface_site <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay2 = "SOS1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_SOS1 with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = SOS1_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}





#### K27 vs K55
# 定义结合界面位点
K27_Binding_interface_site <- c( 21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay2 = "K55",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K27_vs_K55 with 2BIs.pdf",
  interface_sites1 = K27_Binding_interface_site,
  interface_sites2 = K55_Binding_interface_site
)


if (!is.null(result)) {
  print(result$p)}







#### K13 vs K55
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay2 = "K55",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_K55 with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = K55_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K19 vs K55
# 定义结合界面位点
K19_Binding_interface_site <- c( 68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
K55_Binding_interface_site <- c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay2 = "K55",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_K55 with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = K55_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}






#### K13 vs K27
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay2 = "K27",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_K27 with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = K27_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}



#### K19 vs K27
# 定义结合界面位点
K19_Binding_interface_site <- c( 68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay2 = "K27",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K19_vs_K27 with 2BIs.pdf",
  interface_sites1 = K19_Binding_interface_site,
  interface_sites2 = K27_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}




#### K13 vs K19
# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K19_Binding_interface_site <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)

result <- plot_mutation_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay2 = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251103/mutation_residues_Correlation_K13_vs_K19 with 2BIs.pdf",
  interface_sites1 = K13_Binding_interface_site,
  interface_sites2 = K19_Binding_interface_site
)

if (!is.null(result)) {
  print(result$p)}

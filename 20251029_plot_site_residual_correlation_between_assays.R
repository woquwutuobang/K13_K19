library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ===============================
# 主函数：位点水平的残差相关性分析与绘图
# ===============================
plot_site_residual_correlation_between_assays <- function(
    input1, assay1,
    input2, assay2,
    anno_file,
    wt_aa,
    output_pdf,
    interface_sites = list()
) {
  
  # ========== 内部函数：ddG数据处理 ==========
  ddG_data_process <- function(input, wt_aa) {
    ddG <- fread(input)
    num <- nchar(wt_aa) + 1
    aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
    ddG[, Pos_real := Pos_ref + 1]
    ddG[id != "WT", `:=`(
      wt_codon = substr(id, 1, 1),
      mt_codon = substr(id, nchar(id), nchar(id))
    )]
    ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
    heatmap_tool <- data.table(
      wt_codon = rep(unlist(strsplit(wt_aa, "")), each = 20),
      Pos_real = rep(2:num, each = 20),
      mt_codon = unlist(aa_list)
    )
    ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
    codon <- ddG[Pos_real > 1, unique(wt_codon), by = Pos_real]; setnames(codon, "V1", "codon")
    mean <- ddG[Pos_real > 1, sum(`mean_kcal/mol`/(`std_kcal/mol`^2), na.rm=TRUE)/sum(1/(`std_kcal/mol`^2), na.rm=TRUE), by=Pos_real]; setnames(mean, "V1", "mean")
    abs_mean <- ddG[Pos_real > 1, sum(abs(`mean_kcal/mol`)/(`std_kcal/mol`^2), na.rm=TRUE)/sum(1/(`std_kcal/mol`^2), na.rm=TRUE), by=Pos_real]; setnames(abs_mean, "V1", "abs_mean")
    sigma <- ddG[Pos_real > 1, sqrt(1/sum(1/(`std_kcal/mol`^2), na.rm=TRUE)), by=Pos_real]; setnames(sigma, "V1", "sigma")
    max <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), max(`mean_kcal/mol`), by=Pos_real]; setnames(max, "V1", "max")
    min <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), min(`mean_kcal/mol`), by=Pos_real]; setnames(min, "V1", "min")
    count <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), .N, by=Pos_real]; setnames(count, "N", "count")
    median_ddG <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), median(`mean_kcal/mol`, na.rm=TRUE), by=Pos_real]; setnames(median_ddG, "V1", "median")
    abs_median_ddG <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), median(abs(`mean_kcal/mol`), na.rm=TRUE), by=Pos_real]; setnames(abs_median_ddG, "V1", "abs_median")
    output <- Reduce(function(x, y) merge(x, y, by = "Pos_real", all = TRUE),
                     list(codon, mean, abs_mean, sigma, max, min, count, median_ddG, abs_median_ddG))
    return(output)
  }
  
  # ========== 内部函数：距离效应指数拟合和残差计算 ==========
  fit_distance_exponential <- function(input_file, anno_file, wt_aa, assay_name) {
    
    # 1. 处理ddG数据
    ddG_data <- ddG_data_process(input_file, wt_aa)
    
    # 2. 读取注释文件
    anno <- fread(anno_file)
    anno[, Pos_real := Pos]
    
    # 3. 合并数据
    merged_data <- merge(ddG_data, anno, by = "Pos_real", all = FALSE)
    
    # 4. 提取距离和能量数据
    distance_col <- paste0("scHAmin_ligand_", assay_name)
    x <- merged_data[[distance_col]]
    y <- abs(merged_data$mean)  # 取绝对值
    
    # 去掉 NA
    valid_idx <- complete.cases(x, y)
    x <- x[valid_idx]
    y <- y[valid_idx]
    
    df <- data.frame(
      Pos_real = merged_data$Pos_real[valid_idx],
      x = x, 
      y = y,
      original_y = merged_data$mean[valid_idx]  # 保留原始能量值（有正负号）
    )
    
    cat("有效数据点:", nrow(df), "\n")
    
    # 5. 初始参数优化
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]
      b <- params[2]
      sum((y - a * exp(b * x))^2, na.rm = TRUE)
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
    
    # 6. 拟合指数模型
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), data = df, start = list(a = opt_params[1], b = opt_params[2])),
      error = function(e) {
        cat("模型拟合失败:", e$message, "\n")
        NULL
      }
    )
    
    # 7. 计算拟合值和残差
    if (!is.null(fit_model)) {
      # 拟合值和残差
      df$fitted_value <- predict(fit_model, newdata = df)
      df$residual <- df$original_y - df$fitted_value  # 用原始值（有正负号）计算残差
      
      # 拟合参数
      fit_summary <- summary(fit_model)
      coefs <- fit_summary$coefficients
      a_val <- coefs["a", "Estimate"]
      b_val <- coefs["b", "Estimate"]
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
      
      cat("拟合成功: a =", round(a_val, 3), "b =", round(b_val, 3), p_text, "\n")
      cat("残差统计: 均值 =", round(mean(df$residual, na.rm = TRUE), 3), 
          "标准差 =", round(sd(df$residual, na.rm = TRUE), 3), "\n")
      
    } else {
      # 如果拟合失败，设置默认值
      df$fitted_value <- NA
      df$residual <- NA
      cat("未进行拟合，残差设为NA\n")
    }
    
    # 8. 添加元数据
    df$distance_type <- distance_col
    df$assay <- assay_name
    
    # 9. 返回结果
    return(list(
      data = df,
      model = fit_model,
      summary = if (!is.null(fit_model)) summary(fit_model) else NULL
    ))
  }
  
  # ========== Step 1: 分别计算两个 assay 的残差 ==========
  cat("正在处理", assay1, "...\n")
  result1 <- fit_distance_exponential(input1, anno_file, wt_aa, assay1)
  
  cat("正在处理", assay2, "...\n")
  result2 <- fit_distance_exponential(input2, anno_file, wt_aa, assay2)
  
  # 提取残差数据
  data1 <- result1$data
  data2 <- result2$data
  
  # ========== Step 2: 读取原始数据并合并 ==========
  ddG_site1 <- ddG_data_process(input1, wt_aa)
  colnames(ddG_site1)[3:10] <- paste0(colnames(ddG_site1)[3:10], "_", assay1)
  setnames(data1, "original_y", paste0("mean_", assay1))
  merge1 <- merge(ddG_site1, data1, by = c("Pos_real", paste0("mean_", assay1)), all = FALSE)
  
  ddG_site2 <- ddG_data_process(input2, wt_aa)
  colnames(ddG_site2)[3:10] <- paste0(colnames(ddG_site2)[3:10], "_", assay2)
  setnames(data2, "original_y", paste0("mean_", assay2))
  merge2 <- merge(ddG_site2, data2, by = c("Pos_real", paste0("mean_", assay2)), all = FALSE)
  
  # ========== Step 3: 合并两者 ==========
  merge1 <- merge1[, c("Pos_real", "residual")]
  merge2 <- merge2[, c("Pos_real", "residual")]
  
  setnames(merge1, "residual", paste0("residual_", assay1))
  setnames(merge2, "residual", paste0("residual_", assay2))
  
  data_plot <- merge(merge1, merge2, by = "Pos_real", all = FALSE)
  
  # ========== Step 4: 过滤结合界面 ==========
  if (length(interface_sites) > 0) {
    all_interface_sites <- unique(unlist(interface_sites))
    data_plot_filtered <- data_plot[!Pos_real %in% all_interface_sites]
  } else {
    data_plot_filtered <- data_plot
  }
  
  cat("过滤前数据点:", nrow(data_plot), "\n")
  cat("过滤后数据点:", nrow(data_plot_filtered), "\n")
  if (length(interface_sites) > 0) {
    cat("排除的结合界面位点数量:", length(all_interface_sites), "\n")
  }
  cat("剩余位点数量:", nrow(data_plot_filtered), "\n")
  
  # ========== Step 5: 计算相关系数 ==========
  residual_col1 <- paste0("residual_", assay1)
  residual_col2 <- paste0("residual_", assay2)
  
  df <- data_plot_filtered %>%
    select(residual1 = all_of(residual_col1),
           residual2 = all_of(residual_col2)) %>%
    filter(!is.na(residual1) & !is.na(residual2))
  
  cor_value <- cor(df$residual1, df$residual2, method = "pearson")
  cat("Pearson correlation =", round(cor_value, 3), "\n")
  
  # ========== Step 6: 绘图 ==========
  p <- ggplot(df, aes(x = residual1, y = residual2)) +
    geom_point(alpha = 0.6, color = "#75C2F6", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", fill = "#FFB0A5", alpha = 0.3) +
    stat_cor(method = "pearson", 
             label.x = min(df$residual1, na.rm = TRUE), 
             label.y = max(df$residual2, na.rm = TRUE) - 0.5, 
             size = 2.5,
             label.sep = "\n") +
    labs(
      x = paste0("Residual ", assay1, " (kcal/mol)"),
      y = paste0("Residual ", assay2, " (kcal/mol)"),
      title = paste0("Correlation of Residuals between ", assay1, " and ", assay2,
                     "\n(Per Site Analysis, Excluding Binding Interface Sites)")
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.border = element_rect(color = "black", linewidth = 0.8)
    ) +
    coord_cartesian(xlim = c(-1.2, 1.6), ylim = c(-1.2, 1.2))
  
  ggsave(output_pdf, p, device = cairo_pdf, width = 4, height = 4)
  
  # 返回结果
  return(list(
    plot = p,
    correlation = cor_value,
    data = df,
    filtered_data = data_plot_filtered
  ))
}

# ===============================
# 使用示例
# ===============================


######K13_SOS1
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

K13_sites <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
SOS1_sites <- c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153)

# 使用函数
result <- plot_site_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay1 = "SOS1",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay2 = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa,
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5b_Correlation_of_Residuals_per_sites SOS1 VS K13.pdf",
  interface_sites = list(K13 = K13_sites, SOS1 = SOS1_sites)
)

# 显示结果
print(result$plot)
cat("最终相关系数:", round(result$correlation, 3), "\n")





###### K19 VS RALGDS
# 定义WT序列和结合界面位点
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

K19_sites <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
RALGDS_sites <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)

# 使用函数
result <- plot_site_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K19",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay2 = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa,
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5b_Correlation_of_Residuals_per_sites K19 VS RALGDS.pdf",
  interface_sites = list(K19 = K19_sites, RALGDS = RALGDS_sites)
)

# 显示结果
print(result$plot)
cat("最终相关系数:", round(result$correlation, 3), "\n")


######K13 VS RAF1
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

K13_sites <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_sites <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

# 使用函数
result <- plot_site_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa,
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5b_Correlation_of_Residuals_per_sites K13 VS RAF1.pdf",
  interface_sites = list(K13 = K13_sites, RAF1 = RAF1_sites)
)

# 显示结果
print(result$plot)
cat("最终相关系数:", round(result$correlation, 3), "\n")




######K27 VS RAF1
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

K27_sites <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
RAF1_sites <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

# 使用函数
result <- plot_site_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay1 = "K27",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa,
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5b_Correlation_of_Residuals_per_sites K27 VS RAF1.pdf",
  interface_sites = list(K27 = K27_sites, RAF1 = RAF1_sites)
)

# 显示结果
print(result$plot)
cat("最终相关系数:", round(result$correlation, 3), "\n")




######K13 VS K19
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

K13_sites <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
K19_sites <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)

# 使用函数
result <- plot_site_residual_correlation_between_assays(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay1 = "K13",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay2 = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa,
  output_pdf = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/f5b_Correlation_of_Residuals_per_sites K13 VS K19.pdf",
  interface_sites = list(K13 = K13_sites, K19 = K19_sites)
)

# 显示结果
print(result$plot)
cat("最终相关系数:", round(result$correlation, 3), "\n")

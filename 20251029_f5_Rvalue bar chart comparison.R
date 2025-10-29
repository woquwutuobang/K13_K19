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
      a <- params[1]; b <- params[2]
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
# 计算所有binder的残差
# ===============================
# 定义所有binder
binders <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K13", "K19")
binder_groups <- list(
  BI1 = c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG"),
  BI2 = c("K13", "K19")
)

# 存储每个binder的处理结果
binder_data <- list()

# 逐个处理每个binder
for (binder in binders) {
  cat("正在处理binder:", binder, "\n")
  
  # 读取原始数据
  ddG_file <- paste0("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_", binder, ".txt")
  ddG_data <- fread(ddG_file)
  
  # 处理Pos_real
  ddG_data <- ddG_data[, Pos_real := Pos + 1]
  ddG_data <- ddG_data[, c(1, 20, 23)]  # 选择需要的列
  colnames(ddG_data)[2] <- paste0(colnames(ddG_data)[2], "_", binder)
  
  # 计算残差
  result <- calculate_distance_residuals(
    input = ddG_file,
    assay_sele = binder,
    anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  )
  
  if (!is.null(result)) {
    # 重命名original_y列
    setnames(result, "original_y", paste0("mean_kcal/mol_", binder))
    
    # 合并原始数据和残差数据
    merged_data <- merge(ddG_data, result, 
                         by = c("Pos_real", paste0("mean_kcal/mol_", binder)), 
                         all = FALSE)
    
    # 重命名残差列
    setnames(merged_data, "residual", paste0("residual_", binder))
    
    # 选择需要的列
    final_cols <- c("Pos_real", "id", "Pos", paste0("residual_", binder))
    binder_data[[binder]] <- merged_data[, .SD, .SDcols = final_cols]
    
    cat("✅", binder, "处理完成，数据维度:", dim(binder_data[[binder]]), "\n")
  } else {
    cat("❌", binder, "处理失败\n")
  }
  cat("========================================\n")
}

# 逐个合并所有binder的数据
cat("=== 开始逐个合并数据 ===\n")
all_residuals <- binder_data[[1]]

for (i in 2:length(binders)) {
  binder <- binders[i]
  if (binder %in% names(binder_data)) {
    cat("正在合并:", binder, "\n")
    cat("  合并前维度:", dim(all_residuals), "\n")
    
    all_residuals <- merge(all_residuals, binder_data[[binder]], 
                           by = c("Pos_real", "id", "Pos"), 
                           all = TRUE)
    
    cat("  合并后维度:", dim(all_residuals), "\n\n")
  }
}

cat("最终合并数据维度:", dim(all_residuals), "\n")
cat("最终列名:", names(all_residuals), "\n")

# 定义结合界面位点
binding_interface_sites <- list(
  RAF1 = c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153),  
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),  
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),   
  RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), 
  PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73), 
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)  
)

# 计算相关性
correlation_results <- data.frame()

successful_binders <- names(binder_data)
cat("成功处理的binder:", successful_binders, "\n")

for (i in 1:(length(successful_binders)-1)) {
  for (j in (i+1):length(successful_binders)) {
    binder1 <- successful_binders[i]
    binder2 <- successful_binders[j]
    
    # 获取两个binder的结合界面位点
    interface_sites1 <- binding_interface_sites[[binder1]]
    interface_sites2 <- binding_interface_sites[[binder2]]
    
    # 过滤掉两个binder各自的结合界面位点
    df_cor <- all_residuals[!Pos_real %in% c(interface_sites1, interface_sites2)]
    
    # 选择两个binder的残差列
    residual_col1 <- paste0("residual_", binder1)
    residual_col2 <- paste0("residual_", binder2)
    
    df_cor <- df_cor[, .SD, .SDcols = c(residual_col1, residual_col2)]
    df_cor <- df_cor[complete.cases(df_cor)]
    
    if (nrow(df_cor) > 5) {
      cor_value <- cor(df_cor[[1]], df_cor[[2]], method = "pearson")
      
      # 确定分组
      group1 <- ifelse(binder1 %in% binder_groups$BI1, "BI1", "BI2")
      group2 <- ifelse(binder2 %in% binder_groups$BI1, "BI1", "BI2")
      
      if (group1 == "BI1" & group2 == "BI1") {
        comparison_group <- "BI1 vs BI1"
      } else if (group1 == "BI2" & group2 == "BI2") {
        comparison_group <- "BI2 vs BI2"
      } else {
        comparison_group <- "BI1 vs BI2"
      }
      
      correlation_results <- rbind(correlation_results, data.frame(
        Binder1 = binder1,
        Binder2 = binder2,
        Correlation = cor_value,
        Group = comparison_group,
        N_points = nrow(df_cor),
        stringsAsFactors = FALSE
      ))
      
      cat(binder1, "vs", binder2, ": r =", round(cor_value, 3), 
          "(n =", nrow(df_cor), ")", "Group:", comparison_group, "\n")
    }
  }
}


# ===============================
# 绘制分组柱状图
# ===============================

# 按分组排序
correlation_results$Group <- factor(correlation_results$Group, 
                                    levels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2"))

# 创建比较标签
correlation_results$Comparison <- paste0(correlation_results$Binder1, " vs ", correlation_results$Binder2)

# 按分组和相关性值排序
correlation_results <- correlation_results[order(correlation_results$Group, correlation_results$Correlation), ]
correlation_results$Comparison <- factor(correlation_results$Comparison, 
                                         levels = correlation_results$Comparison)

# 绘制柱状图
p <- ggplot(correlation_results, aes(x = Comparison, y = Correlation, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(Correlation, 3)), 
            vjust = -0.5, size = 2.5, color = "black") +
  scale_fill_manual(values = c("BI1 vs BI1" = "#F4270C", 
                               "BI1 vs BI2" = "#F4AD0C", 
                               "BI2 vs BI2" = "#1B38A6")) +
  labs(
    title = "Correlation of Residuals between Different Bindersn\n(Per Mutation Analysis, Excluding Binding Interface Sites)",
    x = "Binder Comparison",
    y = "Pearson Correlation Coefficient (r)",
    fill = "Comparison Group"
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
  ylim(min(correlation_results$Correlation, na.rm = TRUE) - 0.1, 
       max(correlation_results$Correlation, na.rm = TRUE) + 0.1)

print(p)

# 保存图表
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251029/Binder_Residual_Correlations_Grouped_per_mutation 2.pdf",
       p, device = cairo_pdf, width = 10, height = 6)

# 输出汇总统计
cat("\n=== 相关性汇总统计 ===\n")
summary_stats <- correlation_results %>%
  group_by(Group) %>%
  summarise(
    Count = n(),
    Mean_Correlation = round(mean(Correlation), 3),
    SD_Correlation = round(sd(Correlation), 3),
    Min_Correlation = round(min(Correlation), 3),
    Max_Correlation = round(max(Correlation), 3)
  )

print(summary_stats)

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
binders <- c("RAF1","RALGDS","PI3KCG","SOS1", "K55", "K27","K13", "K19")

binder_groups <- list(
  BI1 = c("RAF1","RALGDS","PI3KCG","SOS1", "K55", "K27"),
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


binding_interface_sites <- list(RAF1 = c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
                                RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), 
                                PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73),
                                SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153),
                                K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),  
                                K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
                                K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
                                K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))


# ===============================
# 计算相关性（包含p值计算）
# ===============================
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
      # 计算相关系数和p值
      cor_test <- cor.test(df_cor[[1]], df_cor[[2]], method = "pearson")
      cor_value <- cor_test$estimate
      p_value <- cor_test$p.value
      
      # 添加显著性标记
      significance <- case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      )
      
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
        P_value = p_value,
        Significance = significance,
        Group = comparison_group,
        N_points = nrow(df_cor),
        stringsAsFactors = FALSE
      ))
      
      cat(binder1, "vs", binder2, ": r =", round(cor_value, 3), 
          "p =", round(p_value, 4), significance,
          "(n =", nrow(df_cor), ")", "Group:", comparison_group, "\n")
    }
  }
}




###################==================================================
##### plot2---chatgpt version

# ============================================================
# Correlation Sorting & Visualization Script (Multi-layer logic)
# ===============================

# ===============================
# 1️⃣ 分类定义
# ===============================
binders <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

non_darpin <- c("RAF1", "RALGDS", "PI3KCG", "SOS1")
darpin     <- c("K55", "K27", "K13", "K19")

active     <- c("RAF1", "RALGDS", "PI3KCG",  "K55")
inactive   <- c("SOS1","K27", "K13", "K19")

# 定义BI1和BI2
BI1 <- c("RAF1", "RALGDS", "PI3KCG", "SOS1","K55", "K27")
BI2 <- c("K13", "K19")

# ===============================
# 2️⃣ 数据准备
# ===============================
# 创建 Comparison 列
correlation_results$Comparison <- paste0(correlation_results$Binder1, " vs ", correlation_results$Binder2)

# 分类：DARPin / non-DARPin
correlation_results$Type1 <- ifelse(correlation_results$Binder1 %in% darpin, "DARPin", "non-DARPin")
correlation_results$Type2 <- ifelse(correlation_results$Binder2 %in% darpin, "DARPin", "non-DARPin")

# 分类：Active / Inactive
correlation_results$Act1 <- ifelse(correlation_results$Binder1 %in% active, "Active", "Inactive")
correlation_results$Act2 <- ifelse(correlation_results$Binder2 %in% active, "Active", "Inactive")

# ===============================
# 3️⃣ 层级标签生成
# ===============================
# 生成comparison_type分类 (第一层)
correlation_results$comparison_type <- apply(correlation_results, 1, function(x) {
  binder1 <- x["Binder1"]
  binder2 <- x["Binder2"]
  
  if (binder1 %in% BI1 && binder2 %in% BI1) {
    return("BI1_vs_BI1")
  } else if (binder1 %in% BI2 && binder2 %in% BI2) {
    return("BI2_vs_BI2")
  } else {
    return("BI1_vs_BI2")
  }
})

# 一级：DARPin 分类 (第二层)
correlation_results$DARPinClass <- with(
  correlation_results,
  ifelse(Type1 == "non-DARPin" & Type2 == "non-DARPin", "nonDARPin_vs_nonDARPin",
         ifelse(Type1 == "DARPin" & Type2 == "DARPin", "DARPin_vs_DARPin", "nonDARPin_vs_DARPin"))
)

# 二级：Activity 分类 (第三层)
correlation_results$ActivityClass <- with(
  correlation_results,
  ifelse(Act1 == "Active" & Act2 == "Active", "Active_vs_Active",
         ifelse(Act1 == "Inactive" & Act2 == "Inactive", "Inactive_vs_Inactive", "Active_vs_Inactive"))
)

# 组合排序标签（三层组合）
correlation_results$SortKey <- paste(correlation_results$comparison_type, 
                                     correlation_results$DARPinClass, 
                                     correlation_results$ActivityClass, sep = "_")

# ===============================
# 4️⃣ 三层层级排序
# ===============================
correlation_results <- correlation_results[order(
  # 第一层：comparison_type (BI1_vs_BI1, BI1_vs_BI2, BI2_vs_BI2)
  factor(correlation_results$comparison_type, levels = c(
    "BI1_vs_BI1",
    "BI1_vs_BI2", 
    "BI2_vs_BI2"
  )),
  # 第二层：DARPin分类
  factor(correlation_results$DARPinClass, levels = c(
    "nonDARPin_vs_nonDARPin",
    "nonDARPin_vs_DARPin",
    "DARPin_vs_DARPin"
  )),
  # 第三层：活性分类
  factor(correlation_results$ActivityClass, levels = c(
    "Active_vs_Active",
    "Active_vs_Inactive", 
    "Inactive_vs_Inactive"
  )),
  # 最后按比较名称排序以确保一致性
  factor(correlation_results$Comparison, levels = unique(correlation_results$Comparison))
), ]

# 重新设置因子顺序
correlation_results$Comparison <- factor(
  correlation_results$Comparison,
  levels = correlation_results$Comparison
)

# 设置comparison_type因子顺序
correlation_results$comparison_type <- factor(
  correlation_results$comparison_type, 
  levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2")
)

# ===============================
# 5️⃣ Group 因子顺序保持一致
# ===============================
if ("Group" %in% colnames(correlation_results)) {
  correlation_results$Group <- factor(
    correlation_results$Group, 
    levels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")
  )
}

# ===============================
# 6️⃣ 绘图（添加显著性标记）
# ===============================
p <- ggplot(correlation_results, aes(x = Comparison, y = Correlation, fill = comparison_type)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(Correlation, 3)),
            vjust = -0.5, size = 2.5, color = "black") +
  # 添加显著性标记
  geom_text(aes(label = Significance), 
            vjust = -1.5, size = 3, color = "#F1DD10") +
  scale_fill_manual(values = c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C",
    "BI2_vs_BI2" = "#1B38A6"
  )) +
  labs(
    title = "Correlation of Residuals between Different Binders\n(Per Mutation Analysis, Excluding Binding Interface Sites)",
    x = "Binder Comparison",
    y = "Pearson Correlation Coefficient (r)",
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
  ylim(
    min(correlation_results$Correlation, na.rm = TRUE) - 0.1,
    max(correlation_results$Correlation, na.rm = TRUE) + 0.1
  )

print(p)

# ===============================
# 7️⃣ 保存图表
# ===============================
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251106/Binder_Residual_Correlations_three_layer_sorted_per_mutation_with_significance delete 2BI.pdf",
  plot = p, device = cairo_pdf, width = 12, height = 8
)

# ===============================
# 保存相关性结果
# ===============================
fwrite(correlation_results, "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251103/correlation_results_three_layer_sorted_with_significance.csv")

cat("分析完成！相关性结果和图形已保存。\n")
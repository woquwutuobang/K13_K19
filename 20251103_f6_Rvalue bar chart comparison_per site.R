library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

#####
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

# ===============================
# 函数：距离效应指数拟合和残差计算
# ===============================
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
  
  # 计算残差
  result <- fit_distance_exponential(
    input_file = paste0("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_", binder, ".txt"),
    anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
    wt_aa = wt_aa,
    assay_name = binder
  )
  
  if (!is.null(result)) {
    # 提取残差数据
    residual_data <- result$data[, c("Pos_real", "residual")]
    setnames(residual_data, "residual", paste0("residual_", binder))
    binder_data[[binder]] <- residual_data
    cat("✅", binder, "处理完成，数据维度:", dim(residual_data), "\n")
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
                           by = "Pos_real", 
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

# 计算相关性
correlation_results <- data.frame()
successful_binders <- names(binder_data)
cat("成功处理的binder:", successful_binders, "\n")

# 确保 all_residuals 是 data.table
if (!is.data.table(all_residuals)) {
  cat("将 all_residuals 转换为 data.table\n")
  all_residuals <- as.data.table(all_residuals)
}

cat("all_residuals 的类型:", class(all_residuals), "\n")
cat("all_residuals 的列名:", names(all_residuals), "\n")

for (i in 1:(length(successful_binders)-1)) {
  for (j in (i+1):length(successful_binders)) {
    binder1 <- successful_binders[i]
    binder2 <- successful_binders[j]
    
    cat("处理:", binder1, "vs", binder2, "\n")
    
    # 检查两个binder是否都在binding_interface_sites中
    if (!(binder1 %in% names(binding_interface_sites))) {
      cat("  ❌ 警告:", binder1, "不在binding_interface_sites中\n")
      next
    }
    if (!(binder2 %in% names(binding_interface_sites))) {
      cat("  ❌ 警告:", binder2, "不在binding_interface_sites中\n")
      next
    }
    
    # 获取两个binder的结合界面位点
    interface_sites1 <- binding_interface_sites[[binder1]]
    interface_sites2 <- binding_interface_sites[[binder2]]
    
    # 使用安全的data.table语法过滤
    df_cor <- all_residuals[!Pos_real %in% c(interface_sites1, interface_sites2)]
    
    cat("  过滤后剩余位点:", nrow(df_cor), "个\n")
    
    # 选择两个binder的残差列
    residual_col1 <- paste0("residual_", binder1)
    residual_col2 <- paste0("residual_", binder2)
    
    # 确保列存在
    if (!(residual_col1 %in% names(df_cor)) | !(residual_col2 %in% names(df_cor))) {
      cat("  ❌ 跳过", binder1, "vs", binder2, ": 残差列不存在\n")
      next
    }
    
    # 提取需要的列
    df_cor_subset <- df_cor[, .(res1 = get(residual_col1), res2 = get(residual_col2))]
    df_cor_subset <- df_cor_subset[complete.cases(df_cor_subset)]
    
    cat("  去除NA后剩余:", nrow(df_cor_subset), "个数据点\n")
    
    if (nrow(df_cor_subset) > 5) {
      cor_value <- cor(df_cor_subset$res1, df_cor_subset$res2, method = "pearson")
      
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
        N_points = nrow(df_cor_subset),
        stringsAsFactors = FALSE
      ))
      
      cat("  ✅", binder1, "vs", binder2, ": r =", round(cor_value, 3), 
          "(n =", nrow(df_cor_subset), ")", "Group:", comparison_group, "\n")
    } else {
      cat("  ❌ 跳过", binder1, "vs", binder2, ": 数据点不足 (n =", nrow(df_cor_subset), ")\n")
    }
    cat("  ---\n")
  }
}

# 检查是否有结果
if (nrow(correlation_results) > 0) {
  cat("\n=== 相关性计算完成 ===\n")
  cat("总共计算了", nrow(correlation_results), "对binder的相关性\n")}
  
  
  
  
  # ===============================
  # 绘制分组柱状图（使用已定义的 binders & binder_groups）
  # ===============================
  
  # 创建比较标签
  correlation_results$Comparison <- paste0(correlation_results$Binder1, " vs ", correlation_results$Binder2)
  
  # 指定柱子的顺序（按 binders 顺序展开）
  comparison_levels <- c()
  for (i in binders) {
    for (j in binders) {
      comparison_levels <- c(comparison_levels, paste0(i, " vs ", j))
    }
  }
  
  # 只保留在数据中的组合并按顺序排列
  comparison_levels <- comparison_levels[comparison_levels %in% correlation_results$Comparison]
  correlation_results$Comparison <- factor(correlation_results$Comparison, levels = comparison_levels)
  
  # Group 因子顺序保持不变
  correlation_results$Group <- factor(
    correlation_results$Group, 
    levels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")
  )
  
  # 绘制柱状图
  p <- ggplot(correlation_results, aes(x = Comparison, y = Correlation, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = round(Correlation, 3)), 
              vjust = -0.5, size = 2.5, color = "black") +
    scale_fill_manual(values = c(
      "BI1 vs BI1" = "#F4270C", 
      "BI1 vs BI2" = "#F4AD0C", 
      "BI2 vs BI2" = "#1B38A6"
    )) +
    labs(
      title = "Correlation of Residuals between Different Binders\n(Per Mutation Analysis, Excluding Binding Interface Sites)",
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
  ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251103/Binder_Residual_Correlations_Grouped_per_site_fixed_order.pdf",
         p, device = cairo_pdf, width = 10, height = 6)
  
  
  
  
  
  
  ###################==================================================
  ##### plot2
  
  # ============================================================
  # Correlation Sorting & Visualization Script (Multi-layer logic)
  # ============================================================
  
  
  # ===============================
  # 1️⃣ 分类定义
  # ===============================
  binders <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  
  non_darpin <- c("RAF1", "RALGDS", "PI3KCG", "SOS1")
  darpin     <- c("K55", "K27", "K13", "K19")
  
  active     <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55")
  inactive   <- c("K27", "K13", "K19")
  
  # ===============================
  # 2️⃣ 数据准备
  # ===============================
  # 假设 correlation_results 已经存在，包含至少以下列：
  # Binder1, Binder2, Correlation, Group
  
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
  # 一级：DARPin 分类
  correlation_results$DARPinClass <- with(
    correlation_results,
    ifelse(Type1 == "non-DARPin" & Type2 == "non-DARPin", "nonDARPin_vs_nonDARPin",
           ifelse(Type1 == "DARPin" & Type2 == "DARPin", "DARPin_vs_DARPin", "nonDARPin_vs_DARPin"))
  )
  
  # 二级：Activity 分类
  correlation_results$ActivityClass <- with(
    correlation_results,
    ifelse(Act1 == "Active" & Act2 == "Active", "Active_vs_Active",
           ifelse(Act1 == "Inactive" & Act2 == "Inactive", "Inactive_vs_Inactive", "Active_vs_Inactive"))
  )
  
  # 组合排序标签（DARPin 层次 → Activity 层次）
  correlation_results$SortKey <- paste(correlation_results$DARPinClass, correlation_results$ActivityClass, sep = "_")
  
  # 明确设定因子顺序
  correlation_results$SortKey <- factor(
    correlation_results$SortKey,
    levels = c(
      "nonDARPin_vs_nonDARPin_Active_vs_Active",
      "nonDARPin_vs_nonDARPin_Active_vs_Inactive",
      "nonDARPin_vs_nonDARPin_Inactive_vs_Inactive",
      "nonDARPin_vs_DARPin_Active_vs_Active",
      "nonDARPin_vs_DARPin_Active_vs_Inactive",
      "nonDARPin_vs_DARPin_Inactive_vs_Inactive",
      "DARPin_vs_DARPin_Active_vs_Active",
      "DARPin_vs_DARPin_Active_vs_Inactive",
      "DARPin_vs_DARPin_Inactive_vs_Inactive"
    )
  )
  
  # ===============================
  # 4️⃣ 仅按前两层逻辑排序（DARPin + Activity）
  # ===============================
  
  correlation_results <- correlation_results[order(
    factor(correlation_results$DARPinClass, levels = c(
      "nonDARPin_vs_nonDARPin",
      "nonDARPin_vs_DARPin",
      "DARPin_vs_DARPin"
    )),
    factor(correlation_results$ActivityClass, levels = c(
      "Active_vs_Active",
      "Active_vs_Inactive",
      "Inactive_vs_Inactive"
    )),
    factor(correlation_results$Comparison, levels = unique(correlation_results$Comparison))
  ), ]
  
  # 重新设置因子顺序
  correlation_results$Comparison <- factor(
    correlation_results$Comparison,
    levels = correlation_results$Comparison
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
  # 6️⃣ 绘图
  # ===============================
  p <- ggplot(correlation_results, aes(x = Comparison, y = Correlation, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = round(Correlation, 3)),
              vjust = -0.5, size = 2.5, color = "black") +
    scale_fill_manual(values = c(
      "BI1 vs BI1" = "#F4270C",
      "BI1 vs BI2" = "#F4AD0C",
      "BI2 vs BI2" = "#1B38A6"
    )) +
    labs(
      title = "Correlation of Residuals between Different Binders\n(Per Mutation Analysis, Excluding Binding Interface Sites)",
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
    ylim(
      min(correlation_results$Correlation, na.rm = TRUE) - 0.1,
      max(correlation_results$Correlation, na.rm = TRUE) + 0.1
    )
  
  print(p)
  
  # ===============================
  # 7️⃣ 保存图表
  # ===============================
  ggsave(
    filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251103/Binder_Residual_Correlations_Grouped_per_site_multilayer_sorted 2.pdf",
    plot = p, device = cairo_pdf, width = 10, height = 6
  )
  
  
  
  
  #####===========================================
  ### plot 3
  
  
  
  # ===============================
  # 绘制分组柱状图（使用已定义的 binders & binder_groups）
  # ===============================
  
  # 定义分类
  non_darpin <- c("RAF1", "RALGDS", "PI3KCG", "SOS1")
  darpin <- c("K55", "K27", "K13", "K19")
  active <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55")
  inactive <- c("K27", "K13", "K19")
  
  # 整体排序框架
  base_order <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  
  # 创建比较标签
  correlation_results$Comparison <- paste0(correlation_results$Binder1, " vs ", correlation_results$Binder2)
  
  # 创建排序函数
  get_sort_key <- function(b1, b2) {
    # 第一层：DARPin分类
    if (b1 %in% non_darpin & b2 %in% non_darpin) {
      darpin_cat <- 1  # nonDARPin vs nonDARPin
    } else if ((b1 %in% non_darpin & b2 %in% darpin) | (b1 %in% darpin & b2 %in% non_darpin)) {
      darpin_cat <- 2  # nonDARPin vs DARPin
    } else {
      darpin_cat <- 3  # DARPin vs DARPin
    }
    
    # 第二层：活性分类
    if (b1 %in% active & b2 %in% active) {
      activity_cat <- 1  # active vs active
    } else if ((b1 %in% active & b2 %in% inactive) | (b1 %in% inactive & b2 %in% active)) {
      activity_cat <- 2  # active vs inactive
    } else {
      activity_cat <- 3  # inactive vs inactive
    }
    
    # 第三层：在base_order中的位置
    b1_index <- which(base_order == b1)
    b2_index <- which(base_order == b2)
    
    # 返回排序键
    return(paste0(
      sprintf("%01d", darpin_cat),
      sprintf("%01d", activity_cat),
      sprintf("%02d", b1_index),
      sprintf("%02d", b2_index)
    ))
  }
  
  # 为每个比较添加排序键
  correlation_results$SortKey <- sapply(1:nrow(correlation_results), function(i) {
    get_sort_key(correlation_results$Binder1[i], correlation_results$Binder2[i])
  })
  
  # 按排序键排序
  correlation_results <- correlation_results[order(correlation_results$SortKey), ]
  
  # 设置Comparison因子水平
  correlation_results$Comparison <- factor(correlation_results$Comparison, levels = correlation_results$Comparison)
  
  # Group 因子顺序保持不变
  correlation_results$Group <- factor(
    correlation_results$Group, 
    levels = c("BI1 vs BI1", "BI1 vs BI2", "BI2 vs BI2")
  )
  
  # 绘制柱状图
  p <- ggplot(correlation_results, aes(x = Comparison, y = Correlation, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = round(Correlation, 3)), 
              vjust = -0.5, size = 2.5, color = "black") +
    scale_fill_manual(values = c(
      "BI1 vs BI1" = "#F4270C", 
      "BI1 vs BI2" = "#F4AD0C", 
      "BI2 vs BI2" = "#1B38A6"
    )) +
    labs(
      title = "Correlation of Residuals between Different Binders\n(Per Mutation Analysis, Excluding Binding Interface Sites)",
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
  ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251103/Binder_Residual_Correlations_Grouped_per_mutation_fixed_order.pdf",
         p, device = cairo_pdf, width = 10, height = 6)
  
  
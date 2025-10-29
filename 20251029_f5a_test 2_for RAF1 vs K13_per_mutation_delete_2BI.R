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
# 使用示例
# ===============================

# RAF1 assay
result_raf1 <- calculate_distance_residuals(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

# K13 assay
result_k13 <- calculate_distance_residuals(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

assay_sele1 <- "RAF1"
assay_sele2 <- "K13"

# 首先检查原始数据的列名
cat("=== 检查原始数据列名 ===\n")
names(result_raf1)
names(result_k13)

#### 处理RAF1数据 ####
ddG_raf1 <- fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt")
cat("ddG_raf1原始列名:\n")
print(names(ddG_raf1))

ddG_raf1 <- ddG_raf1[, Pos_real := Pos + 1]
ddG_raf1 <- ddG_raf1[, c(1, 20, 23)]
colnames(ddG_raf1)[2] <- paste0(colnames(ddG_raf1)[2], "_", assay_sele1)

cat("ddG_raf1处理后列名:\n")
print(names(ddG_raf1))

# 重命名result_raf1的列
setnames(result_raf1, "original_y", paste0("mean_kcal/mol_", assay_sele1))

cat("result_raf1列名:\n")
print(names(result_raf1))

# 合并RAF1数据
merge_raf1 <- merge(ddG_raf1, result_raf1, by = c("Pos_real", paste0("mean_kcal/mol_", assay_sele1)), all = FALSE)
cat("merge_raf1列名:\n")
print(names(merge_raf1))

#### 处理K13数据 ####
ddG_k13 <- fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt")
cat("ddG_k13原始列名:\n")
print(names(ddG_k13))

ddG_k13 <- ddG_k13[, Pos_real := Pos + 1]
ddG_k13 <- ddG_k13[, c(1, 20, 23)]
colnames(ddG_k13)[2] <- paste0(colnames(ddG_k13)[2], "_", assay_sele2)

setnames(result_k13, "original_y", paste0("mean_kcal/mol_", assay_sele2))

cat("result_k13列名:\n")
print(names(result_k13))

merge_k13 <- merge(ddG_k13, result_k13, by = c("Pos_real", paste0("mean_kcal/mol_", assay_sele2)), all = FALSE)
cat("merge_k13列名:\n")
print(names(merge_k13))

#### 最终合并 ####
cat("=== 最终合并 ===\n")
cat("merge_raf1列名:\n")
print(names(merge_raf1))
cat("merge_k13列名:\n")
print(names(merge_k13))

# 使用明确的列名进行合并
common_cols <- intersect(names(merge_raf1), names(merge_k13))
cat("共同列:", common_cols, "\n")

# 重命名残差列以避免冲突
setnames(merge_raf1, "residual", "residual_RAF1")
setnames(merge_k13, "residual", "residual_K13")

data_plot <- merge(merge_raf1, merge_k13, 
                   by = c("Pos_real", "id", "Pos"), 
                   all = FALSE)

cat("最终data_plot列名:\n")
print(names(data_plot))

# 检查残差列是否存在
if (!all(c("residual_K13", "residual_RAF1") %in% names(data_plot))) {
  cat("❌ 错误: 残差列不存在！可用列名为:\n")
  print(names(data_plot))
  stop("请检查列名设置")
}

# 定义结合界面位点
K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)

# 过滤结合界面位点
data_plot_filtered <- data_plot[!Pos_real %in% c(K13_Binding_interface_site, RAF1_Binding_interface_site)]

cat("过滤统计:\n")
cat("过滤前数据点:", nrow(data_plot), "\n")
cat("过滤后数据点:", nrow(data_plot_filtered), "\n")
cat("排除的结合界面位点数量:", length(unique(c(K13_Binding_interface_site, RAF1_Binding_interface_site))), "\n")

# 选择需要的列并去掉NA
df <- data_plot_filtered %>%
  select(residual_K13, residual_RAF1) %>%
  filter(!is.na(residual_K13) & !is.na(residual_RAF1))

cat("去除NA后数据点:", nrow(df), "\n")

# 计算相关系数
cor_value <- cor(df$residual_K13, df$residual_RAF1, method = "pearson")
cat("Pearson correlation =", round(cor_value, 3), "\n")

# 绘图
p <- ggplot(df, aes(x = residual_K13, y = residual_RAF1)) +
  geom_point(alpha = 0.6, color = "#75C2F6", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", fill = "#FFB0A5", alpha = 0.3) +
  stat_cor(method = "pearson", 
           label.x = min(df$residual_K13, na.rm = TRUE), 
           label.y = max(df$residual_RAF1, na.rm = TRUE) - 0.5, 
           size = 2.5,
           label.sep = "\n") +
  labs(
    x = "Residual K13 (kcal/mol)",
    y = "Residual RAF1 (kcal/mol)",
    title = "Correlation of Residuals between K13 and RAF1\n(Excluding Binding Interface Sites)"
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
  coord_cartesian(xlim = c(-1.7, 2.3), ylim = c(-2, 2.1))

print(p)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251028/f5a_Correlation_of_Residuals.pdf",
       p, device = cairo_pdf, width = 4, height = 4)
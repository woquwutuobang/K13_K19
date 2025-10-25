library(data.table)
library(ggplot2)
library(dplyr)



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
# 使用示例
# ===============================

# 示例：RAF1 assay
result_raf1 <- fit_distance_exponential(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  wt_aa = "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
  assay_name = "RAF1"
)

# 示例：K27 assay
result_k27 <- fit_distance_exponential(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  wt_aa = "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
  assay_name = "K27"
)


raf1_data <- result_raf1$data
k27_data <- result_k27$data

assay_sele1<-"RAF1"
assay_sele2<-"K27"


ddG_site_raf1<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",wt_aa)
#names(ddG_site_raf1)

colnames(ddG_site_raf1)[3:10] <- paste0(colnames(ddG_site_raf1)[3:10], "_", assay_sele1)

setnames(raf1_data, "original_y",  paste0("mean_", assay_sele1))
#names(raf1_data)

merge_raf1<-merge(ddG_site_raf1,raf1_data,by=c("Pos_real","mean_RAF1"),all=F)

####

ddG_site_k27<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",wt_aa)
#names(ddG_site_k27)

colnames(ddG_site_k27)[3:10] <- paste0(colnames(ddG_site_k27)[3:10], "_", assay_sele2)

setnames(k27_data, "original_y",  paste0("mean_", assay_sele2))
#names(merge_raf1)

merge_k27<-merge(ddG_site_k27,k27_data,by=c("Pos_real","mean_K27"),all=F)


merge_raf1<-merge_raf1[,c(1,14)]
merge_k27<-merge_k27[,c(1,14)]

data_plot<-merge(merge_raf1,merge_k27,by=c("Pos_real"),suffixes = c(paste0("_", assay_sele2), paste0("_", assay_sele1)),all=F)


names(data_plot)





###plot 润色
# 提取数据并去掉 NA
df <- data_plot %>%
  select(residual_K27, residual_RAF1) %>%
  filter(!is.na(residual_K27) & !is.na(residual_RAF1))

# 计算 Pearson 相关系数
cor_value <- cor(df$residual_K27, df$residual_RAF1, method = "pearson")
cat("Pearson correlation =", round(cor_value, 3), "\n")

# 绘图
p <- ggplot(df, aes(x = residual_K27, y = residual_RAF1)) +
  geom_point(alpha = 0.6, color = "#75C2F6", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", fill = "#FFB0A5", alpha = 0.3) +
  stat_cor(method = "pearson", 
           label.x = min(df$residual_K27), 
           label.y = max(df$residual_RAF1) - 0.5, 
           size = 2.5,  # 调整统计文本大小
           label.sep = "\n") +
  labs(
    x = "Residual K27 (kcal/mol)",
    y = "Residual RAF1 (kcal/mol)",
    title = "Correlation of Residuals between K27 and RAF1"
  ) +
  theme_bw(base_size = 8) +  # 增加基础字体大小
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    panel.border = element_rect(color = "black", linewidth = 0.8)
  ) +
  coord_cartesian(xlim = c(-1, 1.6), ylim = c(-1, 1.2))

print(p)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure5/20251025/f5b Correlation of Residuals between K27 and RAF1 per sites .pdf",
       p, device = cairo_pdf, width = 4, height = 4)

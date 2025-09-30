library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# ===============================
# 通用函数：距离效应拟合（指数模型）
# ===============================
plot_energy_distance_decay_expfit <- function(input, assay_sele, anno_file,
                                              x_cols = c("scHAmin_ligand_K13", "scHAmin_ligand_RAF1"),
                                              titles = c("Distance to K13", "Distance to RAF1"),
                                              x_range = c(0, 35), y_range = c(0, 3),
                                              save_path = NULL) {
  
  # 1. 读入数据
  data <- fread(input)
  data <- data[, c(3, 20:22)]
  colnames(data)[2:4] <- paste0(colnames(data)[2:4], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno_final <- merge(anno, data, by = "Pos", all = TRUE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  plot_list <- list()
  
  # 2. 遍历所有 x_col
  for (i in seq_along(x_cols)) {
    xvector <- anno_final[[x_cols[i]]]
    yvector <- abs(anno_final[[y_col]])  # 能量取绝对值
    
    df <- data.frame(x = xvector, y = yvector)
    df <- df[complete.cases(df), ]
    
    # 初始参数优化
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    initial_guess <- c(a = 1, b = -0.1)
    opt_params <- tryCatch(
      optim(initial_guess, residual_sum_of_squares, x = df$x, y = df$y)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    # 拟合指数模型
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), data = df, start = list(a = opt_params[1], b = opt_params[2])),
      error = function(e) NULL
    )
    
    # 拟合结果和注释
    fit_df <- data.frame()
    annotation_text <- NULL
    
    if (!is.null(fit_model)) {
      x_seq <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length.out = 200)
      y_fit <- predict(fit_model, newdata = data.frame(x = x_seq))
      fit_df <- data.frame(x = x_seq, y = y_fit)
      
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
      
      annotation_text <- paste0("a = ", a_val, ", b = ", b_val, "\n", p_text)
    }
    
    # 计算密度数据
    density_data <- density(df$x, na.rm = TRUE)
    density_df <- data.frame(
      x = density_data$x,
      density = density_data$y
    )
    
    # 主图
    p_main <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.1, size = 0.7, color = "#48B3AF") +
      coord_cartesian(xlim = x_range, ylim = y_range) +
      theme_classic(base_size = 13) +
      labs(
        x = titles[i],
        y = paste0("ΔΔG (", assay_sele, " binding)")
      )
    
    if (nrow(fit_df) > 0) {
      p_main <- p_main + geom_line(data = fit_df, aes(x = x, y = y),
                                   inherit.aes = FALSE, color = "#9A3F3F", linewidth = 0.8)
    }
    if (!is.null(annotation_text)) {
      p_main <- p_main + annotate("text",
                                  x = Inf, y = Inf,
                                  hjust = 1.1, vjust = 1.5,
                                  label = annotation_text,
                                  size = 4, color = "black")
    }
    
    # 密度图 - 添加坐标轴信息
    p_density <- ggplot(density_df, aes(x = x, y = density)) +
      geom_area(fill = "#48B3AF", alpha = 0.4) +
      coord_cartesian(xlim = x_range) +
      theme_classic(base_size = 10) +  # 使用 theme_classic 替代 theme_void
      theme(
        plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5),  # 调整边距
        axis.title.x = element_blank(),  # 移除 x 轴标题（主图已有）
        axis.text.x = element_blank(),   # 移除 x 轴刻度文字
        axis.ticks.x = element_blank(),  # 移除 x 轴刻度线
        axis.title.y = element_text(size = 8),  # 减小 y 轴标题字体
        axis.text.y = element_text(size = 7),   # 减小 y 轴刻度文字
        panel.grid = element_blank()     # 确保没有网格线
      ) +
      labs(
        x = NULL,
        y = "Density"
      )
    
    # 组合主图和密度图
    p_combined <- p_main / p_density + 
      plot_layout(heights = c(3, 1))
    
    plot_list[[i]] <- p_combined
  }
  
  # 3. 拼图
  final_plot <- wrap_plots(plot_list, ncol = 2)
  
  # 4. 保存
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, device = cairo_pdf, height = 6, width = 8)
  }
  
  return(final_plot)
}

# ===============================
# ===============================

### RAF1

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/RAF1_expfit.pdf"
)

# 显示图像
print(plot_result)




### SOS1

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay_sele = "SOS1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/SOS1_expfit.pdf"
)

# 显示图像
print(plot_result)




### K55

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay_sele = "K55",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/K55_expfit.pdf"
)

# 显示图像
print(plot_result)


### K27

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay_sele = "K27",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/K27_expfit.pdf"
)

# 显示图像
print(plot_result)




### RALGDS

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay_sele = "RALGDS",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/RALGDS_expfit.pdf"
)

# 显示图像
print(plot_result)


### PIK3CG

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay_sele = "PIK3CG",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/PIK3CG_expfit.pdf"
)

# 显示图像
print(plot_result)


### K13

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/K13_expfit.pdf"
)

# 显示图像
print(plot_result)


### K19

plot_result <- plot_energy_distance_decay_expfit(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  save_path = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/ddG_decay_law/results/20250930/K19_expfit.pdf"
)

# 显示图像
print(plot_result)

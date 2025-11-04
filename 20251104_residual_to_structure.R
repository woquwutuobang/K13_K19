library(data.table)
library(bio3d)

# ===============================
# 函数：计算assay的残差
# ===============================
calculate_residuals_for_assay <- function(
    input_file,
    anno_file,
    wt_aa,
    assay_name
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
    
    # 8. 返回残差数据
    return(df[, c("Pos_real", "residual")])
  }
  
  # ========== 计算残差 ==========
  cat("正在处理", assay_name, "...\n")
  residuals <- fit_distance_exponential(input_file, anno_file, wt_aa, assay_name)
  
  return(residuals)
}

# ===============================
# 函数：将残差绝对值赋给PDB结构的B因子
# ===============================
assign_residuals_to_pdb <- function(
    residual_table,
    input_PDB,
    chain_KRAS = "A",
    Pos_correction = 0,
    output_PDB_file
) {
  
  # 读取输入结构
  input_structure <- read.pdb(input_PDB)
  output_PDB <- input_structure
  
  # 将指定链的所有 B 因子设为 0
  res_indices <- output_PDB$atom$resno[output_PDB$atom$chain == chain_KRAS]
  for (i in unique(res_indices)) {
    output_PDB$atom$b[output_PDB$atom$resno == i & 
                        output_PDB$atom$chain == chain_KRAS] <- 0
  }
  
  # 遍历残差表格，将绝对值写入 B 因子
  for (i in residual_table$Pos_real) {
    val <- residual_table$residual[residual_table$Pos_real == i]
    if (length(val) > 1) val <- unique(val)
    if (length(val) == 1 && !is.na(val)) {
      # 取绝对值
      val <- abs(val)
      output_PDB$atom$b[output_PDB$atom$resno == (i + Pos_correction) & 
                          output_PDB$atom$chain == chain_KRAS] <- val
    }
  }
  
  # 写入新结构文件
  write.pdb(output_PDB, file = output_PDB_file)
  cat("PDB文件已保存到:", output_PDB_file, "\n")
}

# ===============================
# 使用示例
# ===============================

## RAF1
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb"  

# 1. 计算RAF1的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "RAF1"
)

# 2. 将RAF1残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/RAF1_residuals.pdb"
)



## RALGDS
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1lfd.pdb"  

# 1. 计算RALGDS的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "RALGDS"
)

# 2. 将RALGDS残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "B",
  Pos_correction = 200,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/RALGDS_residuals.pdb"
)






## PI3KCG
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1he8.pdb"  

# 1. 计算PI3KCG的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "PI3KCG"
)

# 2. 将PI3KCG残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "B",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/PI3KCG_residuals.pdb"
)





## SOS1--Q
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1nvw.pdb"  

# 1. 计算SOS1的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "SOS1"
)

# 2. 将SOS1残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "Q",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/SOS1_Q_residuals.pdb"
)




## SOS1--R
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1nvw.pdb"  

# 1. 计算SOS1的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "SOS1"
)

# 2. 将SOS1残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "R",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/SOS1_R_residuals.pdb"
)





## K55
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mla.pdb"  

# 1. 计算K55的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "K55"
)

# 2. 将K55残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K55_residuals.pdb"
)




## K27
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mlb.pdb"  

# 1. 计算K27的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "K27"
)

# 2. 将K27残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K27_residuals.pdb"
)





## K13
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h46.pdb"  

# 1. 计算K13的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "K13"
)

# 2. 将K13残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K13_residuals.pdb"
)





## K19
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"


input_pdb <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h47.pdb"  

# 1. 计算K19的残差
residuals <- calculate_residuals_for_assay(
  input_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  anno_file = anno_file,
  wt_aa = wt_aa,
  assay_name = "K19"
)

# 2. 将K19残差赋给PDB结构
assign_residuals_to_pdb(
  residual_table = residuals,
  input_PDB = input_pdb,
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K19_residuals.pdb"
)


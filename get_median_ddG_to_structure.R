library(data.table)
library(krasddpcams)
library(dplyr)
library(bio3d)



###function 1
krasddpcams__get_median_ddG <- function(ddG, assay_sele) {
  # 读取ddG数据
  ddG <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
  
  # 计算每个Pos_real的中位数
  median_ddG <- ddG[Pos_real > 1,
                    .(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE)),
                    by = "Pos_real"]
  
  return(median_ddG)
}


####function 2

krasddpcams__median_ddG_structure <- function(
    median_ddG_table,
    input_PDB,
    chain_KRAS = "A",
    Pos_correction = 0,
    output_PDB_file
) {
  library(bio3d)
  
  # 读取输入结构
  input_structure <- read.pdb(input_PDB)
  output_PDB <- input_structure
  
  # 先将指定链的所有 B 因子设为 0
  res_indices <- output_PDB$atom$resno[output_PDB$atom$chain == chain_KRAS]
  for (i in unique(res_indices)) {
    output_PDB$atom$b[output_PDB$atom$resno == i & 
                        output_PDB$atom$chain == chain_KRAS] <- 0
  }
  
  # 遍历 median_ddG 表格中的每一个残基位置，并写入 B 因子
  for (i in median_ddG_table$Pos_real) {
    val <- median_ddG_table$median_ddG[median_ddG_table$Pos_real == i]
    # 如果存在多个值，取唯一值（避免 length>1 错误）
    if (length(val) > 1) val <- unique(val)
    if (length(val) == 1 && !is.na(val)) {
      output_PDB$atom$b[output_PDB$atom$resno == (i + Pos_correction) & 
                          output_PDB$atom$chain == chain_KRAS] <- val
    }
  }
  
  # 写入新结构文件
  write.pdb(output_PDB, file = output_PDB_file)
}





####调用

##K13
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h46.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/6h46_K13_median_ddG_use_0901_data.pdb"
)




##K19
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h47.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/6h47_K19_median_ddG_use_0901_data.pdb"
)




##RAF1
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/6vjj_RAF1_median_ddG_use_0901_data.pdb"
)





##SOS1-Q
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay_sele = "SOS1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1nvw.pdb",
  chain_KRAS = "Q",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/1nvw_SOS1_median_ddG_use_0901_data.pdb"
)



##SOS1-R
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
  assay_sele = "SOS1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1nvw.pdb",
  chain_KRAS = "R",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/1nvw_SOS1_R_median_ddG_use_0901_data.pdb"
)





##K55
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  assay_sele = "K55"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mla.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/5mla_K55_median_ddG_use_0901_data.pdb"
)







##K27
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay_sele = "K27"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mlb.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/5mlb_K27_median_ddG_use_0901_data.pdb"
)






##RALGDS
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  assay_sele = "RALGDS"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1lfd.pdb",
  chain_KRAS = "B",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/1lfd_RALGDS_median_ddG_use_0901_data.pdb"
)






##PIK3CG
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  assay_sele = "PIK3CG"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1he8.pdb",
  chain_KRAS = "B",
  Pos_correction = 0,
  output_PDB_file = "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG  endow different type structure/results/20250904/1he8_PIK3CG_median_ddG_use_0901_data.pdb"
)

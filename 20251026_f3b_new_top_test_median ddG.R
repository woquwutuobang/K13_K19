library(data.table)
library(krasddpcams)
library(ggplot2)

# Define parameters
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
anno[, Pos_real := Pos]
anno <- anno[, c(3:50)]


# Function to process single assay data - 修改为计算median值
get_ddG_median_data <- function(ddG_path, wt_aa) {
  ddG <- fread(ddG_path)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", wt_codon := substr(id, 1, 1)]
  ddG[id != "WT", mt_codon := substr(id, nchar(id), nchar(id))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]]))
  )
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate median ddG for each position
  median_ddG <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), 
                    .(median = median(`mean_kcal/mol`, na.rm = TRUE),
                      mad = mad(`mean_kcal/mol`, na.rm = TRUE)), 
                    by = "Pos_real"]
  
  median_ddG[, Pos := Pos_real]
  
  return(median_ddG)
}



assay1 = "K13"
assay2 = "RAF1"
ddG_path1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt"
ddG_path2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt"

# Load data for both assays
ddG_data1 <- get_ddG_median_data(ddG_path1, wt_aa)
ddG_data2 <- get_ddG_median_data(ddG_path2, wt_aa)

# Merge data
ddG_merge <- merge(ddG_data1, ddG_data2, 
                   by = "Pos_real", all = FALSE,
                   suffixes = c(paste0("_", assay1), paste0("_", assay2)))

ddG_plot_data <- merge(anno, ddG_merge, by = "Pos_real", all = FALSE)




# 提取列名
x_col <- paste0("scHAmin_ligand_", assay1)
y_col <- paste0("scHAmin_ligand_", assay2)

# 修改列名：mean -> median
median_col1 <- paste0("median_", assay1)
median_col2 <- paste0("median_", assay2)

# 移除NA值
ddG_plot_data <- ddG_plot_data[!is.na(ddG_plot_data[[median_col1]]) & 
                                 !is.na(ddG_plot_data[[median_col2]]), ]

# 计算共同颜色范围
range1 <- range(ddG_plot_data[[median_col1]], na.rm = TRUE)
range2 <- range(ddG_plot_data[[median_col2]], na.rm = TRUE)

# 取两者范围大的那一边作为共同范围
common_min <- min(range1[1], range2[1])
common_max <- max(range1[2], range2[2])
common_limits <- c(common_min, common_max)

cat(assay1, "ΔΔG median range:", range1, "\n")
cat(assay2, "ΔΔG median range:", range2, "\n")
cat("Common color scale limits:", common_limits, "\n")

# 创建第一个图：用assay1的median着色
p1 <- ggplot(ddG_plot_data, aes(x = .data[[x_col]], y = .data[[y_col]], 
                                color = .data[[median_col1]])) +
  geom_point(size = 2, alpha = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", size = 0.3) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray50", size = 0.3) +
  scale_color_gradient2(
    low = "#1B38A6",
    mid = "grey90",
    high = "#F4270C",
    midpoint = 0,
    name = bquote(Delta * Delta * G[b] ~ "(" * .(assay1) * ", kcal/mol)"),
    limits = common_limits
  ) +
  labs(
    #title = paste0("Distance vs ΔΔG: ", assay1),
    x = paste0("Distance to ", assay1, " (Å)"),
    y = paste0("Distance to ", assay2, " (Å)")
  ) +
  coord_cartesian(xlim = c(0, 35), ylim = c(0, 35)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.2, "cm"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5)
  )

# 创建第二个图：用assay2的median着色
p2 <- ggplot(ddG_plot_data, aes(x = .data[[x_col]], y = .data[[y_col]], 
                                color = .data[[median_col2]])) +
  geom_point(size = 2, alpha = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", size = 0.3) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray50", size = 0.3) +
  scale_color_gradient2(
    low = "#1B38A6",
    mid = "grey90",
    high = "#F4270C",
    midpoint = 0,
    name = bquote(Delta * Delta * G[b] ~ "(" * .(assay2) * ", kcal/mol)"),
    limits = common_limits
  ) +
  labs(
    #title = paste0("Distance vs ΔΔG: ", assay2),
    x = paste0("Distance to ", assay1, " (Å)"),
    y = paste0("Distance to ", assay2, " (Å)")
  ) +
  coord_cartesian(xlim = c(0, 35), ylim = c(0, 35)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.2, "cm"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5)
  )

# 显示图形
print(p1)
print(p2)  

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/distance vs K13 VS RAF1 color K13_median.pdf", p1, width = 5, height = 4.2, device = cairo_pdf)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251026/distance vs K13 VS RAF1 color RAF1_median.pdf", p2, width = 5, height = 4.2, device = cairo_pdf)

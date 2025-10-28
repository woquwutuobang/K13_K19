library(data.table)
library(ggplot2)
library(ggnewscale)
library(patchwork)
library(dplyr)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

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

ddG_RAF1 <- ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt", wt_aa)
ddG_K27 <- ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt", wt_aa)

ddG_merged <- merge(
  ddG_RAF1,
  ddG_K27,
  by = c("Pos_real", "codon"),
  all = FALSE,  
  suffixes = c("_RAF1", "_K27")
)

### 结合位点及变构位点及核苷酸交换口袋位点信息
K27_Binding_interface_site <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
GTP_Binding_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
K27_allosteric_site <- c(16, 17, 28, 57, 60, 117, 119,10, 63, 76, 84, 90, 144, 151, 155, 156)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)



### p1---------------all sites
p1<-ggplot(ddG_merged, aes(x = median_RAF1, y = median_K27)) +
  geom_point(color = "grey70", size = 3.5) +
  
  # 纵向误差条 - 不带封顶横线
  geom_errorbar(aes(ymin = median_K27 - sigma_K27, ymax = median_K27 + sigma_K27), 
                width = 0, color = "grey30", size = 0.2) +
  
  # 横向误差条 - 不带封顶横线
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey30", size = 0.2) +
  
  # 固定坐标轴范围
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K27 (kcal/mol)") +
  # title = "Scatter plot of median ΔΔGb: RAF1 vs K13") +
  
  theme_minimal(base_size = 11) +
  theme(
    # 去掉背景网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加坐标轴线
    axis.line = element_line(color = "black", size = 0.5),
    # 可选：添加坐标轴刻度线
    axis.ticks = element_line(color = "black", size = 0.5),
    # X轴刻度值旋转90度
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p1)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251027/scatter plot all sites.pdf", p1,device = cairo_pdf, width = 4, height = 4)



#### p2-------K13 BI
p2 <- ggplot(ddG_merged[Pos_real %in% K27_Binding_interface_site], 
             aes(x = median_RAF1, y = median_K27)) +
  geom_point(color = "#F4270C", size = 3.5) +
  geom_errorbar(aes(ymin = median_K27 - sigma_K27, ymax = median_K27 + sigma_K27), 
                width = 0, color = "grey30", size = 0.2) +
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K27 (kcal/mol)")+
  #title = "K13 Binding Interface Sites") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p2)
#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/scatter plot K13 BI.pdf", p2,device = cairo_pdf, width = 4, height = 4)


### p3---------RAF1 BI
p3 <- ggplot(ddG_merged[Pos_real %in% RAF1_Binding_interface_site], 
             aes(x = median_RAF1, y = median_K27)) +
  geom_point(color = "#1B38A6", size = 3.5) +
  geom_errorbar(aes(ymin = median_K27 - sigma_K27, ymax = median_K27 + sigma_K27), 
                width = 0, color = "grey30", size = 0.2) +
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K27 (kcal/mol)")+
  #title = "RAF1 Binding Interface Sites") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p3)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/scatter plot RAF1 BI.pdf", p3,device = cairo_pdf, width = 4, height = 4)

combined_plot_label <- p1 + p2 + p3 +  
  plot_layout(ncol = 3) 


print(combined_plot_label)

# 保存组合图
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251027/combined_scatter_plots_label_allosites site p1_p2_p3 RAF1 VS K27.pdf", 
       combined_plot_label, device = cairo_pdf, width = 15, height = 5)




core_surface <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv")
core_surface <- core_surface[, Pos_real := Pos]
core_surface_threshold <- 0.25
#surface_threshold_sasa <- 0.26
core_surface <- core_surface %>%
  mutate(type = case_when(
    RASA <= core_surface_threshold ~ "core",
    RASA > core_surface_threshold ~ "surface"
  )) %>%
  filter(!is.na(type))

### 将core/surface信息合并到ddG_merged中
ddG_merged <- merge(ddG_merged, core_surface, by = "Pos_real", all.x = TRUE)

### 移除两个结合界面的位点
remaining_sites <- ddG_merged[!Pos_real %in% c(K27_Binding_interface_site, RAF1_Binding_interface_site)]

# 添加GTP结合口袋标记
remaining_sites[, is_gtp := Pos_real %in% GTP_Binding_pocket]

# 在remaining_sites中添加变构位点分类
remaining_sites[, allosteric_type := "Other"]
remaining_sites[Pos_real %in% K27_allosteric_site & !Pos_real %in% RAF1_allosteric_site, 
                allosteric_type := "K27_allosteric_only"]
remaining_sites[Pos_real %in% RAF1_allosteric_site & !Pos_real %in% K27_allosteric_site, 
                allosteric_type := "RAF1_allosteric_only"]
remaining_sites[Pos_real %in% K27_allosteric_site & Pos_real %in% RAF1_allosteric_site, 
                allosteric_type := "Both_allosteric"]

# 创建颜色映射
allosteric_colors <- c(
  "Other" = "grey70",
  "K27_allosteric_only" = "#FFB0A5",  
  "RAF1_allosteric_only" = "#F4AD0C",  
  "Both_allosteric" = "#C68EFD"        
)

# 函数：计算OR值和创建标签文本
calculate_or_label <- function(positions, k27_allo_sites, raf1_allo_sites) {
  is_k27_allosteric <- positions %in% k27_allo_sites
  is_raf1_allosteric <- positions %in% raf1_allo_sites
  contingency_table <- table(is_k27_allosteric, is_raf1_allosteric)
  fisher_test <- fisher.test(contingency_table)
  
  or_value <- round(fisher_test$estimate, 2)
  if (fisher_test$p.value < 0.001) {
    p_label <- "p < 0.001"
  } else if (fisher_test$p.value < 0.01) {
    p_label <- "p < 0.01" 
  } else if (fisher_test$p.value < 0.05) {
    p_label <- "p < 0.05"
  } else {
    p_label <- paste0("p = ", round(fisher_test$p.value, 3))
  }
  
  return(paste0("OR = ", or_value, "\n", p_label))
}

# 函数：创建散点图
create_scatter_plot <- function(data, title, or_label) {
  ggplot(data, aes(x = median_RAF1, y = median_K27)) +
    # 先画非GTP结合口袋的点（圆形）
    geom_point(data = data[is_gtp == FALSE], 
               aes(color = allosteric_type, shape = "Non-GTP pocket"), 
               size = 3.5, alpha = 1) +
    # 再画GTP结合口袋的点（三角形）
    geom_point(data = data[is_gtp == TRUE], 
               aes(color = allosteric_type, shape = "GTP pocket"), 
               size = 3.5, alpha = 1) +
    # 误差条
    geom_errorbar(aes(ymin = median_K27 - sigma_K27, ymax = median_K27 + sigma_K27), 
                  width = 0, color = "grey30", size = 0.2, alpha = 1) +
    geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                   height = 0, color = "grey30", size = 0.2, alpha = 1) +
    # 颜色标度
    scale_color_manual(values = allosteric_colors,
                       name = "Allosteric Site Type",
                       breaks = c("K27_allosteric_only", "RAF1_allosteric_only", "Both_allosteric", "Other"),
                       labels = c("K27 allosteric only", "RAF1 allosteric only", "Both allosteric", "Other")) +
    # 形状标度
    scale_shape_manual(name = "GTP Binding Pocket",
                       values = c("Non-GTP pocket" = 16, "GTP pocket" = 17),
                       labels = c("GTP pocket" = "Yes", "Non-GTP pocket" = "No")) +
    coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
    labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
         y = "Median ΔΔGb K27 (kcal/mol)",
         title = title) +
    # 添加OR和p值标签
    annotate("text", x = -0.5, y = 1.7, label = or_label, 
             hjust = 0, vjust = 1, size = 3.3, color = "black") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}

# 计算整体的OR值（用于p4）
all_positions <- unique(remaining_sites$Pos_real)
or_label_all <- calculate_or_label(all_positions, K27_allosteric_site, RAF1_allosteric_site)

# p4: 所有位点
p4 <- create_scatter_plot(remaining_sites, "All sites outside interfaces", or_label_all)

# p5: Core位点
core_sites <- remaining_sites[type == "core"]
core_positions <- unique(core_sites$Pos_real)
or_label_core <- calculate_or_label(core_positions, K27_allosteric_site, RAF1_allosteric_site)
p5 <- create_scatter_plot(core_sites, "Core sites outside interfaces", or_label_core)

# p6: Surface位点
surface_sites <- remaining_sites[type == "surface"]
surface_positions <- unique(surface_sites$Pos_real)
or_label_surface <- calculate_or_label(surface_positions, K27_allosteric_site, RAF1_allosteric_site)
p6 <- create_scatter_plot(surface_sites, "Surface sites outside interfaces", or_label_surface)

# 打印各个OR值
cat("All sites OR:", or_label_all, "\n")
cat("Core sites OR:", or_label_core, "\n")
cat("Surface sites OR:", or_label_surface, "\n")

# 打印各个组的位点数量
cat("Core sites count:", length(core_positions), "\n")
cat("Surface sites count:", length(surface_positions), "\n")

print(p4)
print(p5)
print(p6)

# 组合图形
combined_plot <- p4 + p5 + p6 + 
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

# 保存组合图
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251027/core_surface_scatter_plots K27 VS RAF1.pdf", 
       combined_plot, device = cairo_pdf, width = 15, height = 5.5)

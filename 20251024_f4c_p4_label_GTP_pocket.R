library(data.table)
library(ggplot2)
library(ggnewscale)
library(patchwork)


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
  sigma <- ddG[Pos_real > 1, sqrt(1/sum(1/(`std_kcal/mol`^2), na.rm=TRUE)), by=Pos_real]; setnames(sigma, "V1", "sigma") ##通过精度加权平均计算得到的每个氨基酸位置上所有突变测量的合并标准差
  max <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), max(`mean_kcal/mol`), by=Pos_real]; setnames(max, "V1", "max")
  min <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), min(`mean_kcal/mol`), by=Pos_real]; setnames(min, "V1", "min")
  count <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), .N, by=Pos_real]; setnames(count, "N", "count")
  median_ddG <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), median(`mean_kcal/mol`, na.rm=TRUE), by=Pos_real]; setnames(median_ddG, "V1", "median")
  abs_median_ddG <- ddG[Pos_real > 1 & !is.na(`mean_kcal/mol`), median(abs(`mean_kcal/mol`), na.rm=TRUE), by=Pos_real]; setnames(abs_median_ddG, "V1", "abs_median")
  
  output <- Reduce(function(x, y) merge(x, y, by = "Pos_real", all = TRUE),
                   list(codon, mean, abs_mean, sigma, max, min, count, median_ddG, abs_median_ddG))
  return(output)
}


ddG_RAF1<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",wt_aa)
ddG_K13<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",wt_aa)
#names(ddG_RAF1)


ddG_merged <- merge(
  ddG_RAF1,
  ddG_K13,
  by = c("Pos_real", "codon"),
  all = FALSE,  
  suffixes = c("_RAF1", "_K13")
)

###结合位点及变构位点及核苷酸交换口袋位点信息

K13_Binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_Binding_interface_site <- c(21, 25, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
GTP_Binding_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
K13_allosteric_site <- c(15, 16, 17, 35, 145, 10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)



### p1---------------all sites
p1<-ggplot(ddG_merged, aes(x = median_RAF1, y = median_K13)) +
  geom_point(color = "grey40", size = 2) +
  
  # 纵向误差条 - 不带封顶横线
  geom_errorbar(aes(ymin = median_K13 - sigma_K13, ymax = median_K13 + sigma_K13), 
                width = 0, color = "grey90", size = 0.3) +
  
  # 横向误差条 - 不带封顶横线
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey90", size = 0.3) +
  
  # 固定坐标轴范围
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Scatter plot of median ΔΔGb: RAF1 vs K13") +
  
  theme_minimal(base_size = 8) +
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

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/scatter plot all sites.pdf", p1,device = cairo_pdf, width = 4, height = 4)



#### p2-------K13 BI
p2 <- ggplot(ddG_merged[Pos_real %in% K13_Binding_interface_site], 
             aes(x = median_RAF1, y = median_K13)) +
  geom_point(color = "#F4270C", size = 2) +
  geom_errorbar(aes(ymin = median_K13 - sigma_K13, ymax = median_K13 + sigma_K13), 
                width = 0, color = "grey90", size = 0.3) +
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey90", size = 0.3) +
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "K13 Binding Interface Sites") +
  theme_minimal(base_size = 8) +
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
             aes(x = median_RAF1, y = median_K13)) +
  geom_point(color = "#1B38A6", size = 2) +
  geom_errorbar(aes(ymin = median_K13 - sigma_K13, ymax = median_K13 + sigma_K13), 
                width = 0, color = "grey90", size = 0.3) +
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey90", size = 0.3) +
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "RAF1 Binding Interface Sites") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p3)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/scatter plot RAF1 BI.pdf", p3,device = cairo_pdf, width = 4, height = 4)




### p4---------remove 2 BI
remaining_sites <- ddG_merged[!Pos_real %in% c(K13_Binding_interface_site, RAF1_Binding_interface_site)]

# 添加GTP结合口袋标记
remaining_sites[, is_gtp := Pos_real %in% GTP_Binding_pocket]

# Fisher检验计算OR值和p值
all_positions <- unique(remaining_sites$Pos_real)
is_k13_allosteric <- all_positions %in% K13_allosteric_site
is_raf1_allosteric <- all_positions %in% RAF1_allosteric_site
contingency_table <- table(is_k13_allosteric, is_raf1_allosteric)
fisher_test <- fisher.test(contingency_table)

# 提取OR值和p值
or_value <- round(fisher_test$estimate, 3)
p_value <- ifelse(fisher_test$p.value < 0.001, "< 0.001", 
                  ifelse(fisher_test$p.value < 0.01, "< 0.01", 
                         round(fisher_test$p.value, 3)))

# 创建标签文本
label_text <- paste0("OR(K13 allosteric sites vs RAF1 allosteric site) = ", or_value, 
                     "\np-value ", p_value)

p4 <- ggplot(remaining_sites, aes(x = median_RAF1, y = median_K13)) +
  geom_point(aes(color = is_gtp), size = 2) +
  geom_errorbar(aes(ymin = median_K13 - sigma_K13, ymax = median_K13 + sigma_K13), 
                width = 0, color = "grey90", size = 0.3) +
  geom_errorbarh(aes(xmin = median_RAF1 - sigma_RAF1, xmax = median_RAF1 + sigma_RAF1), 
                 height = 0, color = "grey90", size = 0.3) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#F1DD10"),
                     labels = c("FALSE" = "Other", "TRUE" = "GTP Binding Pocket"),
                     name = "Site Type") +
  coord_cartesian(xlim = c(-0.7, 2.3), ylim = c(-0.7, 1.8)) +
  labs(x = "Median ΔΔGb RAF1 (kcal/mol)",
       y = "Median ΔΔGb K13 (kcal/mol)",
       title = "Sites outside the two interfaces") +
  # 添加OR和p值标签
  annotate("text", x = -0.5, y = 1.7, label = label_text, 
           hjust = 0, vjust = 1, size = 2.5, color = "black") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p4)

#ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/scatter plot remove K13_RAF1 BI sites.pdf", p4,device = cairo_pdf, width = 5, height = 4)




# 组合四个图形，ncol=4
combined_plot_nolabel <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 4) 
 

print(combined_plot_nolabel)

# 保存组合图
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/combined_scatter_plots_nolabel_allosites.pdf", 
       combined_plot_nolabel, device = cairo_pdf, width = 17, height = 4)

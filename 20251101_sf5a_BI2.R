# KRAS Allosteric Site Analysis
# This script analyzes the relationship between mutation effects and distance from binding partners

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)
library(dplyr)
library(krasddpcams)

# Wild-type KRAS sequence (residues 2-188)
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Secondary structure elements
secondary_structure <- data.frame(
  xstart = c(3, 15, 38, 51, 67, 77, 87, 109, 127, 139, 148),
  xend = c(9, 24, 44, 57, 73, 84, 104, 115, 136, 143, 166),
  col = c("β1", "α1", "β2", "β3", "α2", "β4", "α3", "β5", "α4", "β6", "α5")
)

# File paths for ddG data
ddG_k13 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt"
assay_k13 <- "K13"
ddG_k19 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt"
assay_k19 <- "K19"

# Load annotation data
anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv")
anno[, Pos_real := Pos]

# 使用与参考代码相同的ddG数据处理函数
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

# 处理K13和K19的数据
ddG_K13 <- ddG_data_process(ddG_k13, wt_aa)
ddG_K19 <- ddG_data_process(ddG_k19, wt_aa)

assay_data <- list(ddG_K13, ddG_K19)
assay_names <- c("K13", "K19")

# 合并数据并添加binding_type
data_plot <- rbindlist(lapply(seq_along(assay_data), function(i) {
  df <- assay_data[[i]]
  assay <- assay_names[i]
  df <- merge(df, anno, by = "Pos_real", all = TRUE)
  df[, assay_name := assay]
  df[, binding_type := "allosteric site"]
  sc_col <- paste0("scHAmin_ligand_", assay)
  if (sc_col %in% names(df)) {
    df[get(sc_col) < 5, binding_type := "binding site"]
  }
  return(df)
}), fill = TRUE)

# 计算阈值（与参考代码相同的逻辑）
calc_threshold <- function(df) {
  df[binding_type == "binding site", sum(abs(median)/sigma^2, na.rm=TRUE)/sum(1/sigma^2, na.rm=TRUE)]
}

merged_list <- split(data_plot, data_plot$assay_name)
threshold1 <- mean(sapply(c("K13","K19"), function(a) calc_threshold(merged_list[[a]])))

cat("Threshold for K13 and K19:", threshold1, "\n")

# 位点分类（与参考代码相同的逻辑）
data_plot[, site_type := NA_character_]
for (assay in assay_names) {
  reg_threshold <- threshold1  # 使用K13/K19的阈值
  idx <- data_plot$assay_name == assay
  df_sub <- data_plot[idx]
  df_sub[, binding_type_gtp_included := binding_type]
  if ("GXPMG_scHAmin_ligand_RAF1" %in% names(df_sub)) {
    df_sub[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  }
  df_sub[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  df_sub[binding_type_gtp_included == "GTP binding site", site_type := "Other GTP pocket site"]
  sc_col <- paste0("scHAmin_ligand_", assay)
  df_sub[binding_type_gtp_included == "GTP binding site" & abs(median) > reg_threshold &
           (!is.null(df_sub[[sc_col]]) & get(sc_col) >= 5),
         site_type := "Allosteric GTP pocket site"]
  df_sub[binding_type_gtp_included == "allosteric site" & abs(median) > reg_threshold,
         site_type := "Major allosteric site"]
  data_plot[idx, binding_type_gtp_included := df_sub$binding_type_gtp_included]
  data_plot[idx, site_type := df_sub$site_type]
  
  # 打印site_type信息
  cat("=== Assay:", assay, "===\n")
  site_counts <- df_sub[, .N, by = site_type]
  print(site_counts)
  
  # 打印具体位置信息
  for (stype in c("Allosteric GTP pocket site", "Major allosteric site", "Binding interface site")) {
    positions <- df_sub[site_type == stype, unique(Pos_real)]
    if (length(positions) > 0) {
      cat(stype, "positions:", paste(sort(positions), collapse = ", "), "\n")
    }
  }
  cat("\n")
}

# 为了与绘图代码兼容，重命名site_type levels
data_plot[, site_type_plot := fcase(
  site_type == "Binding interface site", "binding_interface",
  site_type == "Allosteric GTP pocket site", "allosteric_gtp_pocket", 
  site_type == "Other GTP pocket site", "other_gtp_pocket",
  site_type == "Major allosteric site", "major_allosteric",
  default = "other"
)]

# Map secondary structure information
rects_dt <- as.data.table(secondary_structure)

# Assign beta strands
data_plot[Pos_real >= rects_dt[col == "β1", xstart] & Pos_real <= rects_dt[col == "β1", xend], 
          colors_type := "β1"]
data_plot[Pos_real >= rects_dt[col == "β2", xstart] & Pos_real <= rects_dt[col == "β2", xend], 
          colors_type := "β2"]
data_plot[Pos_real >= rects_dt[col == "β3", xstart] & Pos_real <= rects_dt[col == "β3", xend], 
          colors_type := "β3"]
data_plot[Pos_real >= rects_dt[col == "β4", xstart] & Pos_real <= rects_dt[col == "β4", xend], 
          colors_type := "β4"]
data_plot[Pos_real >= rects_dt[col == "β5", xstart] & Pos_real <= rects_dt[col == "β5", xend], 
          colors_type := "β5"]
data_plot[Pos_real >= rects_dt[col == "β6", xstart] & Pos_real <= rects_dt[col == "β6", xend], 
          colors_type := "β6"]

# Assign alpha helices
data_plot[Pos_real >= rects_dt[col == "α1", xstart] & Pos_real <= rects_dt[col == "α1", xend], 
          colors_type := "α1"]
data_plot[Pos_real >= rects_dt[col == "α2", xstart] & Pos_real <= rects_dt[col == "α2", xend], 
          colors_type := "α2"]
data_plot[Pos_real >= rects_dt[col == "α3", xstart] & Pos_real <= rects_dt[col == "α3", xend], 
          colors_type := "α3"]
data_plot[Pos_real >= rects_dt[col == "α4", xstart] & Pos_real <= rects_dt[col == "α4", xend], 
          colors_type := "α4"]
data_plot[Pos_real >= rects_dt[col == "α5", xstart] & Pos_real <= rects_dt[col == "α5", xend], 
          colors_type := "α5"]

# Classify secondary structure shapes
data_plot[, shape := "other"]
data_plot[colors_type %chin% c("β1", "β2", "β3", "β4", "β5", "β6"), shape := "beta_strand"]
data_plot[colors_type %chin% c("α1", "α2", "α3", "α4", "α5"), shape := "alpha_helix"]

# Filter and factor data
data_plot <- data_plot[Pos_real > 1 & count > 9.5, ]
data_plot[, site_type_plot := factor(site_type_plot,
                                     levels = c("binding_interface", "allosteric_gtp_pocket", 
                                                "other_gtp_pocket", "major_allosteric", "other"))]

data_plot[, assay_name := factor(assay_name, levels = c("K13", "K19"))]

# Calculate distances for each assay
assay_list <- c("K13", "K19")
allosteric_sites <- list()
for (assay_i in assay_list) {
  data_plot[assay_name == assay_i, distance_bp := get(paste0("scHAmin_ligand_", assay_i))]
}

# Create the plot - 使用median来绘图
allosteric_plot <- ggplot() +
  # Points with error bars
  geom_point(data = data_plot,
             aes(x = distance_bp, y = abs_median, color = site_type_plot, shape = shape),
             size = 0.45) +
  geom_pointrange(data = data_plot,
                  aes(x = distance_bp, y = abs_median, ymin = abs_median - sigma, ymax = abs_median + sigma,
                      color = site_type_plot, shape = shape),
                  size = 0.45) +
  
  # Reference lines
  geom_hline(yintercept = threshold1, linetype = 2, size = 0.1) +
  geom_vline(xintercept = 5, linetype = 2, size = 0.1) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.1) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.1) +
  
  # Label allosteric sites
  geom_text_repel(data = data_plot[site_type_plot == "major_allosteric", ],
                  aes(x = distance_bp, y = abs_median, label = Pos_real),
                  nudge_y = 0.05, color = "#F4AD0C", size = 8 * 0.35) +
  geom_text_repel(data = data_plot[site_type_plot == "allosteric_gtp_pocket", ],
                  aes(x = distance_bp, y = abs_median, label = Pos_real),
                  nudge_y = 0.05, color ="#1B38A6", size = 8 * 0.35) +
  
  # Styling
  scale_color_manual(
    values = c("#F4270C", "#1B38A6", "#75C2F6", "#F4AD0C", "gray"),
    labels = c("Binding Interface", "Allosteric GTP Pocket", 
               "Other GTP Pocket", "Major Allosteric", "Others")
  ) +
  xlab("Distance to Binding Partner (Å)") +
  ylab("Median |ΔΔG| (kcal/mol)") +  # 修改y轴标签
  labs(color = "Site Type", shape = "Secondary Structure") +
  facet_wrap(~ assay_name, ncol = 3) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    text = element_text(size = 8),
    legend.position = "right",
    strip.text = element_text(size = 8),
    strip.background = element_rect(colour = "white", fill = "white"),
    panel.spacing = unit(0.2, "mm"),
    legend.text = element_text(size = 8),
    plot.margin = margin(0, 1, 0, 1, "mm"),
    legend.margin = margin(0, 0, 0, -2, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.key.height = unit(4, "mm")
  )

# Display and save plot
print(allosteric_plot)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf5/20251101/K13_K19_allosteric_sites_with_energy_vs_distance_plot_median.pdf",
       plot = allosteric_plot,
       device = cairo_pdf,
       height = 5,
       width = 8,
       dpi = 300)

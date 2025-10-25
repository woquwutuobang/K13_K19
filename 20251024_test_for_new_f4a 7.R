library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

# ======================= 基本设定 =======================
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
core_threshold_sasa <- 0.25
surface_threshold_sasa <- 0.26

anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
anno[, Pos_real := Pos]
core_surface <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/KRAS_WT_166_monomer_get_rasa_20250701_2.csv")
core_surface[, Pos_real := Pos]

# ======================= 定义函数 =======================
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

# ======================= 数据导入 =======================
ddG_RAF1<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",wt_aa)
ddG_SOS1<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",wt_aa)
ddG_K55<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",wt_aa)
ddG_K27<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",wt_aa)
ddG_RALGDS<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",wt_aa)
ddG_PIK3CG<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",wt_aa)
ddG_K13<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",wt_aa)
ddG_K19<-ddG_data_process("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",wt_aa)

assay_data <- list(ddG_RAF1, ddG_SOS1, ddG_K55, ddG_K27, ddG_RALGDS, ddG_PIK3CG, ddG_K13, ddG_K19)
assay_names <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K13", "K19")

# ======================= 数据合并与注释 =======================
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

calc_threshold <- function(df) {
  df[binding_type == "binding site", sum(abs(median)/sigma^2, na.rm=TRUE)/sum(1/sigma^2, na.rm=TRUE)]
}

merged_list <- split(data_plot, data_plot$assay_name)
threshold1 <- mean(sapply(c("K13","K19"), function(a) calc_threshold(merged_list[[a]])))
threshold2 <- mean(sapply(setdiff(assay_names, c("K13","K19")), function(a) calc_threshold(merged_list[[a]])))

data_plot[, site_type := NA_character_]
for (assay in assay_names) {
  reg_threshold <- if (assay %in% c("K13","K19")) threshold1 else threshold2
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
  data_plot[idx, site_type := df_sub$site_type]
}

# ======================= 计算 core/surface =======================
data_plot <- merge(data_plot, core_surface[, .(Pos_real, RASA)], by = "Pos_real", all.x = TRUE)
data_plot[, core_surface_type := fifelse(
  RASA <= core_threshold_sasa, "core",
  fifelse(RASA >= surface_threshold_sasa, "surface", NA_character_)
)]




# =============================
# 绘制KRAS突变能量图 + 距离渐变图（含5Å标注）
# =============================


# ===== 1️⃣ 数据准备 =====
# 提取 RAF1 与 K13 结合距离
pos_mapping <- data_plot[Pos_real >= 2 & Pos_real <= 166,
                         .(Pos_real, scHAmin_ligand_RAF1, scHAmin_ligand_K13)]

# 判断哪个结合位点更近（负号表示靠RAF1）
pos_mapping[, side := ifelse(scHAmin_ligand_RAF1 < scHAmin_ligand_K13, "RAF1", "K13")]

# 定义相对距离（越小越靠近对应界面）
pos_mapping[, rel_distance := ifelse(side == "RAF1",
                                     scHAmin_ligand_RAF1,
                                     scHAmin_ligand_K13)]

# ===== 2️⃣ 生成排序顺序 =====
pos_mapping[, min_distance := pmin(scHAmin_ligand_RAF1, scHAmin_ligand_K13)]
pos_left  <- pos_mapping[side == "RAF1"][order(min_distance)]
pos_right <- pos_mapping[side == "K13"][order(min_distance)]
pos_order <- unique(c(pos_left$Pos_real, rev(pos_right$Pos_real)))
pos_mapping[, Pos_real_factor := factor(Pos_real, levels = pos_order)]
pos_mapping[, color_value := pmin(scHAmin_ligand_RAF1, scHAmin_ligand_K13)]

# ===== 3️⃣ 合并 core/surface 信息 =====
unique_core_surface <- unique(data_plot[, .(Pos_real, core_surface_type)])
pos_mapping <- merge(pos_mapping, unique_core_surface, by = "Pos_real", all.x = TRUE)

# ===== 4️⃣ 计算两端5Å内残基范围 =====
raf1_close <- pos_mapping[scHAmin_ligand_RAF1 < 5 & side == "RAF1", Pos_real]
k13_close  <- pos_mapping[scHAmin_ligand_K13 < 5 & side == "K13", Pos_real]

# 找到这些残基在排序中的位置（用于画虚线）
raf1_cut_pos <- max(which(levels(pos_mapping$Pos_real_factor) %in% raf1_close))
k13_cut_pos  <- min(which(levels(pos_mapping$Pos_real_factor) %in% k13_close))

# ===== 5️⃣ 绘制下方“距离渐变图” =====
p_bottom <- ggplot(pos_mapping, aes(x = Pos_real_factor, y = 0)) +
  geom_tile(aes(fill = color_value), width = 0.9, height = 0.6) +
  # core/surface 点
  geom_point(data = pos_mapping[!is.na(core_surface_type)],
             aes(x = Pos_real_factor, y = 0, color = core_surface_type),
             size = 2, shape = 16) +
  # 添加虚线分隔 (两端5Å内区域)
  geom_vline(xintercept = raf1_cut_pos + 0.5, 
             linetype = "solid", color = "#F4270C", linewidth = 2) +
  geom_vline(xintercept = k13_cut_pos - 0.5, 
             linetype = "solid", color = "#F4270C", linewidth = 2) +
  # 添加注释文字
  annotate("text", x = raf1_cut_pos - 5, y = 0.4, 
           label = "RAF1 (<5Å)", size = 2, hjust = 1) +
  annotate("text", x = k13_cut_pos + 5, y = 0.4, 
           label = "K13 (<5Å)", size = 2, hjust = 0) +
  scale_fill_gradient2(
    low = "white", mid = "black", high = "white",
    midpoint = 0,
    name = "Proximity to Binding Sites"
  ) +
  scale_color_manual(
    values = c("core" = "#C68EFD", "surface" = "#007A20"),
    name = "Core/Surface"
  ) +
  labs(
    x = "Residues ordered by proximity (RAF1 → K13)",
    y = NULL,
    title = "KRAS residues: distance-based ordering with Core/Surface annotation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.box = "horizontal"
  )

# ===== 6️⃣ 绘制上方主热图 =====
top_data <- copy(data_plot)
top_data[, Pos_real_factor := factor(Pos_real, levels = pos_order)]
assay_order <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K13", "K19")
top_data$assay_name <- factor(top_data$assay_name, levels = rev(assay_order))

p_top <- ggplot(top_data, aes(x = Pos_real_factor, y = assay_name, fill = median)) +
  geom_tile(color = "gray20", size = 0.1) +
  scale_fill_gradient2(
    low = "#1B38A6", mid = "gray", high = "#F4270C", midpoint = 0,
    name = expression(median~Delta*Delta*"Gb (kcal/mol)")
  ) +
  geom_point(data = top_data[site_type %in% c("Allosteric GTP pocket site","Major allosteric site")],
             aes(x = Pos_real_factor, y = assay_name, color = site_type),
             size = 1.5) +
  scale_color_manual(values = c("Allosteric GTP pocket site" = "#75C2F6", 
                                "Major allosteric site" = "#FF0066")) +
  labs(
    x = NULL,
    y = "Assay",
    title = "Mutation energetic effects (ΔΔGb)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

# ===== 7️⃣ 上下图拼接与导出 =====
combined_plot <- p_top / p_bottom + plot_layout(heights = c(4, 1))
combined_plot

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/f4a_min_5A_dash_label 3.pdf",
       combined_plot, device = cairo_pdf, width = 45, height = 6)





# =============================
# 绘制KRAS突变能量图 + 距离渐变图（5Å区域矩形框标注）
# =============================


# ===== 1️⃣ 数据准备 =====
pos_mapping <- data_plot[Pos_real >= 2 & Pos_real <= 166,
                         .(Pos_real, scHAmin_ligand_RAF1, scHAmin_ligand_K13)]

pos_mapping[, side := ifelse(scHAmin_ligand_RAF1 < scHAmin_ligand_K13, "RAF1", "K13")]
pos_mapping[, rel_distance := ifelse(side == "RAF1",
                                     scHAmin_ligand_RAF1,
                                     scHAmin_ligand_K13)]

# ===== 2️⃣ 生成排序顺序 =====
pos_mapping[, min_distance := pmin(scHAmin_ligand_RAF1, scHAmin_ligand_K13)]
pos_left  <- pos_mapping[side == "RAF1"][order(min_distance)]
pos_right <- pos_mapping[side == "K13"][order(min_distance)]
pos_order <- unique(c(pos_left$Pos_real, rev(pos_right$Pos_real)))
pos_mapping[, Pos_real_factor := factor(Pos_real, levels = pos_order)]
pos_mapping[, color_value := pmin(scHAmin_ligand_RAF1, scHAmin_ligand_K13)]

# ===== 3️⃣ 合并 core/surface 信息 =====
unique_core_surface <- unique(data_plot[, .(Pos_real, core_surface_type)])
pos_mapping <- merge(pos_mapping, unique_core_surface, by = "Pos_real", all.x = TRUE)

# ===== 4️⃣ 计算两端5Å内残基范围 =====
raf1_close <- pos_mapping[scHAmin_ligand_RAF1 < 5 & side == "RAF1", Pos_real]
k13_close  <- pos_mapping[scHAmin_ligand_K13 < 5 & side == "K13", Pos_real]

raf1_cut_min <- min(which(levels(pos_mapping$Pos_real_factor) %in% raf1_close))
raf1_cut_max <- max(which(levels(pos_mapping$Pos_real_factor) %in% raf1_close))
k13_cut_min  <- min(which(levels(pos_mapping$Pos_real_factor) %in% k13_close))
k13_cut_max  <- max(which(levels(pos_mapping$Pos_real_factor) %in% k13_close))


# ===== 5️⃣ 绘制下方“距离渐变图” =====
p_bottom <- ggplot(pos_mapping, aes(x = Pos_real_factor, y = 0)) +
  geom_tile(aes(fill = color_value), width = 0.9, height = 1.1) +
  
  # Core/surface 注释点
  geom_point(data = pos_mapping[!is.na(core_surface_type)],
             aes(x = Pos_real_factor, y = 0, color = core_surface_type),
             size = 2, shape = 16) +
  
  # ===== 添加方框标出两端 <5Å 区域 =====
annotate("rect",
         xmin = 0.5,
         xmax = raf1_cut_pos + 0.5,
         ymin = -0.5,
         ymax = 0.5,
         color = "#F4270C", fill = NA, linewidth = 1.5, linetype = "solid") +
  annotate("rect",
           xmin = k13_cut_pos - 0.5,
           xmax = length(levels(pos_mapping$Pos_real_factor)) + 0.5,
           ymin = -0.5,
           ymax = 0.5,
           color = "#F4270C", fill = NA, linewidth = 1.5, linetype = "solid") +
  
  # ===== 添加方框文字标签（放在方框正上方） =====
annotate("text",
         x = (raf1_cut_pos + 1) / 2,
         y = 0.65,
         label = "RAF1 (<5Å)",
         size = 2,
         hjust = 0.5
         ) +
  annotate("text",
           x = (k13_cut_pos + length(levels(pos_mapping$Pos_real_factor))) / 2,
           y = 0.65,
           label = "K13 (<5Å)",
           size = 2,
           hjust = 0.5
           ) +
  
  # ===== 颜色与图例设定 =====
scale_fill_gradient2(
  low = "white", mid = "black", high = "white",
  midpoint = 0,
  name = "Proximity to Binding Sites"
) +
  scale_color_manual(
    values = c("core" = "#C68EFD", "surface" = "#007A20"),
    name = "Core/Surface"
  ) +
  
  # ===== 图形主题设置 =====
labs(
  x = "Residues ordered by proximity (RAF1 → K13)",
  y = NULL
  #title = "KRAS residues: distance-based ordering with Core/Surface annotation"
) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.box = "horizontal"
  )

# ===== 6️⃣ 绘制上方主热图 =====
top_data <- copy(data_plot)
top_data[, Pos_real_factor := factor(Pos_real, levels = pos_order)]
assay_order <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K13", "K19")
top_data$assay_name <- factor(top_data$assay_name, levels = rev(assay_order))

p_top <- ggplot(top_data, aes(x = Pos_real_factor, y = assay_name, fill = median)) +
  geom_tile(color = "gray20", size = 0.1) +
  scale_fill_gradient2(
    low = "#1B38A6", mid = "gray", high = "#F4270C", midpoint = 0,
    name = expression(median~Delta*Delta*"Gb (kcal/mol)")
  ) +
  geom_point(data = top_data[site_type %in% c("Allosteric GTP pocket site","Major allosteric site")],
             aes(x = Pos_real_factor, y = assay_name, color = site_type),
             size = 1.5) +
  scale_color_manual(values = c("Allosteric GTP pocket site" = "#75C2F6", 
                                "Major allosteric site" = "#FF0066")) +
  labs(
    x = NULL,
    y = "Assay",
    title = "Mutation energetic effects (ΔΔGb)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

# ===== 7️⃣ 上下图拼接与导出 =====
combined_plot <- p_top / p_bottom + plot_layout(heights = c(4, 1))
combined_plot

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/f4a_min_5A_box 6.pdf",
       combined_plot, device = cairo_pdf, width = 45, height = 6)

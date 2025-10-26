library(ggplot2)
library(reshape2)
library(dplyr)

wt_aa <- unlist(strsplit("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM", ""))
positions <- 1:length(wt_aa)

K13_binding_interface <- c(63,68,87,88,90,91,92,94,95,96,97,98,99,101,102,105,106,107,129,133,136,137,138) 
K19_binding_interface <- c(68,87,88,90,91,92,94,95,97,98,99,101,102,105,107,108,125,129,133,136,137) 
RAF1_binding_interface <- c(21,25,29,31,33,36,37,38,39,40,41,67,71) 
K27_binding_interface <- c(21,24,25,27,29,30,31,32,33,34,35,36,38,39,40,41,43,52,54,67,70,71)

df <- data.frame(
  Position = rep(positions, 4),
  Residue = rep(wt_aa, 4),
  Partner = rep(c("DARPin K13", "DARPin K19", "RAF1-RBD","K27"), each = length(wt_aa)),
  Contact = FALSE
)

# 设置接触位点
df$Contact[df$Partner == "DARPin K13" & df$Position %in% K13_binding_interface] <- TRUE
df$Contact[df$Partner == "DARPin K19" & df$Position %in% K19_binding_interface] <- TRUE
df$Contact[df$Partner == "RAF1-RBD" & df$Position %in% RAF1_binding_interface] <- TRUE
df$Contact[df$Partner == "K27" & df$Position %in% K27_binding_interface] <- TRUE

# 添加是否为共同接触点（出现在2个以上配体中）
all_contacts <- list(K13_binding_interface, K19_binding_interface, RAF1_binding_interface,K27_binding_interface)
contact_freq <- table(unlist(all_contacts))
common_contact_positions <- as.integer(names(contact_freq[contact_freq >= 2]))
df$Common <- df$Position %in% common_contact_positions & df$Contact

# 获取所有接触位点的并集
all_contact_positions <- unique(unlist(all_contacts))

# 提取用于画图的子集
df_plot <- df %>%
  filter(Position %in% all_contact_positions) %>%
  mutate(ResLabel = paste0(Residue, Position))  

# 构建唯一的 Position → ResLabel 映射
residue_labels <- df_plot %>%
  distinct(Position, Residue) %>%
  arrange(Position) %>%
  mutate(ResLabel = paste0(Residue, Position)) %>%
  pull(ResLabel)

# 然后给 ResLabel 设置 factor 顺序（避免重复）
df_plot$ResLabel <- paste0(df_plot$Residue, df_plot$Position)
df_plot$ResLabel <- factor(df_plot$ResLabel, levels = residue_labels)

# 重新设置 Fill 列（你原本的逻辑）
df_plot$Fill <- ifelse(df_plot$Common, "Common contact", 
                       ifelse(df_plot$Contact, "Contact", "None"))
df_plot$Fill <- factor(df_plot$Fill, levels = c("None", "Contact", "Common contact"))

# 关键修改：设置Partner因子的顺序，确保从上到下为RAF1/K27/K13/K19
df_plot$Partner <- factor(df_plot$Partner, 
                          levels = c("RAF1-RBD", "K27", "DARPin K13", "DARPin K19"))

# ========== 修改部分：创建BI数据，将BI1和BI2放在同一行 ==========
# 定义BI1 (RAF1/K27结合位点) 和 BI2 (K13/K19结合位点)
BI1_positions <- unique(c(RAF1_binding_interface, K27_binding_interface))
BI2_positions <- unique(c(K13_binding_interface, K19_binding_interface))

# 创建BI数据框 - 只创建一行
bi_df <- data.frame(
  Position = all_contact_positions,
  ResLabel = paste0(wt_aa[all_contact_positions], all_contact_positions),
  Partner = "BI",
  Fill = "None"
)

# 设置填充颜色：如果同时属于BI1和BI2，使用混合颜色；否则分别使用BI1或BI2颜色
bi_df$Fill <- ifelse(
  bi_df$Position %in% BI1_positions & bi_df$Position %in% BI2_positions,
  "Both",
  ifelse(
    bi_df$Position %in% BI1_positions,
    "BI1",
    ifelse(
      bi_df$Position %in% BI2_positions,
      "BI2",
      "None"
    )
  )
)

# 设置ResLabel为因子
bi_df$ResLabel <- factor(bi_df$ResLabel, levels = residue_labels)

# 合并原始数据和BI数据
combined_df <- bind_rows(df_plot, bi_df)
combined_df$Partner <- factor(combined_df$Partner, 
                              levels = c("RAF1-RBD", "K27", "DARPin K13", "DARPin K19", "BI"))

# 更新填充颜色的级别和颜色映射
combined_df$Fill <- factor(combined_df$Fill, 
                           levels = c("None", "Contact", "Common contact", "BI1", "BI2", "Both"))

# ========== 绘制图形 ==========
ggplot(combined_df, aes(x = ResLabel, y = Partner, fill = Fill)) +
  geom_tile(color = "grey90", width = 1, height = 1) +  # 正方形格子
  scale_fill_manual(
    values = c(
      "None" = alpha("white", 1),           
      "Contact" = alpha("grey70", 0.8),    
      "Common contact" = alpha("black", 0.8), 
      "BI1" = alpha("#F4270C", 0.8),          # 红色 - BI1 RAF1/K27 并集
      "BI2" = alpha("#1B38A6", 0.8),          # 蓝色 - BI2 K13/K19 并集
      "Both" = alpha("#C68EFD", 1)          # 紫色 - 同时属于BI1和BI2
    ),
    name = "Contact Type",
    labels = c("None", "Contact", "Common contact", "BI1 (RAF1/K27)", "BI2 (K13/K19)", "Both")
  ) +
  scale_y_discrete(limits = rev(levels(combined_df$Partner))) +  # 反转Y轴顺序
  coord_fixed(ratio = 1) +  # 保持正方形
  theme_minimal(base_size = 8) +
  labs(x = "Residue", y = "Binding Partner", fill = "") +  # 添加图例标题
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251025/Common_and_unique_structural_contacts_with_BI 2.pdf", 
       width = 10, height = 5.5, dpi = 300)  # 调整高度

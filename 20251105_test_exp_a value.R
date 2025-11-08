library(bio3d)

calculate_interface_contacts_simple <- function(mapping_df, cutoff = 5) {
  results <- data.frame(
    PDB = character(),
    Binder = character(),
    Interface_contacts = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(mapping_df)) {
    pdb_file <- mapping_df$pdb[i]
    kras_chain <- mapping_df$kras_chain[i]
    binder_chain <- mapping_df$binder_chain[i]
    
    cat("Processing:", pdb_file, "\n")
    
    pdb <- read.pdb(pdb_file, rm.alt = TRUE)
    
    # 多model只取第一个
    if (!is.null(pdb$mdl)) pdb <- trim.pdb(pdb, model = 1)
    
    # 自动选择binder链
    if (is.na(binder_chain) | binder_chain == "") {
      binder_chain <- unique(pdb$atom$chain)[unique(pdb$atom$chain) != kras_chain][1]
    }
    
    sel_kras <- atom.select(pdb, chain = kras_chain)
    sel_binder <- atom.select(pdb, chain = binder_chain)
    
    if (length(sel_kras$xyz) == 0 | length(sel_binder$xyz) == 0) {
      warning(paste("⚠️ Chain not found for file:", pdb_file))
      next
    }
    
    # 提取坐标矩阵 (n x 3)
    xyz_kras <- matrix(pdb$xyz[sel_kras$xyz], ncol=3, byrow=TRUE)
    xyz_binder <- matrix(pdb$xyz[sel_binder$xyz], ncol=3, byrow=TRUE)
    
    # 计算距离矩阵 (KRAS原子 x binder原子)
    dists <- as.matrix(dist(rbind(xyz_kras, xyz_binder)))
    n_kras <- nrow(xyz_kras)
    n_binder <- nrow(xyz_binder)
    
    # KRAS x Binder 子矩阵
    dists_sub <- dists[1:n_kras, (n_kras+1):(n_kras+n_binder)]
    
    # 判断接触
    contact <- dists_sub < cutoff
    
    # 接触的KRAS残基数
    kras_residues <- pdb$atom$resno[sel_kras$atom]
    contact_residues <- unique(kras_residues[rowSums(contact) > 0])
    
    results <- rbind(results, data.frame(
      PDB = basename(pdb_file),
      Binder = binder_chain,
      Interface_contacts = length(contact_residues)
    ))
  }
  
  return(results)
}



# =====================
# Example usage
# =====================
mapping <- data.frame(
  pdb = c("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h46.pdb", 
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h47.pdb", 
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb",
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1lfd.pdb",
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1he8.pdb",
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/1nvw.pdb",
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mla.pdb",
          "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/5mlb.pdb"),
  kras_chain = c("A", "A", "A","B","B","Q","A","A"),
  binder_chain = c("B", "B", "B","A","A","S","B","B")
)

contacts_df <- calculate_interface_contacts_simple(mapping, cutoff=5)
print(contacts_df)




#assay_names<-c("K13","K19","RAF1","RALGDS","PI3KCG","SOS1","K55","K27")
#interface_contacts_count<-c(40,40,34,30,19,36,31,25)
#a_value<-c(0.365,0.298,0.753,0.979,0.909,0.731,0.891,0.510)


assay_names <- c("K13", "K19", "RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27")
interface_contacts_count <- c(40, 40, 34, 30, 19, 36, 31, 25)
a_value <- c(0.365, 0.298, 0.753, 0.979, 0.909, 0.731, 0.891, 0.510)

# 创建数据框
data <- data.frame(
  assay = assay_names,
  contacts = interface_contacts_count,
  a_value = a_value
)

# 计算相关性
cor_test <- cor.test(data$contacts, data$a_value, method = "pearson")
cor_value <- cor_test$estimate
p_value <- cor_test$p.value

# 添加显著性标记
significance <- ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(p_value < 0.05, "*", "NS")))

# 创建散点图
library(ggplot2)

p <- ggplot(data, aes(x = contacts, y = a_value)) +
  geom_point(size = 4, color = "steelblue", alpha = 0.7) +
  # 添加assay名称标签
  geom_text(aes(label = assay), 
            vjust = -0.8, hjust = 0.5, size = 3.5, 
            color = "darkred", fontface = "bold") +
  # 添加趋势线
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              linetype = "dashed", alpha = 0.2) +
  # 添加相关性信息
  annotate("text", x = min(data$contacts), y = max(data$a_value),
           label = paste0("r = ", round(cor_value, 3), 
                          "\np = ", round(p_value, 4),
                          " ", significance),
           hjust = 0, vjust = 1, size = 4.5,
           color = "darkgreen", fontface = "bold") +
  labs(
    title = "Correlation between Interface Contacts Count and a-value",
    x = "Interface Contacts Count",
    y = "a-value",
    caption = paste("Pearson correlation analysis, n =", nrow(data))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(face = "italic")
  )

print(p)

# 保存图片
ggsave("interface_contacts_a_value_correlation.pdf", 
       p, width = 10, height = 8, device = cairo_pdf)
ggsave("interface_contacts_a_value_correlation.png", 
       p, width = 10, height = 8, dpi = 300)

# 输出详细的相关性结果
cat("=== 相关性分析结果 ===\n")
cat("变量: Interface Contacts Count vs a-value\n")
cat("样本数量: ", nrow(data), "\n")
cat("Pearson相关系数 (r): ", round(cor_value, 4), "\n")
cat("p值: ", format.pval(p_value, digits = 4), "\n")
cat("显著性: ", significance, "\n")
cat("95%置信区间: [", 
    round(cor_test$conf.int[1], 4), ", ", 
    round(cor_test$conf.int[2], 4), "]\n")

# 显示数据摘要
cat("\n=== 数据摘要 ===\n")
print(data)

# 线性回归模型详情
cat("\n=== 线性回归模型 ===\n")
lm_model <- lm(a_value ~ contacts, data = data)
summary(lm_model)





################========================================
## 接触面积
#####==================================================


# 数据输入
assay_names <- c("K13", "K19", "RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27")
interface_contacts_area <- c(1916.95, 2121.42, 2228.88, 1331.40, 1918.89, 3264.04, 2269.17, 2416.52)
a_value <- c(0.365, 0.298, 0.753, 0.979, 0.909, 0.731, 0.891, 0.510)

# 创建数据框
data <- data.frame(
  assay = assay_names,
  area = interface_contacts_area,
  a_value = a_value
)

# 计算相关性
cor_test <- cor.test(data$area, data$a_value, method = "pearson")
cor_value <- cor_test$estimate
p_value <- cor_test$p.value

# 添加显著性标记
significance <- ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(p_value < 0.05, "*", "NS")))

# 创建散点图
library(ggplot2)

p <- ggplot(data, aes(x = area, y = a_value)) +
  geom_point(size = 4, color = "steelblue", alpha = 0.7) +
  # 添加assay名称标签
  geom_text(aes(label = assay), 
            vjust = -0.8, hjust = 0.5, size = 3.5, 
            color = "darkred", fontface = "bold") +
  # 添加趋势线
  geom_smooth(method = "lm", se = TRUE, color = "red", 
              linetype = "dashed", alpha = 0.2) +
  # 添加相关性信息
  annotate("text", x = min(data$area), y = max(data$a_value),
           label = paste0("r = ", round(cor_value, 3), 
                          "\np = ", round(p_value, 4),
                          " ", significance),
           hjust = 0, vjust = 1, size = 4.5,
           color = "darkgreen", fontface = "bold") +
  labs(
    title = "Correlation between Interface Contacts Area and a-value",
    x = "Interface Contacts Area (Å²)",
    y = "a-value",
    caption = paste("Pearson correlation analysis, n =", nrow(data))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(face = "italic")
  )

print(p)

# 保存图片
ggsave("interface_contacts_area_a_value_correlation.pdf", 
       p, width = 10, height = 8, device = cairo_pdf)
ggsave("interface_contacts_area_a_value_correlation.png", 
       p, width = 10, height = 8, dpi = 300)

# 输出详细的相关性结果
cat("=== 相关性分析结果 ===\n")
cat("变量: Interface Contacts Area vs a-value\n")
cat("样本数量: ", nrow(data), "\n")
cat("Pearson相关系数 (r): ", round(cor_value, 4), "\n")
cat("p值: ", format.pval(p_value, digits = 4), "\n")
cat("显著性: ", significance, "\n")
cat("95%置信区间: [", 
    round(cor_test$conf.int[1], 4), ", ", 
    round(cor_test$conf.int[2], 4), "]\n")

# 显示数据摘要
cat("\n=== 数据摘要 ===\n")
print(data)

# 线性回归模型详情
cat("\n=== 线性回归模型 ===\n")
lm_model <- lm(a_value ~ area, data = data)
summary(lm_model)


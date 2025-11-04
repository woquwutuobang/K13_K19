

library(data.table)

wt_aa <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Secondary structure information from the paper
rects_sheet <- data.frame(xstart = c(3,38,51,77,109,139),
                          xend = c(9,44,57,84,115,143),
                          col = c("β1","β2","β3","β4","β5","β6"))

rects_helix <- data.frame(xstart = c(15,67,87,127,148),
                          xend = c(24,73,104,136,166),
                          col = c("α1","α2","α3","α4","α5"))

## common allosteric sites
K13_RAF1_common_allosteric_sites <- c(10,15,16,17,35,55,77,79,145,159,163)

## unique allosteric sites
K13_unique_allosteric_sites <- c(19,21,24,53,82,93,151)
RAF1_unique_allosteric_sites <- c(20,28,32,34,54,57,58,60,112,113,114,134,144,146,156)

# 函数：判断位点属于α螺旋还是β折叠
classify_secondary_structure <- function(sites, helix_ranges, sheet_ranges) {
  alpha_sites <- c()
  beta_sites <- c()
  other_sites <- c()
  
  for (site in sites) {
    in_helix <- FALSE
    in_sheet <- FALSE
    
    # 检查是否在α螺旋中
    for (i in 1:nrow(helix_ranges)) {
      if (site >= helix_ranges$xstart[i] && site <= helix_ranges$xend[i]) {
        alpha_sites <- c(alpha_sites, site)
        in_helix <- TRUE
        break
      }
    }
    
    # 检查是否在β折叠中
    for (i in 1:nrow(sheet_ranges)) {
      if (site >= sheet_ranges$xstart[i] && site <= sheet_ranges$xend[i]) {
        beta_sites <- c(beta_sites, site)
        in_sheet <- TRUE
        break
      }
    }
    
    # 如果既不在α螺旋也不在β折叠中
    if (!in_helix && !in_sheet) {
      other_sites <- c(other_sites, site)
    }
  }
  
  return(list(alpha = alpha_sites, beta = beta_sites, other = other_sites))
}

# 分类所有位点
common_classified <- classify_secondary_structure(K13_RAF1_common_allosteric_sites, rects_helix, rects_sheet)
k13_unique_classified <- classify_secondary_structure(K13_unique_allosteric_sites, rects_helix, rects_sheet)
raf1_unique_classified <- classify_secondary_structure(RAF1_unique_allosteric_sites, rects_helix, rects_sheet)

# 输出结果
cat("=== K13/RAF1 共有变构位点分布 ===\n")
cat("α螺旋位点 (", length(common_classified$alpha), "个): ", paste(common_classified$alpha, collapse = ", "), "\n")
cat("β折叠位点 (", length(common_classified$beta), "个): ", paste(common_classified$beta, collapse = ", "), "\n")
cat("其他结构位点 (", length(common_classified$other), "个): ", paste(common_classified$other, collapse = ", "), "\n")
cat("总计: ", length(K13_RAF1_common_allosteric_sites), "个位点\n\n")

cat("=== K13 特有变构位点分布 ===\n")
cat("α螺旋位点 (", length(k13_unique_classified$alpha), "个): ", paste(k13_unique_classified$alpha, collapse = ", "), "\n")
cat("β折叠位点 (", length(k13_unique_classified$beta), "个): ", paste(k13_unique_classified$beta, collapse = ", "), "\n")
cat("其他结构位点 (", length(k13_unique_classified$other), "个): ", paste(k13_unique_classified$other, collapse = ", "), "\n")
cat("总计: ", length(K13_unique_allosteric_sites), "个位点\n\n")

cat("=== RAF1 特有变构位点分布 ===\n")
cat("α螺旋位点 (", length(raf1_unique_classified$alpha), "个): ", paste(raf1_unique_classified$alpha, collapse = ", "), "\n")
cat("β折叠位点 (", length(raf1_unique_classified$beta), "个): ", paste(raf1_unique_classified$beta, collapse = ", "), "\n")
cat("其他结构位点 (", length(raf1_unique_classified$other), "个): ", paste(raf1_unique_classified$other, collapse = ", "), "\n")
cat("总计: ", length(RAF1_unique_allosteric_sites), "个位点\n\n")

# 创建汇总表格
summary_table <- data.frame(
  Category = c("K13/RAF1共有", "K13特有", "RAF1特有"),
  Total = c(length(K13_RAF1_common_allosteric_sites), 
            length(K13_unique_allosteric_sites), 
            length(RAF1_unique_allosteric_sites)),
  Alpha_Helix = c(length(common_classified$alpha), 
                  length(k13_unique_classified$alpha), 
                  length(raf1_unique_classified$alpha)),
  Beta_Sheet = c(length(common_classified$beta), 
                 length(k13_unique_classified$beta), 
                 length(raf1_unique_classified$beta)),
  Other = c(length(common_classified$other), 
            length(k13_unique_classified$other), 
            length(raf1_unique_classified$other))
)

cat("=== 汇总表格 ===\n")
print(summary_table)

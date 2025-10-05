library(krasddpcams)

# Function to merge predicted and observed fitness data for 5 experimental blocks
# This function standardizes fitness data across multiple experimental blocks and assays
krasddpcams__merge_ddGb_ob_pre_fitness_5blocks <- function(
    prediction = "path/to/predicted_phenotypes_all.txt",
    block1_dimsum_df = "path/to/block1_fitness.RData",
    block2_dimsum_df = "path/to/block2_fitness.RData", 
    block3_dimsum_df = "path/to/block3_fitness.RData",
    block4_dimsum_df = "path/to/block4_fitness.RData",
    block5_dimsum_df = "path/to/block5_fitness.RData",
    assay_sele = "RAF1", 
    wt_aa_input = wt_aa) {
  
  # Read prediction data
  pre <- fread(prediction)
  pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
  
  # Load all 5 block datasets
  load(block1_dimsum_df)
  block1 <- as.data.table(all_variants)
  load(block2_dimsum_df)
  block2 <- as.data.table(all_variants)
  load(block3_dimsum_df)
  block3 <- as.data.table(all_variants)
  load(block4_dimsum_df)
  block4 <- as.data.table(all_variants)
  load(block5_dimsum_df)
  block5 <- as.data.table(all_variants)
  
  # Merge all block data
  data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3,
                           block4 = block4, block5 = block5, 
                           idcol = "block", fill = TRUE)
  
  # Calculate weighted fitness intermediates
  data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness/(data_before_nor$sigma)^2
  data_before_nor$one_over_fitness_sigmasquared <- 1/(data_before_nor$sigma)^2
  
  # Calculate stop codon fitness for each block (weighted average)
  calculate_stop_fitness <- function(block_name) {
    dead_fitness <- data_before_nor[STOP == TRUE & block == block_name, ]
    stop_fitness <- sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(stop_fitness)
  }
  
  stop1_fitness <- calculate_stop_fitness("block1")
  stop2_fitness <- calculate_stop_fitness("block2")
  stop3_fitness <- calculate_stop_fitness("block3")
  stop4_fitness <- calculate_stop_fitness("block4")
  stop5_fitness <- calculate_stop_fitness("block5")
  
  # Calculate wild-type fitness for each block (weighted average)
  calculate_wt_fitness <- function(block_name) {
    wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
    wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(wt_fitness_block$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(wt_fitness)
  }
  
  wt1_fitness <- calculate_wt_fitness("block1")
  wt2_fitness <- calculate_wt_fitness("block2")
  wt3_fitness <- calculate_wt_fitness("block3")
  wt4_fitness <- calculate_wt_fitness("block4")
  wt5_fitness <- calculate_wt_fitness("block5")
  
  # Create scaling data for normalization
  scaling_data_fitness <- data.frame(
    block1 = c(stop1_fitness, wt1_fitness),
    block2 = c(stop2_fitness, wt2_fitness),
    block3 = c(stop3_fitness, wt3_fitness),
    block4 = c(stop4_fitness, wt4_fitness),
    block5 = c(stop5_fitness, wt5_fitness)
  )
  
  # Calculate scaling parameters relative to block1
  calculate_scaling_params <- function(target_block) {
    lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
    return(list(
      slope = lm_model$coefficients[[2]],
      intercept = lm_model$coefficients[[1]]
    ))
  }
  
  params_block2 <- calculate_scaling_params("block2")
  params_block3 <- calculate_scaling_params("block3")
  params_block4 <- calculate_scaling_params("block4")
  params_block5 <- calculate_scaling_params("block5")
  
  d2 <- params_block2$slope
  e2 <- params_block2$intercept
  d3 <- params_block3$slope
  e3 <- params_block3$intercept
  d4 <- params_block4$slope
  e4 <- params_block4$intercept
  d5 <- params_block5$slope
  e5 <- params_block5$intercept
  
  # Process prediction data
  pre_nor <- pre_pos
  
  # Extract predicted fitness
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # Assay to phenotype mapping
  assay_sele_df <- data.table(assay = c("K13", "K19", "K27", "K55", "PI3", "RAF1", "RAL", "SOS"), 
                              phenotype_base = c(6, 11, 16, 21, 26, 31, 36, 41))
  
  # Validate assay selection
  if (!assay_sele %in% assay_sele_df$assay) {
    stop("Assay '", assay_sele, "' not found. Available assays: ", 
         paste(assay_sele_df$assay, collapse = ", "))
  }
  
  # Get base phenotype number
  base_pheno <- assay_sele_df[assay == assay_sele, phenotype_base]
  
  # Apply normalization transformations for 5 blocks
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_mean_fitness, mean)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
  pre_nor[phenotype == base_pheno + 3, `:=`(pre_nor_mean_fitness, mean * d4 + e4)]
  pre_nor[phenotype == base_pheno + 4, `:=`(pre_nor_mean_fitness, mean * d5 + e5)]
  
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness_sigma, std)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness_sigma, std * d2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness_sigma, std * d3)]
  pre_nor[phenotype == base_pheno + 3, `:=`(pre_nor_fitness_sigma, std * d4)]
  pre_nor[phenotype == base_pheno + 4, `:=`(pre_nor_fitness_sigma, std * d5)]
  
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness, fitness)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness, fitness * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness, fitness * d3 + e3)]
  pre_nor[phenotype == base_pheno + 3, `:=`(ob_nor_fitness, fitness * d4 + e4)]
  pre_nor[phenotype == base_pheno + 4, `:=`(ob_nor_fitness, fitness * d5 + e5)]
  
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness_sigma, sigma)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness_sigma, sigma * d2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness_sigma, sigma * d3)]
  pre_nor[phenotype == base_pheno + 3, `:=`(ob_nor_fitness_sigma, sigma * d4)]
  pre_nor[phenotype == base_pheno + 4, `:=`(ob_nor_fitness_sigma, sigma * d5)]
  
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness, predicted_fitness)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
  pre_nor[phenotype == base_pheno + 3, `:=`(pre_nor_fitness, predicted_fitness * d4 + e4)]
  pre_nor[phenotype == base_pheno + 4, `:=`(pre_nor_fitness, predicted_fitness * d5 + e5)]
  
  return(pre_nor)
}

# Wild-type KRAS sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# RAF1: Merge predicted and observed fitness data
RAF1_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  block1_dimsum_df = "path/to/RAF_block1_fitness.RData",
  block2_dimsum_df = "path/to/RAF_block2_fitness.RData",
  block3_dimsum_df = "path/to/RAF_block2_2_fitness.RData", 
  block4_dimsum_df = "path/to/RAF_block3_fitness.RData",
  block5_dimsum_df = "path/to/RAF_block3_2_fitness.RData",
  assay_sele = "RAF1",
  wt_aa_input = wt_aa)

# Generate and save RAF1 fitness comparison plots for each block
krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 31)
ggplot2::ggsave("path/to/RAF1_fitness_block1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 32)
ggplot2::ggsave("path/to/RAF1_fitness_block2_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 33)
ggplot2::ggsave("path/to/RAF1_fitness_block2_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 34)
ggplot2::ggsave("path/to/RAF1_fitness_block3_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF1_pre_ob_fitness, phenotypen = 35)
ggplot2::ggsave("path/to/RAF1_fitness_block3_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

# K13: Merge predicted and observed fitness data  
K13_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  block1_dimsum_df = "path/to/K13_block1_fitness.RData",
  block2_dimsum_df = "path/to/K13_block2_fitness.RData",
  block3_dimsum_df = "path/to/K13_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/K13_block3_fitness.RData", 
  block5_dimsum_df = "path/to/K13_block3_2_fitness.RData",
  assay_sele = "K13",
  wt_aa_input = wt_aa)

# Generate and save K13 fitness comparison plots for each block
krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K13_pre_ob_fitness, phenotypen = 6)
ggplot2::ggsave("path/to/K13_fitness_block1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K13_pre_ob_fitness, phenotypen = 7)
ggplot2::ggsave("path/to/K13_fitness_block2_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K13_pre_ob_fitness, phenotypen = 8)
ggplot2::ggsave("path/to/K13_fitness_block2_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K13_pre_ob_fitness, phenotypen = 9)
ggplot2::ggsave("path/to/K13_fitness_block3_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K13_pre_ob_fitness, phenotypen = 10)
ggplot2::ggsave("path/to/K13_fitness_block3_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

# K19: Merge predicted and observed fitness data
K19_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt", 
  block1_dimsum_df = "path/to/K19_block1_fitness.RData",
  block2_dimsum_df = "path/to/K19_block2_fitness.RData",
  block3_dimsum_df = "path/to/K19_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/K19_block3_fitness.RData",
  block5_dimsum_df = "path/to/K19_block3_2_fitness.RData",
  assay_sele = "K19",
  wt_aa_input = wt_aa)

# Generate and save K19 fitness comparison plots for each block
krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K19_pre_ob_fitness, phenotypen = 11)
ggplot2::ggsave("path/to/K19_fitness_block1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K19_pre_ob_fitness, phenotypen = 12)
ggplot2::ggsave("path/to/K19_fitness_block2_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K19_pre_ob_fitness, phenotypen = 13)
ggplot2::ggsave("path/to/K19_fitness_block2_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K19_pre_ob_fitness, phenotypen = 14)
ggplot2::ggsave("path/to/K19_fitness_block3_1.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")

krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = K19_pre_ob_fitness, phenotypen = 15)
ggplot2::ggsave("path/to/K19_fitness_block3_2.pdf", device = cairo_pdf, height = 45, width = 60, units = "mm")
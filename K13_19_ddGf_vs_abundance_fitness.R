library(data.table)
library(krasddpcams)

# Function to merge folding free energy (ddGf) and fitness data for 5 experimental blocks
# This function standardizes abundance fitness data and integrates with folding energy predictions
krasddpcams__merge_ddGf_fitness_5blocks <- function(
    prediction = "path/to/predicted_phenotypes_all.txt",
    folding_ddG = "path/to/folding_weights.txt", 
    block1_dimsum_df = "path/to/block1_fitness.RData",
    block2_dimsum_df = "path/to/block2_fitness.RData", 
    block3_dimsum_df = "path/to/block3_fitness.RData",
    block4_dimsum_df = "path/to/block4_fitness.RData",
    block5_dimsum_df = "path/to/block5_fitness.RData",
    wt_aa_input = wt_aa) {
  
  # Read prediction and folding free energy data
  pre <- fread(prediction)
  folding_ddG <- fread(folding_ddG)
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
  
  # Extract predicted fitness (adjust indices based on actual data structure)
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # Extract additive traits (adjust indices based on actual data structure)
  extract_additive_trait0 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2 - 1])
  }
  extract_additive_trait1 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2])
  }
  
  pre_nor$additive_trait0 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait0)
  pre_nor$additive_trait0 <- as.numeric(pre_nor$additive_trait0)
  pre_nor$additive_trait1 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait1)
  pre_nor$additive_trait1 <- as.numeric(pre_nor$additive_trait1)
  pre_nor[, `:=`(additive_trait, additive_trait0 + additive_trait1)]
  
  # Apply normalization transformations for 5 blocks
  pre_nor[phenotype == 1, `:=`(pre_nor_mean_fitness, mean)]
  pre_nor[phenotype == 2, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_mean_fitness, mean * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_mean_fitness, mean * d5 + e5)]
  
  pre_nor[phenotype == 1, `:=`(pre_nor_fitness_sigma, std)]
  pre_nor[phenotype == 2, `:=`(pre_nor_fitness_sigma, std * d2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_fitness_sigma, std * d3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_fitness_sigma, std * d4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_fitness_sigma, std * d5)]
  
  pre_nor[phenotype == 1, `:=`(ob_nor_fitness, fitness)]
  pre_nor[phenotype == 2, `:=`(ob_nor_fitness, fitness * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(ob_nor_fitness, fitness * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(ob_nor_fitness, fitness * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(ob_nor_fitness, fitness * d5 + e5)]
  
  pre_nor[phenotype == 1, `:=`(ob_nor_fitness_sigma, sigma)]
  pre_nor[phenotype == 2, `:=`(ob_nor_fitness_sigma, sigma * d2)]
  pre_nor[phenotype == 3, `:=`(ob_nor_fitness_sigma, sigma * d3)]
  pre_nor[phenotype == 4, `:=`(ob_nor_fitness_sigma, sigma * d4)]
  pre_nor[phenotype == 5, `:=`(ob_nor_fitness_sigma, sigma * d5)]
  
  pre_nor[phenotype == 1, `:=`(pre_nor_fitness, predicted_fitness)]
  pre_nor[phenotype == 2, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_fitness, predicted_fitness * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_fitness, predicted_fitness * d5 + e5)]
  
  return(pre_nor)
}

# Wild-type KRAS sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Merge folding free energy and fitness data for AbundancePCA
Folding_pre_ob_fitness <- krasddpcams__merge_ddGf_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  folding_ddG = "path/to/weights_Folding.txt",
  block1_dimsum_df = "path/to/Abundance_block1_fitness.RData",
  block2_dimsum_df = "path/to/Abundance_block2_fitness.RData",
  block3_dimsum_df = "path/to/Abundance_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/Abundance_block3_fitness.RData",
  block5_dimsum_df = "path/to/Abundance_block3_2_fitness.RData",
  wt_aa_input = wt_aa)

# Generate 2D density plots showing global epistasis between folding energy and fitness

# Block1: Relationship between folding free energy and fitness
krasddpcams__plot2d_ddGf_fitness(
  pre_nor = Folding_pre_ob_fitness,
  fold_n = 1,
  mochi_parameters = "path/to/linears_weights_Abundance1.txt",
  phenotypen = 1, 
  RT = 0.001987*(273+30), 
  bin_input = 50)
ggplot2::ggsave("path/to/global_epistasis_block1_ddGf_fitness.pdf", 
                device = cairo_pdf, height = 35, width = 60, units = "mm")

# Block2-1: Relationship between folding free energy and fitness
krasddpcams__plot2d_ddGf_fitness(
  pre_nor = Folding_pre_ob_fitness,
  fold_n = 1,
  mochi_parameters = "path/to/linears_weights_Abundance2_1.txt",
  phenotypen = 2,
  RT = 0.001987*(273+30), 
  bin_input = 50)
ggplot2::ggsave("path/to/global_epistasis_block2_1_ddGf_fitness.pdf", 
                device = cairo_pdf, height = 35, width = 60, units = "mm")

# Block2-2: Relationship between folding free energy and fitness
krasddpcams__plot2d_ddGf_fitness(
  pre_nor = Folding_pre_ob_fitness,
  fold_n = 1,
  mochi_parameters = "path/to/linears_weights_Abundance2_2.txt",
  phenotypen = 3,
  RT = 0.001987*(273+30), 
  bin_input = 50)
ggplot2::ggsave("path/to/global_epistasis_block2_2_ddGf_fitness.pdf", 
                device = cairo_pdf, height = 35, width = 60, units = "mm")

# Block3-1: Relationship between folding free energy and fitness
krasddpcams__plot2d_ddGf_fitness(
  pre_nor = Folding_pre_ob_fitness,
  fold_n = 1,
  mochi_parameters = "path/to/linears_weights_Abundance3_1.txt",
  phenotypen = 4,
  RT = 0.001987*(273+30), 
  bin_input = 50)
ggplot2::ggsave("path/to/global_epistasis_block3_1_ddGf_fitness.pdf", 
                device = cairo_pdf, height = 35, width = 60, units = "mm")

# Block3-2: Relationship between folding free energy and fitness
krasddpcams__plot2d_ddGf_fitness(
  pre_nor = Folding_pre_ob_fitness,
  fold_n = 1,
  mochi_parameters = "path/to/linears_weights_Abundance3_2.txt",
  phenotypen = 5,
  RT = 0.001987*(273+30), 
  bin_input = 50)
ggplot2::ggsave("path/to/global_epistasis_block3_2_ddGf_fitness.pdf", 
                device = cairo_pdf, height = 35, width = 60, units = "mm")
library(krasddpcams)
library(data.table)

### 20250920
## Binding free energy (ddGb) analysis
## for 0901 energy model

# Color scheme for consistent visualization
colour_scheme <- list(
  "blue" = "#1B38A6",        # rgb(27, 56, 166)
  "red" = "#F4270C",         # rgb(244, 39, 12)
  "orange" = "#F4AD0C",      # rgb(244, 173, 12)
  "green" = "#09B636",       # rgb(9, 182, 54)
  "yellow" = "#F1DD10",      # rgb(241, 221, 16)
  "purple" = "#6D17A0",      # rgb(109, 23, 160)
  "pink" = "#FFB0A5",        # rgb(255, 176, 165)
  "light orange" = "#FFE4A5", # rgb(255, 228, 165)
  "light blue" = "#9DACE3",   # rgb(157, 172, 227)
  "light green" = "#97E9AD",  # rgb(151, 233, 173)
  "light red" = "#FF6A56",    # rgb(255, 106, 86)
  "dark red" = "#A31300",     # rgb(163, 19, 0)
  "dark blue" = "#0C226F",    # rgb(12, 34, 111)
  "dark green" = "#007A20"    # rgb(0, 122, 32)
)

## S1i: Non-linear relationships (global epistasis) between observed BindingPCA fitness 
## and both free energies of binding and folding

# Function to merge binding free energy (ddGb) and fitness data for 5 experimental blocks
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

## RAF1 Analysis

# Merge RAF1 binding fitness data
RAF1_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  block1_dimsum_df = "path/to/RAF_block1_fitness.RData",
  block2_dimsum_df = "path/to/RAF_block2_fitness.RData",
  block3_dimsum_df = "path/to/RAF_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/RAF_block3_fitness.RData",
  block5_dimsum_df = "path/to/RAF_block3_2_fitness.RData",
  assay_sele = "RAF1",
  wt_aa_input = wt_aa)

# Generate 3D plots showing global epistasis between folding energy, binding energy, and fitness

# Block1: Relationship between folding energy, binding energy, and fitness
Cairo::CairoPDF(file = "path/to/RAF1_global_epistasis_block1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = RAF1_pre_ob_fitness,
  folding_assay = Abundance1,
  binding_assay = Binding1_RAF1,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding1_RAF1.txt", 
  colour_scheme)
dev.off()

# Block2-1: Relationship between folding energy, binding energy, and fitness
Cairo::CairoPDF(file = "path/to/RAF1_global_epistasis_block2_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = RAF1_pre_ob_fitness,
  folding_assay = Abundance2_1,
  binding_assay = Binding2_1_RAF1,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_1_RAF1.txt", 
  colour_scheme)
dev.off()

# Block2-2: Relationship between folding energy, binding energy, and fitness
Cairo::CairoPDF(file = "path/to/RAF1_global_epistasis_block2_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = RAF1_pre_ob_fitness,
  folding_assay = Abundance2_2,
  binding_assay = Binding2_2_RAF1,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_2_RAF1.txt", 
  colour_scheme)
dev.off()

# Block3-1: Relationship between folding energy, binding energy, and fitness
Cairo::CairoPDF(file = "path/to/RAF1_global_epistasis_block3_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = RAF1_pre_ob_fitness,
  folding_assay = Abundance3_1,
  binding_assay = Binding3_1_RAF1,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_1_RAF1.txt", 
  colour_scheme)
dev.off()

# Block3-2: Relationship between folding energy, binding energy, and fitness
Cairo::CairoPDF(file = "path/to/RAF1_global_epistasis_block3_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = RAF1_pre_ob_fitness,
  folding_assay = Abundance3_2,
  binding_assay = Binding3_2_RAF1,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_2_RAF1.txt", 
  colour_scheme)
dev.off()

## K13 Analysis

# Merge K13 binding fitness data
K13_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  block1_dimsum_df = "path/to/K13_block1_fitness.RData",
  block2_dimsum_df = "path/to/K13_block2_fitness.RData",
  block3_dimsum_df = "path/to/K13_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/K13_block3_fitness.RData",
  block5_dimsum_df = "path/to/K13_block3_2_fitness.RData",
  assay_sele = "K13",
  wt_aa_input = wt_aa)

# Generate K13 3D global epistasis plots
Cairo::CairoPDF(file = "path/to/K13_global_epistasis_block1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K13_pre_ob_fitness,
  folding_assay = Abundance1,
  binding_assay = Binding1_K13,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding1_K13.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K13_global_epistasis_block2_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K13_pre_ob_fitness,
  folding_assay = Abundance2_1,
  binding_assay = Binding2_1_K13,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_1_K13.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K13_global_epistasis_block2_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K13_pre_ob_fitness,
  folding_assay = Abundance2_2,
  binding_assay = Binding2_2_K13,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_2_K13.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K13_global_epistasis_block3_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K13_pre_ob_fitness,
  folding_assay = Abundance3_1,
  binding_assay = Binding3_1_K13,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_1_K13.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K13_global_epistasis_block3_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K13_pre_ob_fitness,
  folding_assay = Abundance3_2,
  binding_assay = Binding3_2_K13,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_2_K13.txt", 
  colour_scheme)
dev.off()

## K19 Analysis

# Merge K19 binding fitness data
K19_pre_ob_fitness <- krasddpcams__merge_ddGb_ob_pre_fitness_5blocks(
  prediction = "path/to/predicted_phenotypes_all.txt",
  block1_dimsum_df = "path/to/K19_block1_fitness.RData",
  block2_dimsum_df = "path/to/K19_block2_fitness.RData",
  block3_dimsum_df = "path/to/K19_block2_2_fitness.RData",
  block4_dimsum_df = "path/to/K19_block3_fitness.RData",
  block5_dimsum_df = "path/to/K19_block3_2_fitness.RData",
  assay_sele = "K19",
  wt_aa_input = wt_aa)

# Generate K19 3D global epistasis plots
Cairo::CairoPDF(file = "path/to/K19_global_epistasis_block1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K19_pre_ob_fitness,
  folding_assay = Abundance1,
  binding_assay = Binding1_K19,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding1_K19.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K19_global_epistasis_block2_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K19_pre_ob_fitness,
  folding_assay = Abundance2_1,
  binding_assay = Binding2_1_K19,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_1_K19.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K19_global_epistasis_block2_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K19_pre_ob_fitness,
  folding_assay = Abundance2_2,
  binding_assay = Binding2_2_K19,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding2_2_K19.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K19_global_epistasis_block3_1_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K19_pre_ob_fitness,
  folding_assay = Abundance3_1,
  binding_assay = Binding3_1_K19,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_1_K19.txt", 
  colour_scheme)
dev.off()

Cairo::CairoPDF(file = "path/to/K19_global_epistasis_block3_2_ddGf_ddGb_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(
  binding_input = K19_pre_ob_fitness,
  folding_assay = Abundance3_2,
  binding_assay = Binding3_2_K19,
  RT = 0.001987*(273+30),
  mochi_parameters = "path/to/linears_weights_Binding3_2_K19.txt", 
  colour_scheme)
dev.off()
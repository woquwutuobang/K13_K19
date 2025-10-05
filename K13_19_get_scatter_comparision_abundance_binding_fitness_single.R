# KRAS DARPins Binding Fitness Analysis
# This script analyzes and visualizes the fitness effects of mutations on KRAS binding to DARPins K13 and K19

library(ggplot2)
library(data.table)
library(dplyr)
library(wlab.block)
library(krasddpcams)

# Wild-type KRAS sequence (residues 2-189)
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Color scheme for consistent visualization
colour_scheme <- list(
  "blue" = "#1B38A6",
  "red" = "#F4270C", 
  "orange" = "#F4AD0C",
  "green" = "#09B636",
  "yellow" = "#F1DD10",
  "purple" = "#6D17A0",
  "pink" = "#FFB0A5",
  "light_orange" = "#FFE4A5",
  "light_blue" = "#9DACE3",
  "light_green" = "#97E9AD",
  "light_red" = "#FF6A56",
  "dark_red" = "#A31300",
  "dark_blue" = "#0C226F",
  "dark_green" = "#007A20"
)

# Function to merge two DiMSum data frames
merge_dimsum_data <- function(merge_1, merge_2) {
  a1 <- as.character(substitute(merge_1))
  a2 <- as.character(substitute(merge_2))
  merge_1[, assay := a1]
  merge_2[, assay := a2]
  output <- rbind(merge_1, merge_2)
  return(output)
}

# Function to identify mutation positions
identify_mutation_positions <- function(input, wt_aa) {
  output <- input
  output[, AA_Pos1 := which(unlist(strsplit(aa_seq, "")) != 
                              unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[1], aa_seq]
  output[, AA_Pos2 := which(unlist(strsplit(aa_seq, "")) != 
                              unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[2], aa_seq]
  
  for (i in 1:188) {
    output[AA_Pos1 == i, mt1 := substr(aa_seq, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos2 == i, mt2 := substr(aa_seq, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos1 == i, wtcodon1 := substr(wt_aa, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos2 == i, wtcodon2 := substr(wt_aa, i, i)]
  }
  
  output[, codon1 := substr(aa_seq, AA_Pos1, AA_Pos1)]
  output[, codon2 := substr(aa_seq, AA_Pos2, AA_Pos2)]
  output[, AA_Pos1 := AA_Pos1 + 1]  # Convert to 1-based indexing
  output[, AA_Pos2 := AA_Pos2 + 1]
  output[, mt1 := paste0(wtcodon1, AA_Pos1, codon1)]
  output[, mt2 := paste0(wtcodon2, AA_Pos2, codon2)]
  
  return(output)
}

# Function to create scatter plot for K13 binding
plot_binding_fitness_k13 <- function(input, assay_sele, anno, colour_scheme) {
  input_abundance <- input[assay == "stab", ]
  input_abundance_single <- krasddpcams__nor_overlap_single_mt_fitness(input_abundance)
  input_binding <- input[assay == assay_sele, ]
  input_binding_single <- krasddpcams__nor_overlap_single_mt_fitness(input_binding)
  
  input_long <- rbind(input_abundance_single, input_binding_single)
  input_dc <- dcast(input_long, nt_seq + aa_seq + Nham_aa + 
                      AA_Pos1 + wtcodon1 ~ assay, 
                    value.var = c("nor_fitness_nooverlap", "nor_fitness_nooverlap_sigma", 
                                  "nor_gr_nooverlap", "nor_gr_nooverlap_sigma"), 
                    drop = TRUE)
  
  input_single_pos <- input_dc
  input_single_pos[, position := AA_Pos1]
  input_single_pos[, WT_AA := wtcodon1]
  
  anno_single <- merge(input_single_pos, anno, by.x = c("position", "WT_AA"), 
                       by.y = c("Pos_real", "codon"), all = TRUE)
  anno_single[, K13_type_bs := "others"]
  anno_single[get(paste0("scHAmin_ligand_", assay_sele)) < 5, 
              K13_type_bs := "binding_interface"]
  
  ggplot() +
    geom_point(data = anno_single[position > 1 & K13_type_bs == "others", ],
               aes(x = nor_fitness_nooverlap_stab, 
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = "black", alpha = 0.6, size = 0.1) +
    geom_point(data = anno_single[position > 1 & K13_type_bs == "binding_interface", ],
               aes(x = nor_fitness_nooverlap_stab, 
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = colour_scheme[["red"]], alpha = 0.6, size = 0.5) +
    theme_classic2() +
    labs(color = NULL) +
    xlab("Protein Abundance Fitness") +
    ylab("K13 Binding Fitness") +
    theme(text = element_text(size = 8),
          legend.position = "right",
          axis.text = element_text(size = 8),
          legend.key.height = unit(3.1, "mm"),
          legend.key.width = unit(3.1, "mm"),
          plot.margin = margin(0, 0, 0, 0)) +
    coord_fixed(ratio = 1.6)
}

# Function to create scatter plot for K19 binding  
plot_binding_fitness_k19 <- function(input, assay_sele, anno, colour_scheme) {
  input_abundance <- input[assay == "stab", ]
  input_abundance_single <- krasddpcams__nor_overlap_single_mt_fitness(input_abundance)
  input_binding <- input[assay == assay_sele, ]
  input_binding_single <- krasddpcams__nor_overlap_single_mt_fitness(input_binding)
  
  input_long <- rbind(input_abundance_single, input_binding_single)
  input_dc <- dcast(input_long, nt_seq + aa_seq + Nham_aa + 
                      AA_Pos1 + wtcodon1 ~ assay, 
                    value.var = c("nor_fitness_nooverlap", "nor_fitness_nooverlap_sigma", 
                                  "nor_gr_nooverlap", "nor_gr_nooverlap_sigma"), 
                    drop = TRUE)
  
  input_single_pos <- input_dc
  input_single_pos[, position := AA_Pos1]
  input_single_pos[, WT_AA := wtcodon1]
  
  anno_single <- merge(input_single_pos, anno, by.x = c("position", "WT_AA"), 
                       by.y = c("Pos_real", "codon"), all = TRUE)
  anno_single[, K19_type_bs := "others"]
  anno_single[get(paste0("scHAmin_ligand_", assay_sele)) < 5, 
              K19_type_bs := "binding_interface"]
  
  ggplot() +
    geom_point(data = anno_single[position > 1 & K19_type_bs == "others", ],
               aes(x = nor_fitness_nooverlap_stab, 
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = "black", alpha = 0.6, size = 0.1) +
    geom_point(data = anno_single[position > 1 & K19_type_bs == "binding_interface", ],
               aes(x = nor_fitness_nooverlap_stab, 
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = colour_scheme[["red"]], alpha = 0.6, size = 0.5) +
    theme_classic2() +
    labs(color = NULL) +
    xlab("Protein Abundance Fitness") +
    ylab("K19 Binding Fitness") +
    theme(text = element_text(size = 8),
          legend.position = "right", 
          axis.text = element_text(size = 8),
          legend.key.height = unit(3.1, "mm"),
          legend.key.width = unit(3.1, "mm"),
          plot.margin = margin(0, 0, 0, 0)) +
    coord_fixed(ratio = 1.5)
}

# Load and normalize experimental data
stability_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/abundance_block1.RData",
  block2_dimsum_df = "path/to/abundance_block2.RData", 
  block3_dimsum_df = "path/to/abundance_block3.RData"
)

K13_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/K13_block1.RData",
  block2_dimsum_df = "path/to/K13_block2.RData",
  block3_dimsum_df = "path/to/K13_block3.RData"
)

K19_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/K19_block1.RData",
  block2_dimsum_df = "path/to/K19_block2.RData",
  block3_dimsum_df = "path/to/K19_block3.RData"
)

# Simplify variable names
stab <- stability_nor_df
K13 <- K13_nor_df  
K19 <- K19_nor_df

# Merge abundance and binding data
all_data_K13 <- merge_dimsum_data(stab, K13)
all_data_K19 <- merge_dimsum_data(stab, K19)

# Identify mutation positions
all_data_pos_K13 <- identify_mutation_positions(all_data_K13, wt_aa)
all_data_pos_K19 <- identify_mutation_positions(all_data_K19, wt_aa)

# Load annotation data
anno <- fread("path/to/annotation_data.csv")
anno[, Pos_real := Pos]

# Generate binding fitness plots
plot_K13 <- plot_binding_fitness_k13(
  input = all_data_pos_K13,
  assay_sele = "K13",
  anno = anno,
  colour_scheme = colour_scheme
)

plot_K19 <- plot_binding_fitness_k19(
  input = all_data_pos_K19, 
  assay_sele = "K19",
  anno = anno,
  colour_scheme = colour_scheme
)

# Display plots
print(plot_K13)
print(plot_K19)

# Save plots
ggsave("K13_binding_fitness.pdf", plot = plot_K13, device = cairo_pdf, 
       height = 4, width = 4)
ggsave("K19_binding_fitness.pdf", plot = plot_K19, device = cairo_pdf,
       height = 4, width = 4)
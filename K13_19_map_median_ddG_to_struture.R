# KRAS Median ΔΔG Structure Mapping
# This script calculates median ΔΔG values and maps them to protein structures

library(data.table)
library(krasddpcams)
library(dplyr)
library(bio3d)

# Function 1: Calculate median ΔΔG per position
krasddpcams__get_median_ddG <- function(ddG, assay_sele) {
  # Read ddG data
  ddG <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
  
  # Calculate median for each Pos_real
  median_ddG <- ddG[Pos_real > 1,
                    .(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE)),
                    by = "Pos_real"]
  
  return(median_ddG)
}



# Function 2: Map median ΔΔG to structure B-factors
krasddpcams__median_ddG_structure <- function(
    median_ddG_table,
    input_PDB,
    chain_KRAS = "A",
    Pos_correction = 0,
    output_PDB_file
) {
  library(bio3d)
  
  # Read input structure
  input_structure <- read.pdb(input_PDB)
  output_PDB <- input_structure
  
  # Set all B-factors to 0 for specified chain
  res_indices <- output_PDB$atom$resno[output_PDB$atom$chain == chain_KRAS]
  for (i in unique(res_indices)) {
    output_PDB$atom$b[output_PDB$atom$resno == i & 
                        output_PDB$atom$chain == chain_KRAS] <- 0
  }
  
  # Iterate through each residue position in median_ddG table and write to B-factors
  for (i in median_ddG_table$Pos_real) {
    val <- median_ddG_table$median_ddG[median_ddG_table$Pos_real == i]
    # If multiple values exist, take unique value (avoid length>1 error)
    if (length(val) > 1) val <- unique(val)
    if (length(val) == 1 && !is.na(val)) {
      output_PDB$atom$b[output_PDB$atom$resno == (i + Pos_correction) & 
                          output_PDB$atom$chain == chain_KRAS] <- val
    }
  }
  
  # Write new structure file
  write.pdb(output_PDB, file = output_PDB_file)
}



# Function 3: Special function for RALGDS
krasddpcams__median_ddG_structure_for_RALGDS <- function(
    median_ddG_table,
    input_PDB,
    chain_KRAS = "A",
    Pos_correction = 0,
    output_PDB_file
) {
  library(bio3d)
  
  # Read input structure
  input_structure <- read.pdb(input_PDB)
  output_PDB <- input_structure
  
  # Set all B-factors to 0 for specified chain
  res_indices <- output_PDB$atom$resno[output_PDB$atom$chain == chain_KRAS]
  for (i in unique(res_indices)) {
    output_PDB$atom$b[output_PDB$atom$resno == i & 
                        output_PDB$atom$chain == chain_KRAS] <- 0
  }
  
  # Iterate through each residue position in median_ddG table and write to B-factors
  for (i in median_ddG_table$Pos_real) {
    val <- median_ddG_table$median_ddG[median_ddG_table$Pos_real == i]
    # If multiple values exist, take unique value (avoid length>1 error)
    if (length(val) > 1) val <- unique(val)
    if (length(val) == 1 && !is.na(val)) {
      # Add 200 to Pos_real
      residue_number <- i + 200 + Pos_correction
      output_PDB$atom$b[output_PDB$atom$resno == residue_number & 
                          output_PDB$atom$chain == chain_KRAS] <- val
    }
  }
  
  # Write new structure file
  write.pdb(output_PDB, file = output_PDB_file)
}

# =====================
# Process binding partners
# =====================

## K13
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/6h46.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "path/to/K13_median_ddG.pdb"
)

## K19
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_K19.txt",
  assay_sele = "K19"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/6h47.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "path/to/K19_median_ddG.pdb"
)

## RAF1
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_RAF1.txt",
  assay_sele = "RAF1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/6vjj.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "path/to/RAF1_median_ddG.pdb"
)

## SOS1-Q
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_SOS.txt",
  assay_sele = "SOS1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/1nvw.pdb",
  chain_KRAS = "Q",
  Pos_correction = 0,
  output_PDB_file = "path/to/SOS1_Q_median_ddG.pdb"
)

## SOS1-R
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_SOS.txt",
  assay_sele = "SOS1"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/1nvw.pdb",
  chain_KRAS = "R",
  Pos_correction = 0,
  output_PDB_file = "path/to/SOS1_R_median_ddG.pdb"
)

## K55
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_K55.txt",
  assay_sele = "K55"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/5mla.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "path/to/K55_median_ddG.pdb"
)

## K27
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_K27.txt",
  assay_sele = "K27"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/5mlb.pdb",
  chain_KRAS = "A",
  Pos_correction = 0,
  output_PDB_file = "path/to/K27_median_ddG.pdb"
)

## RALGDS
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_RAL.txt",
  assay_sele = "RALGDS"
)

krasddpcams__median_ddG_structure_for_RALGDS(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/1lfd.pdb",
  chain_KRAS = "B",
  Pos_correction = 0,
  output_PDB_file = "path/to/RALGDS_median_ddG.pdb"
)

## PIK3CG
median_by_pos <- krasddpcams__get_median_ddG(
  ddG = "path/to/weights_Binding_PI3.txt",
  assay_sele = "PIK3CG"
)

krasddpcams__median_ddG_structure(
  median_ddG_table = median_by_pos,
  input_PDB = "path/to/1he8.pdb",
  chain_KRAS = "B",
  Pos_correction = 0,
  output_PDB_file = "path/to/PIK3CG_median_ddG.pdb"
)
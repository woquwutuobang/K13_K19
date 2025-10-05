# KRAS Binding Hotspot Analysis
# This script identifies hotspot residues in KRAS binding interfaces based on ddG values

library(data.table)
library(krasddpcams)
library(dplyr)

# Function to calculate median ddG values per position
calculate_median_ddG <- function(ddG, assay_sele) {
  # Read ddG data
  ddG_data <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
  
  # Calculate median ddG for each position
  median_ddG <- ddG_data[Pos_real > 1,
                         .(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE)),
                         by = "Pos_real"]
  
  return(median_ddG)
}

# Function to identify hotspot residues
identify_hotspots <- function(median_data, interface_positions, sd_multiplier = 1.5) {
  # Calculate mean and standard deviation
  mean_ddG <- mean(median_data$median_ddG, na.rm = TRUE)
  sd_ddG <- sd(median_data$median_ddG, na.rm = TRUE)
  
  # Set threshold (mean + n * standard deviation)
  threshold <- mean_ddG + sd_multiplier * sd_ddG
  
  # Identify hotspots (interface residues with ddG above threshold)
  median_data <- median_data %>%
    mutate(hotspot = ifelse(Pos_real %in% interface_positions & median_ddG >= threshold, 
                            TRUE, FALSE))
  
  return(median_data)
}

# Define binding interface residues for DARPins
interface_residues <- c(88, 91, 87, 129, 90, 133, 94, 137, 95, 68, 136, 99, 102, 101, 107, 98)

# Analyze K13 binding hotspots
cat("Analyzing K13 binding hotspots...\n")
k13_median <- calculate_median_ddG(
  ddG = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13"
)

k13_hotspots <- identify_hotspots(k13_median, interface_residues, sd_multiplier = 1.5)

# Display K13 interface residues analysis
cat("K13 Interface Residues Analysis:\n")
k13_interface_analysis <- k13_hotspots %>%
  filter(Pos_real %in% interface_residues) %>%
  select(Pos_real, median_ddG, hotspot)

print(k13_interface_analysis)

# Extract K13 hotspot positions
k13_hotspot_positions <- k13_hotspots %>%
  filter(hotspot == TRUE) %>%
  pull(Pos_real)

cat("\nK13 Hotspot Residues:", paste(sort(k13_hotspot_positions), collapse = ", "), "\n\n")

# Analyze K19 binding hotspots
cat("Analyzing K19 binding hotspots...\n")
k19_median <- calculate_median_ddG(
  ddG = "path/to/weights_Binding_K19.txt", 
  assay_sele = "K19"
)

k19_hotspots <- identify_hotspots(k19_median, interface_residues, sd_multiplier = 1.5)

# Display K19 interface residues analysis
cat("K19 Interface Residues Analysis:\n")
k19_interface_analysis <- k19_hotspots %>%
  filter(Pos_real %in% interface_residues) %>%
  select(Pos_real, median_ddG, hotspot)

print(k19_interface_analysis)

# Extract K19 hotspot positions
k19_hotspot_positions <- k19_hotspots %>%
  filter(hotspot == TRUE) %>%
  pull(Pos_real)

cat("\nK19 Hotspot Residues:", paste(sort(k19_hotspot_positions), collapse = ", "), "\n")

# Summary report
cat("\n=== SUMMARY ===\n")
cat("K13 Hotspots:", length(k13_hotspot_positions), "residues\n")
cat("K19 Hotspots:", length(k19_hotspot_positions), "residues\n")
cat("Shared Interface Residues:", length(interface_residues), "residues\n")




## plot step

#########20251002
## K13
## open "C:\\Users\\36146\\OneDrive - USTC\\DryLab\\Data_analysis_scripts\\median ddG endow different type structure\\results\\20250904_median for structure\\abs_median\\6h46_K13_abs_median_ddG_use_0901_data.pdb"
## set bgColor white;show #1/A surface;color bfactor #1/A range -2,2;color #1/B grey;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;hide solvent;transparency 50 cartoon;
## save "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20251002_BI_structure/K13.cxs"
## hide #1/B cartoon;hide #1/B atoms;sel #1/A:91 ,94 ,95 ,101;label sel residue size 13;label sel residue color #0C226F;zoom
## save "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20251002_BI_structure/K13_label_BI hotpots.cxs"



## K19
## open "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20250904_median for structure/abs_median/6h47_K19_abs_median_ddG_use_0901_data.pdb"
## matchmaker #!2 to #1;show #2/A surface;color bfactor #2/A range -2,2;color #2/B grey;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;hide solvent;transparency 50 cartoon;
## save "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20251002_BI_structure/K19.cxs"
## hide #2/B cartoon;hide #2/B atoms;sel #2/A:91 ,94 ,95,98;label sel residue size 13;label sel residue color #0C226F;zoom
## save "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20251002_BI_structure/K19_label_BI hotpots_2.cxs"
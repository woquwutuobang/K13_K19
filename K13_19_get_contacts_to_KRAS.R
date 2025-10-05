library(bio3d)
library(wlab.block)
library(dplyr)
library(tidyr)
library(ggplot2)

# DARPin-KRAS Binding Interface Analysis
# This script analyzes and visualizes the binding interface residues between DARPin proteins (K13, K19) and KRAS

# Function to calculate minimum interchain distances and identify binding interface residues
calculate_binding_interface <- function(pdb_file, query_chain, target_chain, distance_cutoff = 5) {
  # Calculate minimum distances between chains
  dist_data <- bio3D_minimum_interchain_distances(
    input_file = pdb_file,
    chain_query = query_chain,
    chain_target = target_chain
  )
  
  # Rename columns with specific interaction identifiers
  interaction_name <- paste0(query_chain, "_to_", target_chain)
  colnames(dist_data)[2:3] <- paste0(colnames(dist_data)[2:3], "_", interaction_name)
  colnames(dist_data)[1] <- "Pos_real"
  
  # Classify binding residues based on distance cutoff (5Å)
  dist_data[, bind_status := 'no']
  dist_data[get(paste0("scHAmin_ligand_", interaction_name)) < distance_cutoff, bind_status := 'yes']
  dist_data[bind_status == "yes", site_type := "Binding interface site"]
  
  return(dist_data)
}

### K13 to KRAS binding interface analysis
K13_interface <- calculate_binding_interface(
  input_file = "path/to/6h46.pdb",  # K13-KRAS complex structure
  chain_query = "B",                # DARPin K13 chain
  chain_target = "A"                # KRAS chain
)

# Extract binding interface positions for K13
K13_binding_positions <- K13_interface %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)

print("K13 binding interface positions:")
print(K13_binding_positions)

### K19 to KRAS binding interface analysis  
K19_interface <- calculate_binding_interface(
  input_file = "path/to/6h47.pdb",  # K19-KRAS complex structure
  chain_query = "B",                # DARPin K19 chain
  chain_target = "A"                # KRAS chain
)

# Extract binding interface positions for K19
K19_binding_positions <- K19_interface %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)

print("K19 binding interface positions:")
print(K19_binding_positions)

# Define DARPin binding interface residues from structural analysis
# These residues are critical for DARPin-KRAS interactions based on crystal structures

K13_residues <- c("ASP132", "MET101", "PHE134", "ASP99", "THR100", "HIS68", "SER112",
                  "MET78", "TRP35", "TRP37", "GLN70", "LEU75", "TRP45", "TRP46", "LEU42", "ARG12")

K19_residues <- c("ARG111", "HIS107", "MET101", "HIS68", "THR100", "ASP99", "HIS114",
                  "SER112", "MET78", "TRP35", "GLN70", "TRP37", "ARG12", "LEU42", "LEU75", "TRP46", "TRP45")

# Combine all unique residues from both DARPin interfaces
all_residues <- unique(c(K13_residues, K19_residues))

# Create data frame to track residue presence in each DARPin interface
interface_df <- data.frame(
  residue = all_residues,
  K13 = as.integer(all_residues %in% K13_residues),
  K19 = as.integer(all_residues %in% K19_residues)
)

# Convert to long format for visualization
interface_long <- interface_df %>%
  pivot_longer(cols = c("K13", "K19"), 
               names_to = "DARPin", 
               values_to = "present")

# Classify contact types for visualization
interface_long <- interface_long %>%
  mutate(contact_type = case_when(
    # Residues common to both K13 and K19 interfaces
    residue %in% K13_residues & residue %in% K19_residues ~ "Common contacts",
    # Residues specific to individual DARPins
    present == 1 ~ "Specific contacts",
    # Residues not present in this DARPin interface
    present == 0 ~ "None"
  ))

# Improve DARPin labels for the plot
interface_long <- interface_long %>%
  mutate(DARPin = recode(DARPin,
                         "K13" = "DARPin K13 to KRAS Binding Interface Residues",
                         "K19" = "DARPin K19 to KRAS Binding Interface Residues"))

# Create contact map visualization
contact_plot <- ggplot(interface_long, aes(x = residue, y = DARPin, fill = contact_type)) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(
      "Common contacts" = "#F4270C",      # Red for shared residues
      "Specific contacts" = "#FFB0A5",    # Light red for specific residues
      "None" = "white"                     # White for absent residues
    ),
    name = "Contact Type"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10)
  ) +
  labs(title = "DARPin-KRAS Binding Interface Residue Comparison")

# Display the plot
print(contact_plot)

# Save the contact map
ggsave("path/to/DARPin_KRAS_binding_interface_comparison.pdf", 
       plot = contact_plot,
       width = 8, 
       height = 3, 
       dpi = 300)

# Summary statistics for binding interface analysis
cat("\nDARPin-KRAS Binding Interface Summary:\n")
cat("========================================\n")
cat("Total unique residues across both DARPins:", length(all_residues), "\n")
cat("K13 interface residues:", sum(interface_df$K13), "\n")
cat("K19 interface residues:", sum(interface_df$K19), "\n")

# Calculate common residues between K13 and K19
common_residues <- sum(interface_df$K13 == 1 & interface_df$K19 == 1)
cat("Common residues between K13 and K19:", common_residues, "\n")

# Calculate residue overlap percentages
cat("K13-K19 overlap:", round(common_residues / sum(interface_df$K13) * 100, 1), "% of K13 residues\n")
cat("K19-K13 overlap:", round(common_residues / sum(interface_df$K19) * 100, 1), "% of K19 residues\n")

# Additional structural insights
cat("\nStructural Insights:\n")
cat("• Binding interface defined by residues within 5Å distance cutoff\n")
cat("• Common residues indicate conserved binding motifs\n") 
cat("• Specific residues determine DARPin binding specificity\n")
cat("• Tryptophan residues (TRP) often involved in hydrophobic interactions\n")
cat("• Charged residues (ASP, ARG, HIS) mediate electrostatic interactions\n")
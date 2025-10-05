# KRAS Binding Interface Comparison
# This script creates a contact map comparing binding interfaces of different KRAS partners

library(ggplot2)
library(reshape2)
library(dplyr)

# Wild-type KRAS sequence
wt_aa <- unlist(strsplit("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM", ""))
positions <- 1:length(wt_aa)

# Define binding interface residues for each partner
K13_binding_interface <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138) 
K19_binding_interface <- c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137) 
RAF1_binding_interface <- c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71) 
K27_binding_interface <- c(21, 24, 25, 27, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)

# Create data frame for plotting
df <- data.frame(
  Position = rep(positions, 4),
  Residue = rep(wt_aa, 4),
  Partner = rep(c("DARPin K13", "DARPin K19", "RAF1-RBD", "K27"), each = length(wt_aa)),
  Contact = FALSE
)

# Mark contact residues for each partner
df$Contact[df$Partner == "DARPin K13" & df$Position %in% K13_binding_interface] <- TRUE
df$Contact[df$Partner == "DARPin K19" & df$Position %in% K19_binding_interface] <- TRUE
df$Contact[df$Partner == "RAF1-RBD" & df$Position %in% RAF1_binding_interface] <- TRUE
df$Contact[df$Partner == "K27" & df$Position %in% K27_binding_interface] <- TRUE

# Identify common contact residues (appear in 2 or more partners)
all_contacts <- list(K13_binding_interface, K19_binding_interface, RAF1_binding_interface, K27_binding_interface)
contact_freq <- table(unlist(all_contacts))
common_contact_positions <- as.integer(names(contact_freq[contact_freq >= 2]))

# Mark common contacts
df$Common <- df$Position %in% common_contact_positions & df$Contact

# Get all unique contact positions for plotting
all_contact_positions <- unique(unlist(all_contacts))

# Prepare data for plotting
df_plot <- df %>%
  filter(Position %in% all_contact_positions) %>%
  mutate(ResLabel = paste0(Residue, Position))  

# Create residue label mapping
residue_labels <- df_plot %>%
  distinct(Position, Residue) %>%
  arrange(Position) %>%
  mutate(ResLabel = paste0(Residue, Position)) %>%
  pull(ResLabel)

# Set factor levels for residue labels
df_plot$ResLabel <- paste0(df_plot$Residue, df_plot$Position)
df_plot$ResLabel <- factor(df_plot$ResLabel, levels = residue_labels)

# Define contact types for coloring
df_plot$ContactType <- ifelse(df_plot$Common, "Common Contact", 
                              ifelse(df_plot$Contact, "Unique Contact", "No Contact"))
df_plot$ContactType <- factor(df_plot$ContactType, 
                              levels = c("No Contact", "Unique Contact", "Common Contact"))

# Set partner order (top to bottom in plot)
df_plot$Partner <- factor(df_plot$Partner, 
                          levels = c("RAF1-RBD", "K27", "DARPin K13", "DARPin K19"))

# Create contact map plot
contact_plot <- ggplot(df_plot, aes(x = ResLabel, y = Partner, fill = ContactType)) +
  geom_tile(color = "black", width = 1, height = 1) +
  scale_fill_manual(
    values = c(
      "No Contact" = "white",
      "Unique Contact" = "yellow", 
      "Common Contact" = "orange"
    ),
    name = "Contact Type"
  ) +
  scale_y_discrete(limits = rev(levels(df_plot$Partner))) +  # Reverse Y-axis order
  coord_fixed(ratio = 1) +  # Maintain square aspect ratio
  theme_minimal(base_size = 14) +
  labs(
    x = "Residue Position", 
    y = "Binding Partner", 
    title = "KRAS Binding Interface Comparison",
    subtitle = "Comparison of contact residues across different binding partners"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Display plot
print(contact_plot)

# Save plot
ggsave("KRAS_binding_interface_comparison.pdf",
       plot = contact_plot,
       width = 10, 
       height = 5,
       dpi = 300)


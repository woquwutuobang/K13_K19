library(ggplot2)
library(dplyr)

# KRAS Binding Interface Contact Map Analysis
# This script creates heatmaps showing common and unique structural contacts 
# between KRAS and different binding partners

# Wild-type KRAS sequence
wt_aa <- unlist(strsplit("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM", ""))
positions <- 1:length(wt_aa)

### Analysis 1: 7 binding partners + 2 SOS interfaces

# Define all binding interfaces from structural data
K13_binding_interface <- c(63,68,87,88,90,91,92,94,95,96,97,98,99,101,102,105,106,107,129,133,136,137,138) 
K19_binding_interface <- c(68,87,88,90,91,92,94,95,97,98,99,101,102,105,107,108,125,129,133,136,137) 
RAF1_binding_interface <- c(21,25,29,31,33,36,37,38,39,40,41,67,71) 
SOS1_binding_interface_1 <- c(1, 22, 24, 25, 26, 27, 31, 33, 34, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 61, 64, 65, 66, 67, 70, 149, 153)
SOS1_binding_interface_2 <- c(15, 17, 18, 21, 25, 30, 31, 32, 33, 34, 35, 37, 39, 40, 41, 54, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 73, 95, 102, 103, 105)
K55_binding_interface <- c(5,24,25,31,33,36,37,38,39,40,54,56,64,66,67,70,73,74)
K27_binding_interface <- c(21,24,25,27,29,30,31,32,33,34,35,36,38,39,40,41,43,52,54,67,70,71)
RALGDS_binding_interface <- c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67)
PIK3CG_binding_interface <- c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73)

# Create data frame for all binding partners
df <- data.frame(
  Position = rep(positions, 9),
  Residue = rep(wt_aa, 9),
  Partner = rep(c("DARPin K13", "DARPin K19", "RAF1-RBD", "SOS1-Interface1", 
                  "SOS1-Interface2", "DARPin K55", "DARPin K27", 
                  "RALGDS-RBD", "PIK3CG-RBD"), each = length(wt_aa)),
  Contact = FALSE
)

# Set contact positions for each binding partner
df$Contact[df$Partner == "DARPin K13" & df$Position %in% K13_binding_interface] <- TRUE
df$Contact[df$Partner == "DARPin K19" & df$Position %in% K19_binding_interface] <- TRUE
df$Contact[df$Partner == "RAF1-RBD" & df$Position %in% RAF1_binding_interface] <- TRUE
df$Contact[df$Partner == "SOS1-Interface1" & df$Position %in% SOS1_binding_interface_1] <- TRUE
df$Contact[df$Partner == "SOS1-Interface2" & df$Position %in% SOS1_binding_interface_2] <- TRUE
df$Contact[df$Partner == "DARPin K55" & df$Position %in% K55_binding_interface] <- TRUE
df$Contact[df$Partner == "DARPin K27" & df$Position %in% K27_binding_interface] <- TRUE
df$Contact[df$Partner == "RALGDS-RBD" & df$Position %in% RALGDS_binding_interface] <- TRUE
df$Contact[df$Partner == "PIK3CG-RBD" & df$Position %in% PIK3CG_binding_interface] <- TRUE

# Identify common contact positions (appear in 2 or more partners)
all_contacts <- list(K13_binding_interface, K19_binding_interface, 
                     RAF1_binding_interface, SOS1_binding_interface_1, 
                     SOS1_binding_interface_2, K55_binding_interface,
                     K27_binding_interface, RALGDS_binding_interface,
                     PIK3CG_binding_interface)

contact_freq <- table(unlist(all_contacts))
common_contact_positions <- as.integer(names(contact_freq[contact_freq >= 2]))

df$Common <- df$Position %in% common_contact_positions & df$Contact

# Extract all unique contact positions for plotting
all_contact_positions <- sort(unique(c(K13_binding_interface, K19_binding_interface, 
                                       RAF1_binding_interface, SOS1_binding_interface_1,
                                       SOS1_binding_interface_2, K55_binding_interface,
                                       K27_binding_interface, RALGDS_binding_interface,
                                       PIK3CG_binding_interface)))

# Create plotting subset
df_plot <- df %>%
  filter(Position %in% all_contact_positions) %>%
  mutate(ResLabel = paste0(Residue, Position))

# Create unique residue label mapping
residue_labels <- df_plot %>%
  distinct(Position, Residue) %>%
  arrange(Position) %>%
  mutate(ResLabel = paste0(Residue, Position)) %>%
  pull(ResLabel)

# Set factor levels
df_plot$ResLabel <- factor(df_plot$ResLabel, levels = residue_labels)

# Set partner order as specified
df_plot$Partner <- factor(df_plot$Partner, 
                          levels = c("RAF1-RBD", "SOS1-Interface1", "SOS1-Interface2",
                                     "DARPin K55", "DARPin K27", "RALGDS-RBD", 
                                     "PIK3CG-RBD", "DARPin K13", "DARPin K19"))

# Create fill categories
df_plot$Fill <- ifelse(df_plot$Common, "Common contact", 
                       ifelse(df_plot$Contact, "Contact", "None"))
df_plot$Fill <- factor(df_plot$Fill, levels = c("None", "Contact", "Common contact"))

# Create contact map heatmap
contact_heatmap_9partners <- ggplot(df_plot, aes(x = ResLabel, y = Partner, fill = Fill)) +
  geom_tile(color = "black", width = 1, height = 1) +
  scale_fill_manual(
    values = c(
      "None" = "white",
      "Contact" = "yellow",
      "Common contact" = "orange"
    ),
    name = "Contact Type"
  ) +
  scale_y_discrete(limits = rev(levels(df_plot$Partner))) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "KRAS Residue", 
    y = "Binding Partner", 
    title = "KRAS Binding Interface Contact Map (9 Interfaces)",
    subtitle = "Common contacts (orange) appear in ≥2 binding partners"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the plot
ggsave("path/to/kras_binding_interface_contact_map_9partners.pdf", 
       plot = contact_heatmap_9partners, 
       width = 14, height = 16, dpi = 300)

### Analysis 2: 7 binding partners + 1 SOS interface

# Define binding interfaces (single SOS interface)
SOS1_binding_interface <- c(1, 22, 24, 25, 26, 27, 31, 33, 34, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 61, 64, 65, 66, 67, 70, 149, 153)

# Create data frame for 8 binding partners
df2 <- data.frame(
  Position = rep(positions, 8),
  Residue = rep(wt_aa, 8),
  Partner = rep(c("DARPin K13", "DARPin K19", "RAF1-RBD", "SOS1-Interface", 
                  "DARPin K55", "DARPin K27", "RALGDS-RBD", "PIK3CG-RBD"), 
                each = length(wt_aa)),
  Contact = FALSE
)

# Set contact positions
df2$Contact[df2$Partner == "DARPin K13" & df2$Position %in% K13_binding_interface] <- TRUE
df2$Contact[df2$Partner == "DARPin K19" & df2$Position %in% K19_binding_interface] <- TRUE
df2$Contact[df2$Partner == "RAF1-RBD" & df2$Position %in% RAF1_binding_interface] <- TRUE
df2$Contact[df2$Partner == "SOS1-Interface" & df2$Position %in% SOS1_binding_interface] <- TRUE
df2$Contact[df2$Partner == "DARPin K55" & df2$Position %in% K55_binding_interface] <- TRUE
df2$Contact[df2$Partner == "DARPin K27" & df2$Position %in% K27_binding_interface] <- TRUE
df2$Contact[df2$Partner == "RALGDS-RBD" & df2$Position %in% RALGDS_binding_interface] <- TRUE
df2$Contact[df2$Partner == "PIK3CG-RBD" & df2$Position %in% PIK3CG_binding_interface] <- TRUE

# Identify common contact positions
all_contacts2 <- list(K13_binding_interface, K19_binding_interface, 
                      RAF1_binding_interface, SOS1_binding_interface, 
                      K55_binding_interface, K27_binding_interface, 
                      RALGDS_binding_interface, PIK3CG_binding_interface)

contact_freq2 <- table(unlist(all_contacts2))
common_contact_positions2 <- as.integer(names(contact_freq2[contact_freq2 >= 2]))

df2$Common <- df2$Position %in% common_contact_positions2 & df2$Contact

# Extract all unique contact positions
all_contact_positions2 <- sort(unique(c(K13_binding_interface, K19_binding_interface, 
                                        RAF1_binding_interface, SOS1_binding_interface,
                                        K55_binding_interface, K27_binding_interface, 
                                        RALGDS_binding_interface, PIK3CG_binding_interface)))

# Create plotting subset
df_plot2 <- df2 %>%
  filter(Position %in% all_contact_positions2) %>%
  mutate(ResLabel = paste0(Residue, Position))

# Set factor levels
df_plot2$ResLabel <- factor(df_plot2$ResLabel, levels = residue_labels)
df_plot2$Partner <- factor(df_plot2$Partner, 
                           levels = c("RAF1-RBD", "SOS1-Interface",
                                      "DARPin K55", "DARPin K27", "RALGDS-RBD", 
                                      "PIK3CG-RBD", "DARPin K13", "DARPin K19"))

# Create fill categories
df_plot2$Fill <- ifelse(df_plot2$Common, "Common contact", 
                        ifelse(df_plot2$Contact, "Contact", "None"))
df_plot2$Fill <- factor(df_plot2$Fill, levels = c("None", "Contact", "Common contact"))

# Create second contact map heatmap
contact_heatmap_8partners <- ggplot(df_plot2, aes(x = ResLabel, y = Partner, fill = Fill)) +
  geom_tile(color = "black", width = 1, height = 1) +
  scale_fill_manual(
    values = c(
      "None" = "white",
      "Contact" = "yellow", 
      "Common contact" = "orange"
    ),
    name = "Contact Type"
  ) +
  scale_y_discrete(limits = rev(levels(df_plot2$Partner))) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "KRAS Residue", 
    y = "Binding Partner", 
    title = "KRAS Binding Interface Contact Map (8 Interfaces)",
    subtitle = "Common contacts (orange) appear in ≥2 binding partners"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the second plot
ggsave("path/to/kras_binding_interface_contact_map_8partners.pdf", 
       plot = contact_heatmap_8partners, 
       width = 14, height = 16, dpi = 300)

# Print summary statistics
cat("KRAS Binding Interface Analysis Summary:\n")
cat("========================================\n")
cat("Analysis 1 (9 interfaces):\n")
cat("- Total unique contact positions:", length(all_contact_positions), "\n")
cat("- Common contact positions (≥2 partners):", length(common_contact_positions), "\n")
cat("- Binding partners: RAF1, SOS1 (2 interfaces), K55, K27, RALGDS, PIK3CG, K13, K19\n\n")

cat("Analysis 2 (8 interfaces):\n") 
cat("- Total unique contact positions:", length(all_contact_positions2), "\n")
cat("- Common contact positions (≥2 partners):", length(common_contact_positions2), "\n")
cat("- Binding partners: RAF1, SOS1 (1 interface), K55, K27, RALGDS, PIK3CG, K13, K19\n")
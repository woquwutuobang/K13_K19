# KRAS sequence annotation plot
# This script visualizes various functional and structural features of the KRAS protein

# 1. KRAS reference sequence
KRAS_sequence <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 2. Binding interface sites
K13_binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 
                                101, 102, 105, 106, 107, 129, 133, 136, 137, 138)

RAF1_binding_interface_site <- c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71) 

# 3. GTP binding pocket residues
GTP_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# 4. Functional loops
functional_loop <- data.frame(
  xstart = c(10, 25, 58),
  xend = c(17, 40, 76),
  col = c("P-loop", "switch I", "switch II")
)

# 5. Common core residues
common_core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40,
                          42, 44, 46, 51, 52, 53, 54, 55, 56, 57, 58, 68, 71, 72, 75, 77, 78, 79,
                          80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110,
                          111, 112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141,
                          142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 159, 160, 162, 163)

# 6. Secondary structure elements
rects_sheet <- data.frame(
  xstart = c(3, 38, 51, 77, 109, 139),
  xend = c(9, 44, 57, 84, 115, 143),
  col = c("β1", "β2", "β3", "β4", "β5", "β6")
)

rects_helix <- data.frame(
  xstart = c(15, 67, 87, 127, 148),
  xend = c(24, 73, 104, 136, 166),
  col = c("α1", "α2", "α3", "α4", "α5")
)

# Create data frame for sequence plotting
seq_length <- nchar(KRAS_sequence)
seq_df <- data.frame(
  position = 1:seq_length,
  residue = strsplit(KRAS_sequence, "")[[1]],
  K13_binding = ifelse(1:seq_length %in% K13_binding_interface_site, "#1B38A6", "black"),
  RAF1_binding = ifelse(1:seq_length %in% RAF1_binding_interface_site, "#F4270C", "black")
)

# Generate the plot
library(ggplot2)

ggplot() +
  # Sequence residues (y = 0)
  geom_text(data = seq_df, aes(x = position, y = 0, label = residue, 
                               color = interaction(K13_binding != "black", RAF1_binding != "black")),
            size = 3, family = "mono") +
  scale_color_manual(values = c("black", "#B6F500", "#FF3F33", "#FF2DF1"),
                     labels = c("Other", "K13 binding", "RAF1 binding", "Both"),
                     name = "Binding Sites") +
  
  # GTP binding pocket annotation (y = -0.3 to -0.2)
  geom_rect(data = data.frame(position = GTP_pocket),
            aes(xmin = position - 0.5, xmax = position + 0.5, ymin = -0.3, ymax = -0.2),
            fill = "#1B38A6", alpha = 0.5) +
  annotate("text", x = mean(range(GTP_pocket)), y = -0.25, label = "GTP Pocket", 
           color = "#1B38A6", size = 3, fontface = "bold") +
  
  # Functional loops annotation (y = -0.5 to -0.4)
  geom_rect(data = functional_loop, 
            aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -0.5, ymax = -0.4, fill = col),
            alpha = 0.3) +
  geom_text(data = functional_loop, 
            aes(x = (xstart + xend)/2, y = -0.45, label = col),
            size = 3, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("#9DACE3", "#DE3163", "#F1DD10"), name = "Functional Loops") +
  
  # Common core residues annotation (y = -0.7 to -0.6)
  geom_rect(data = data.frame(position = common_core_residues),
            aes(xmin = position - 0.5, xmax = position + 0.5, ymin = -0.7, ymax = -0.6),
            fill = "#6D17A0", alpha = 0.5) +
  annotate("text", x = mean(range(common_core_residues)), y = -0.65, 
           label = "Core Residues", color = "black", size = 3, fontface = "bold") +
  
  # Secondary structure annotation - β-sheets and α-helices (y = -0.9 to -0.8)
  geom_rect(data = rects_sheet, 
            aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -0.9, ymax = -0.8),
            fill = "#FFE4A5", alpha = 0.7) +
  geom_text(data = rects_sheet, 
            aes(x = (xstart + xend)/2, y = -0.85, label = col),
            size = 3, color = "darkred", fontface = "bold") +
  
  geom_rect(data = rects_helix, 
            aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -0.9, ymax = -0.8),
            fill = "#97E9AD", alpha = 0.7) +
  geom_text(data = rects_helix, 
            aes(x = (xstart + xend)/2, y = -0.85, label = col),
            size = 3, color = "darkred", fontface = "bold") +
  
  # Axes and theme
  scale_x_continuous(breaks = seq(0, seq_length, by = 10), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-1.2, 0.5)) +
  labs(title = "KRAS Protein Sequence Annotation",
       x = "Amino Acid Position",
       y = "",
       caption = "Visualization of KRAS structural and functional features") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        plot.caption = element_text(hjust = 0.5, face = "italic"))

# Save the plot
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/Display _different_sequence_rectangles/results/20250902/sequence structure 250902.pdf", device = cairo_pdf,height = 6,width=20,dpi=300)

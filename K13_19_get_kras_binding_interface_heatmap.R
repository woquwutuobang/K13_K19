# KRAS Binding Interface Mutation Heatmaps
# This script generates heatmaps showing the effects of mutations on KRAS binding interfaces

library(ggplot2)
library(data.table)
library(dplyr)
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

# Function to plot binding interface mutation heatmap
plot_binding_interface_heatmap <- function(ddG, assay_sele, wt_aa, colour_scheme) {
  # Read and process ddG data
  ddG_data <- krasddpcams__read_ddG(ddG, assay_sele)
  ddG_data <- ddG_data[id != "WT", ]
  
  # Define amino acid order for consistent plotting
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  
  # Create heatmap template
  heatmap_template <- data.table(
    wt_codon = rep(unlist(strsplit(wt_aa, "")), each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = unlist(aa_list)
  )
  
  # Merge with ddG data
  input_heatmap <- merge(ddG_data, heatmap_template, 
                         by = c("wt_codon", "Pos_real", "mt_codon"), all = TRUE)
  
  # Set factor levels for amino acids
  input_heatmap[, mt_codon := factor(mt_codon, 
                                     levels = c("D", "E", "R", "H", "K", "S", "T", "N", "Q", 
                                                "C", "G", "P", "A", "V", "I", "L", "M", "F", "W", "Y"))]
  
  # Filter for binding interface residues
  binding_interface_residues <- c(88, 91, 87, 129, 90, 133, 94, 137, 95, 68, 136, 99, 102, 101, 107, 98)
  ddG_plot <- input_heatmap[Pos_real %in% binding_interface_residues, ]
  
  # Create position labels (e.g., "K88")
  ddG_plot[, wt_Pos := paste0(wt_codon, Pos_real)]
  
  # Set factor levels for position labels
  position_levels <- c("K88", "E91", "T87", "Q129", "F90", "L133", "H94", "Y137", 
                       "H95", "R68", "S136", "Q99", "R102", "K101", "E107", "E98")
  ddG_plot[, wt_Pos := factor(wt_Pos, levels = position_levels)]
  
  # Set wild-type values to zero
  ddG_plot[, assay := assay_sele]
  ddG_plot[wt_codon == mt_codon, `mean_kcal/mol` := 0]
  
  # Create heatmap plot
  ggplot() +
    geom_tile(
      data = ddG_plot,
      aes(y = wt_Pos, x = mt_codon, fill = `mean_kcal/mol`)
    ) +
    scale_fill_gradient2(
      low = colour_scheme[["blue"]],
      mid = "gray",
      high = colour_scheme[["red"]],
      na.value = "white",
      limits = c(-1, 3.3),
      oob = scales::squish,
      name = bquote(.(assay_sele) ~ Delta * Delta * "G (kcal/mol)")
    ) +
    geom_text(
      data = ddG_plot[Pos_real > 1 & wt_codon == mt_codon, ],
      aes(x = mt_codon, y = wt_Pos),
      label = "-",
      size = 5 * 5 / 14
    ) +
    xlab("Mutant Amino Acid") +
    ylab("Residue Position") +
    theme_classic2() +
    theme(
      text = element_text(size = 8),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, hjust = 0.5),
      axis.text = element_text(size = 8),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(3.1, "mm"),
      plot.margin = margin(0, 0, 0, 0),
      strip.background = element_rect(color = "white", fill = NULL)
    ) +
    coord_fixed()
}

# Generate and save K13 binding interface heatmap
k13_heatmap <- plot_binding_interface_heatmap(
  ddG = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13", 
  wt_aa = wt_aa, 
  colour_scheme = colour_scheme
)

ggsave("K13_binding_interface_heatmap.pdf", 
       plot = k13_heatmap, 
       device = cairo_pdf, 
       height = 6, 
       width = 6)

# Generate and save K19 binding interface heatmap  
k19_heatmap <- plot_binding_interface_heatmap(
  ddG = "path/to/weights_Binding_K19.txt",
  assay_sele = "K19",
  wt_aa = wt_aa,
  colour_scheme = colour_scheme
)

ggsave("K19_binding_interface_heatmap.pdf",
       plot = k19_heatmap,
       device = cairo_pdf,
       height = 6,
       width = 6)
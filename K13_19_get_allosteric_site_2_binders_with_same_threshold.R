# KRAS Allosteric Site Analysis
# This script analyzes the relationship between mutation effects and distance from binding partners

library(ggplot2)
library(ggpubr)
library(ggrepel)
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

# Secondary structure elements
secondary_structure <- data.frame(
  xstart = c(3, 15, 38, 51, 67, 77, 87, 109, 127, 139, 148),
  xend = c(9, 24, 44, 57, 73, 84, 104, 115, 136, 143, 166),
  col = c("β1", "α1", "β2", "β3", "α2", "β4", "α3", "β5", "α4", "β6", "α5")
)

# File paths for ddG data
ddG_k13 <- "path/to/weights_Binding_K13.txt"
assay_k13 <- "K13"
ddG_k19 <- "path/to/weights_Binding_K19.txt"
assay_k19 <- "K19"

# Load annotation data
anno <- fread("path/to/annotation_data.csv")
anno[, Pos_real := Pos]

# Calculate weighted mean ddG for each assay
weighted_mean_ddG <- list()
weighted_mean_ddG[[assay_k13]] <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(
  ddG = ddG_k13,
  assay_sele = assay_k13
)

weighted_mean_ddG[[assay_k19]] <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(
  ddG = ddG_k19,
  assay_sele = assay_k19
)

# Combine data from both assays
assay_list <- list(assay_k13, assay_k19)
data_plot <- data.table()

for (assay_i in assay_list) {
  data_plot_assay <- merge(weighted_mean_ddG[[assay_i]], anno, by = "Pos_real", all = TRUE)
  
  # Classify binding types
  data_plot_assay[, binding_type := "allosteric_site"]
  data_plot_assay[get(paste0("scHAmin_ligand_", assay_i)) < 5, binding_type := "binding_site"]
  
  data_plot_assay[, binding_type_gtp_included := binding_type]
  data_plot_assay[GXPMG_scHAmin_ligand_RAF1 < 5, 
                  binding_type_gtp_included := "gtp_binding_site"]
  
  data_plot <- rbind(data_plot, data_plot_assay)
}

# Calculate regulatory threshold
reg_threshold <- data_plot[binding_type == "binding_site",
                           sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                           .SDcols = c("mean", "sigma")]

# Classify site types
data_plot[, site_type := "other"]
data_plot[binding_type_gtp_included == "binding_site", site_type := "binding_interface"]
data_plot[binding_type_gtp_included == "gtp_binding_site", site_type := "other_gtp_pocket"]
data_plot[binding_type_gtp_included == "gtp_binding_site" & mean > reg_threshold & 
            binding_type != "binding_site" & count > 9.5, 
          site_type := "allosteric_gtp_pocket"]
data_plot[binding_type_gtp_included == "allosteric_site" & mean > reg_threshold & count > 9.5, 
          site_type := "major_allosteric"]

# Map secondary structure information
rects_dt <- as.data.table(secondary_structure)

# Assign beta strands
data_plot[Pos_real >= rects_dt[col == "β1", xstart] & Pos_real <= rects_dt[col == "β1", xend], 
          colors_type := "β1"]
data_plot[Pos_real >= rects_dt[col == "β2", xstart] & Pos_real <= rects_dt[col == "β2", xend], 
          colors_type := "β2"]
data_plot[Pos_real >= rects_dt[col == "β3", xstart] & Pos_real <= rects_dt[col == "β3", xend], 
          colors_type := "β3"]
data_plot[Pos_real >= rects_dt[col == "β4", xstart] & Pos_real <= rects_dt[col == "β4", xend], 
          colors_type := "β4"]
data_plot[Pos_real >= rects_dt[col == "β5", xstart] & Pos_real <= rects_dt[col == "β5", xend], 
          colors_type := "β5"]
data_plot[Pos_real >= rects_dt[col == "β6", xstart] & Pos_real <= rects_dt[col == "β6", xend], 
          colors_type := "β6"]

# Assign alpha helices
data_plot[Pos_real >= rects_dt[col == "α1", xstart] & Pos_real <= rects_dt[col == "α1", xend], 
          colors_type := "α1"]
data_plot[Pos_real >= rects_dt[col == "α2", xstart] & Pos_real <= rects_dt[col == "α2", xend], 
          colors_type := "α2"]
data_plot[Pos_real >= rects_dt[col == "α3", xstart] & Pos_real <= rects_dt[col == "α3", xend], 
          colors_type := "α3"]
data_plot[Pos_real >= rects_dt[col == "α4", xstart] & Pos_real <= rects_dt[col == "α4", xend], 
          colors_type := "α4"]
data_plot[Pos_real >= rects_dt[col == "α5", xstart] & Pos_real <= rects_dt[col == "α5", xend], 
          colors_type := "α5"]

# Classify secondary structure shapes
data_plot[, shape := "other"]
data_plot[colors_type %chin% c("β1", "β2", "β3", "β4", "β5", "β6"), shape := "beta_strand"]
data_plot[colors_type %chin% c("α1", "α2", "α3", "α4", "α5"), shape := "alpha_helix"]

# Filter and factor data
data_plot <- data_plot[Pos_real > 1 & count > 9.5, ]
data_plot[, site_type := factor(site_type,
                                levels = c("binding_interface", "allosteric_gtp_pocket", 
                                           "other_gtp_pocket", "major_allosteric", "other"))]

data_plot[, assay := factor(assay, levels = c("K13", "K19"))]

# Calculate distances for each assay
allosteric_sites <- list()
for (assay_i in assay_list) {
  data_plot[assay == assay_i, distance_bp := get(paste0("scHAmin_ligand_", assay_i))]
  allosteric_sites[[assay_i]] <- list(
    binding_interface = data_plot[binding_type == "binding_site" & assay == assay_i, Pos_real],
    allosteric_gtp_pocket = data_plot[site_type == "allosteric_gtp_pocket" & assay == assay_i, Pos_real],
    other_gtp_pocket = data_plot[site_type == "other_gtp_pocket" & assay == assay_i, Pos_real],
    major_allosteric = data_plot[site_type == "major_allosteric" & assay == assay_i, Pos_real]
  )
}

# Create the plot
allosteric_plot <- ggplot() +
  # Points with error bars
  geom_point(data = data_plot,
             aes(x = distance_bp, y = mean, color = site_type, shape = shape),
             size = 0.35) +
  geom_pointrange(data = data_plot,
                  aes(x = distance_bp, y = mean, ymin = mean - sigma, ymax = mean + sigma,
                      color = site_type, shape = shape),
                  size = 0.35) +
  
  # Reference lines
  geom_hline(yintercept = reg_threshold, linetype = 2, size = 0.1) +
  geom_vline(xintercept = 5, linetype = 2, size = 0.1) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.1) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.1) +
  
  # Label allosteric sites
  geom_text_repel(data = data_plot[site_type == "major_allosteric", ],
                  aes(x = distance_bp, y = mean, label = Pos_real),
                  nudge_y = 0.05, color = colour_scheme[["orange"]], size = 5 * 0.35) +
  geom_text_repel(data = data_plot[site_type == "allosteric_gtp_pocket", ],
                  aes(x = distance_bp, y = mean, label = Pos_real),
                  nudge_y = 0.05, color = colour_scheme[["blue"]], size = 5 * 0.35) +
  
  # Styling
  scale_color_manual(
    values = c("red", "blue", "lightblue", "orange", "gray"),
    labels = c("Binding Interface", "Allosteric GTP Pocket", 
               "Other GTP Pocket", "Major Allosteric", "Others")
  ) +
  xlab("Distance to Binding Partner (Å)") +
  ylab("Weighted Mean |ΔΔG| (kcal/mol)") +
  labs(color = "Site Type", shape = "Secondary Structure") +
  facet_wrap(~ assay, ncol = 3) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    text = element_text(size = 9),
    legend.position = "right",
    strip.text = element_text(size = 9),
    strip.background = element_rect(colour = "white", fill = "white"),
    panel.spacing = unit(0.2, "mm"),
    legend.text = element_text(size = 5),
    plot.margin = margin(0, 1, 0, 1, "mm"),
    legend.margin = margin(0, 0, 0, -2, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.key.height = unit(4, "mm")
  )

# Display and save plot
print(allosteric_plot)
ggsave("K13_K19_allosteric_analysis.pdf",
       plot = allosteric_plot,
       device = cairo_pdf,
       height = 5,
       width = 8,
       dpi = 300)
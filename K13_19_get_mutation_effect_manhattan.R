# KRAS-DARPins Mutation Effect Manhattan Plots
# This script generates Manhattan plots showing the effects of mutations on K13 and K19 binding

library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)

# Wild-type KRAS sequence (residues 2-189)
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Secondary structure elements for annotation
beta_sheets <- data.frame(
  xstart = c(3, 38, 51, 77, 109, 139),
  xend = c(9, 44, 57, 84, 115, 143),
  col = c("β1", "β2", "β3", "β4", "β5", "β6")
)

alpha_helices <- data.frame(
  xstart = c(15, 67, 87, 127, 148),
  xend = c(24, 73, 104, 136, 166),
  col = c("α1", "α2", "α3", "α4", "α5")
)

# Function to create Manhattan plot for mutation effects
plot_mutation_effects <- function(input, assay_sele, anno, beta_sheets, alpha_helices, wt_aa) {
  # Read and process ddG data
  ddG <- fread(input)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", wt_codon := substr(id, 1, 1)]
  ddG[id != "WT", mt_codon := substr(id, nchar(id), nchar(id))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  
  # Create heatmap template for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_template <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]]))
  )
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_template, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG per position
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  anno_data <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  # Classify binding sites
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding_site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get(paste0("GXPMG_scHAmin_ligand_", assay_sele)) < 5, 
            binding_type_gtp_included := "gtp_binding_site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding_site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "other"]
  data_plot[binding_type_gtp_included == "binding_site", site_type := "binding_interface"]
  data_plot[binding_type_gtp_included == "gtp_binding_site", site_type := "gtp_binding_interface"]
  
  # Classify mutation types
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "other"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # Detailed mutation classification
  data_plot_mutation[Pos %in% data_plot[site_type == "binding_interface", Pos] &
                       allosteric_mutation == TRUE, 
                     mutation_type := "orthosteric_large_effect"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "binding_interface", Pos] &
                       allosteric_mutation == FALSE, 
                     mutation_type := "orthosteric_small_effect"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "gtp_binding_interface", Pos] &
                       allosteric_mutation == TRUE, 
                     mutation_type := "gtp_allosteric"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "gtp_binding_interface", Pos] &
                       allosteric_mutation == FALSE, 
                     mutation_type := "gtp_other"]
  
  data_plot_mutation[!site_type %in% c("gtp_binding_interface", "binding_interface") &
                       allosteric_mutation == TRUE, 
                     mutation_type := "allosteric"]
  
  data_plot_mutation[!site_type %in% c("gtp_binding_interface", "binding_interface") &
                       allosteric_mutation == FALSE, 
                     mutation_type := "other"]
  
  # Set factor levels for consistent coloring
  data_plot_mutation[, mutation_type := factor(mutation_type,
                                               levels = c("orthosteric_large_effect", "orthosteric_small_effect",
                                                          "gtp_allosteric", "gtp_other", "allosteric", "other"))]
  
  # Create the plot
  p <- ggplot() +
    # Secondary structure annotations
    geom_rect(data = beta_sheets, 
              aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5),
              fill = "#FFE4A5", alpha = 0.15) +
    geom_rect(data = alpha_helices, 
              aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5),
              fill = "#97E9AD", alpha = 0.15) +
    
    # Mutation points
    geom_point(data = data_plot_mutation,
               aes(x = Pos_real, y = `mean_kcal/mol`, color = mutation_type), 
               size = 0.9) +
    
    # Styling
    scale_color_manual(values = c("#F4270C", "#FFB0A5", "#1B38A6", "#9DACE3", "#F4AD0C", "gray")) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(expand = c(1 / 188, 11 / 188)) +
    ylab(paste0("Binding ΔΔG (", assay_sele, ")")) +
    xlab("Amino Acid Position") +
    labs(color = "Mutation Type") +
    
    # Secondary structure labels
    annotate("text", x = (beta_sheets$xstart[1] + beta_sheets$xend[1]) / 2, y = 3.1, 
             label = "strand", size = 2.5, vjust = 0, fontface = "bold") +
    annotate("text", x = (alpha_helices$xstart[1] + alpha_helices$xend[1]) / 2, y = 3.1, 
             label = "helix", size = 2.5, vjust = 0, fontface = "bold") +
    
    theme_classic2() +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.position = "none",
          strip.background = element_rect(colour = "black", fill = "white")) +
    coord_fixed(ratio = 10, xlim = c(-0.5, 190), ylim = c(-1.5, 3.5))
  
  return(p)
}

# Function to create combined K13 and K19 plots
plot_combined_darpins <- function() {
  # Generate individual plots
  p_k13 <- plot_mutation_effects(
    input = "path/to/weights_Binding_K13.txt",
    assay_sele = "K13",
    anno = "path/to/annotation_data.csv",
    beta_sheets = beta_sheets,
    alpha_helices = alpha_helices,
    wt_aa = wt_aa
  )
  
  p_k19 <- plot_mutation_effects(
    input = "path/to/weights_Binding_K19.txt",
    assay_sele = "K19",
    anno = "path/to/annotation_data.csv",
    beta_sheets = beta_sheets,
    alpha_helices = alpha_helices,
    wt_aa = wt_aa
  )
  
  # Combine plots vertically
  combined_plot <- p_k13 / p_k19 + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "none")
  
  print(combined_plot)
  return(combined_plot)
}

# Generate and save combined plot
combined_plot <- plot_combined_darpins()
ggsave("K13_K19_mutation_effects.pdf",
       plot = combined_plot,
       device = cairo_pdf,
       height = 8,
       width = 10)
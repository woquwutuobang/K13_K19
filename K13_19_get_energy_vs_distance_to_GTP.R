library(krasddpcams)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(ggrepel)

##=======================================
##Put the output images together into one figure + Their respective figure
##========================================

# Wild-type KRAS sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Input files and assay names
input_files <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt", 
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt"
)

assay_names <- c("K13", "K19", "RAF1", "K55", "K27", "RALGDS", "PI3KCG", "SOS1")

# Annotation file with structural information
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
anno <- fread(anno_file)

# Output directory
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s6/20251011/energy_vs_distance_to_GTP"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Color scheme for site types
site_type_colors <- c(
  "Binding interface site"   = "#F4270C",   # Red
  "Other GTP pocket site"    = "#75C2F6",   # Light blue  
  "Allosteric GTP pocket site" = "#1B38A6", # Dark blue
  "Major allosteric site"    = "#F4AD0C",   # Orange
  "Reminder"                 = "grey70"       # Grey
)

# ================================
# Main Analysis Function
# ================================

analyze_allosteric_sites <- function(input_files, assay_names, anno) {
  stopifnot(length(input_files) == length(assay_names))
  
  # ---- Step 1: Data Preparation ----
  prepare_assay_data <- function(ddG_file, assay_name, anno) {
    cat("Processing assay:", assay_name, "\n")
    
    # Calculate weighted mean absolute ddG values
    weighted_mean_ddG <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(
      ddG = ddG_file,
      assay_sele = assay_name
    )
    
    # Merge with structural annotation
    data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
    
    # Classify binding sites based on distance cutoff (5Å)
    data_plot[, binding_type := "allosteric site"]
    sc_col <- paste0("scHAmin_ligand_", assay_name)
    if (sc_col %in% names(data_plot)) {
      data_plot[get(sc_col) < 5, binding_type := "binding site"]
    } else {
      warning("Column ", sc_col, " not found, cannot classify binding sites")
    }
    return(data_plot)
  }
  
  # ---- Step 2: Load and Process Data ----
  assay_data <- lapply(seq_along(input_files), function(i) {
    prepare_assay_data(input_files[i], assay_names[i], anno)
  })
  names(assay_data) <- assay_names
  
  # ---- Step 3: Calculate Classification Thresholds ----
  calculate_binding_threshold <- function(data_plot) {
    data_plot[binding_type == "binding site",
              sum(abs(mean)/sigma^2, na.rm = TRUE) / sum(1/sigma^2, na.rm = TRUE)]
  }
  
  # Different thresholds for DARPins vs other binders
  threshold_darpins <- mean(sapply(c("K13", "K19"), function(a) 
    calculate_binding_threshold(assay_data[[a]])))
  threshold_others <- mean(sapply(setdiff(assay_names, c("K13", "K19")),
                                  function(a) calculate_binding_threshold(assay_data[[a]])))
  
  cat("Threshold for DARPins (K13/K19):", threshold_darpins, "\n")
  cat("Threshold for other binders:", threshold_others, "\n")
  
  # ---- Step 4: Site Type Classification ----
  classify_site_types <- function(data_plot, assay_name, reg_threshold) {
    data_plot[, binding_type_gtp_included := binding_type]
    
    # Identify GTP binding sites
    if ("GXPMG_scHAmin_ligand_RAF1" %in% names(data_plot)) {
      data_plot[GXPMG_scHAmin_ligand_RAF1 < 5, binding_type_gtp_included := "GTP binding site"]
    }
    
    # Initialize all sites as "Reminder"
    data_plot[, site_type := "Reminder"]
    
    # Binding interface sites (direct contacts)
    data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
    
    # GTP pocket sites
    data_plot[binding_type_gtp_included == "GTP binding site", site_type := "Other GTP pocket site"]
    
    # Allosteric GTP pocket sites (GTP pocket with allosteric effects)
    sc_col <- paste0("scHAmin_ligand_", assay_name)
    data_plot[binding_type_gtp_included == "GTP binding site" &
                mean > reg_threshold &
                (!is.null(data_plot[[sc_col]]) & get(sc_col) >= 5),
              site_type := "Allosteric GTP pocket site"]
    
    # Major allosteric sites (distant from binding interface with large effects)
    data_plot[binding_type_gtp_included == "allosteric site" & mean > reg_threshold,
              site_type := "Major allosteric site"]
    
    return(data_plot)
  }
  
  # ---- Step 5: Combine all data for faceted plot ----
  combined_data <- data.frame()
  
  for (i in seq_along(assay_data)) {
    assay_name <- assay_names[i]
    reg_threshold <- ifelse(assay_name %in% c("K13", "K19"), threshold_darpins, threshold_others)
    data_plot <- classify_site_types(assay_data[[i]], assay_name, reg_threshold)
    
    # Add assay name and threshold to data
    data_plot$assay <- assay_name
    data_plot$reg_threshold <- reg_threshold
    
    # Combine data
    combined_data <- rbind(combined_data, data_plot)
    
    # Update assay_data with classified sites
    assay_data[[i]] <- data_plot
  }
  
  # ---- Step 6: Create combined faceted plot ----
  combined_plot <- ggplot(data = combined_data,
                          aes(x = GXPMG_scHAmin_ligand_RAF1,  # Distance to GTP pocket
                              y = mean,                       # Weighted mean |ddG|
                              color = site_type)) +
    geom_point(size = 1.5, alpha = 0.8) +
    
    # Add threshold lines
    geom_hline(aes(yintercept = reg_threshold), linetype = 2, size = 0.3) +
    geom_vline(xintercept = 5, linetype = 2, size = 0.3) +
    
    # Label major allosteric sites
    geom_text_repel(data = combined_data[combined_data$site_type == "Major allosteric site", ],
                    aes(label = Pos), nudge_y = 0.05,
                    color = site_type_colors["Major allosteric site"],
                    size = 2.5, max.overlaps = 20, segment.size = 0.2) +
    
    # Label allosteric GTP pocket sites  
    geom_text_repel(data = combined_data[combined_data$site_type == "Allosteric GTP pocket site", ],
                    aes(label = Pos), nudge_y = 0.05,
                    color = site_type_colors["Allosteric GTP pocket site"],
                    size = 2.5, max.overlaps = 20, segment.size = 0.2) +
    
    # Color scheme
    scale_color_manual(values = site_type_colors) +
    
    # Axis labels
    labs(x = "Distance to GTP Pocket (Å)",
         y = "Weighted mean |ΔΔG| (kcal/mol)",
         color = "Site Type",
         title = "Allosteric Site Analysis for All Binders") +
    
    # Facet by assay with 4 columns - each panel has its own axes
    facet_wrap(~ assay, ncol = 4, scales = "fixed") +
    
    # Theme settings
    theme_classic2() +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.line = element_line(size = 0.5),  
      axis.ticks = element_line(size = 0.5), 
      text = element_text(size = 8),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      plot.margin = margin(5, 5, 5, 5, "mm"),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      legend.spacing.y = unit(2, "mm"),
      legend.key.height = unit(4, "mm"),
      plot.title = element_text(hjust = 0.5, size = 8),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5)
    ) +
    coord_cartesian(xlim = c(0, 34), ylim = c(0, 2.3))
  
  # Save combined plot
  combined_pdf <- file.path(output_dir, "allosteric_analysis_combined.pdf")
  ggsave(combined_pdf, plot = combined_plot, width = 14, height = 8, units = "in", device = cairo_pdf)
  cat("Saved combined plot:", combined_pdf, "\n")
  
  # ---- Step 7: Also save individual plots (optional) ----
  for (i in seq_along(assay_data)) {
    assay_name <- assay_names[i]
    data_plot <- assay_data[[i]]
    
    # Create individual plot
    p <- ggplot(data = data_plot,
                aes(x = GXPMG_scHAmin_ligand_RAF1,
                    y = mean,
                    color = site_type)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_hline(yintercept = ifelse(assay_name %in% c("K13", "K19"), threshold_darpins, threshold_others), 
                 linetype = 2, size = 0.3) +
      geom_vline(xintercept = 5, linetype = 2, size = 0.3) +
      geom_text_repel(data = data_plot[site_type == "Major allosteric site", ],
                      aes(label = Pos), nudge_y = 0.05,
                      color = site_type_colors["Major allosteric site"],
                      size = 2.5, max.overlaps = 20) +
      geom_text_repel(data = data_plot[site_type == "Allosteric GTP pocket site", ],
                      aes(label = Pos), nudge_y = 0.05,
                      color = site_type_colors["Allosteric GTP pocket site"],
                      size = 2.5, max.overlaps = 20) +
      scale_color_manual(values = site_type_colors) +
      labs(x = "Distance to GTP Pocket (Å)",
           y = "Weighted mean |ΔΔG| (kcal/mol)",
           color = "Site Type",
           title = paste("Allosteric Site Analysis -", assay_name)) +
      theme_classic2() +
      theme(
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        text = element_text(size = 8),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
      ) +
      coord_cartesian(xlim = c(0, 34), ylim = c(0, 2.3))
    
    # Save individual plot
    pdf_file <- file.path(output_dir, paste0("allosteric_analysis_", assay_name, ".pdf"))
    ggsave(pdf_file, plot = p, width = 6, height = 4, units = "in", device = cairo_pdf)
    cat("Saved:", pdf_file, "\n")
  }
  
  return(list(assay_data = assay_data, combined_data = combined_data))
}

# ================================
# Execute Analysis
# ================================
result_data <- analyze_allosteric_sites(input_files, assay_names, anno)

# Print analysis summary
cat("\nAllosteric Site Analysis Summary:\n")
cat("==================================\n")
cat("Assays analyzed:", paste(assay_names, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")
cat("Site classification complete for all assays\n")

# Additional summary statistics
for (assay in assay_names) {
  data <- result_data[[assay]]
  cat("\n", assay, "Site Distribution:\n")
  print(table(data$site_type))

}

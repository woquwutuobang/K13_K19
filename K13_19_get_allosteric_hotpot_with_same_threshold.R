# KRAS Binding Site Classification Analysis
# This script classifies binding sites and allosteric sites across multiple KRAS binding partners

library(krasddpcams)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

# Define wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Define input file paths and assay names
input_files <- c(
  "path/to/weights_Binding_K13.txt",
  "path/to/weights_Binding_K19.txt",
  "path/to/weights_Binding_RAF1.txt",
  "path/to/weights_Binding_K55.txt",
  "path/to/weights_Binding_K27.txt",
  "path/to/weights_Binding_RAL.txt",
  "path/to/weights_Binding_PI3.txt",
  "path/to/weights_Binding_SOS.txt"
)

assay_names <- c("K13", "K19", "RAF1", "K55", "K27", "RALGDS", "PI3KCG", "SOS1")

# Define annotation file path
anno_file <- "path/to/annotations.csv"

# Read annotation data
annotations <- fread(anno_file)

#' Classify binding sites and allosteric sites across multiple assays
#'
#' @param input_files Vector of file paths to ddG data
#' @param assay_names Vector of assay names
#' @param annotations Annotation data table
#' @return List of classified data for each assay

classify_binding_sites <- function(input_files, assay_names, annotations) {
  stopifnot(length(input_files) == length(assay_names))
  
  # Internal function: prepare data for single assay
  prepare_assay_data <- function(ddG_file, assay_name, anno) {
    cat("Processing assay:", assay_name, "\n")
    
    # Calculate weighted mean ΔΔG
    weighted_mean_ddG <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(
      ddG = ddG_file,
      assay_sele = assay_name
    )
    
    # Merge with annotations
    data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
    
    # Classify binding sites
    data_plot[, binding_type := "allosteric site"]
    distance_column <- paste0("scHAmin_ligand_", assay_name)
    
    if (distance_column %in% names(data_plot)) {
      data_plot[get(distance_column) < 5, binding_type := "binding site"]
    } else {
      warning("Column ", distance_column, " not found, skipping binding site classification")
    }
    
    return(data_plot)
  }
  
  # Process all assays
  assay_data <- lapply(seq_along(input_files), function(i) {
    prepare_assay_data(input_files[i], assay_names[i], annotations)
  })
  names(assay_data) <- assay_names
  
  # Calculate thresholds
  calculate_threshold <- function(data_plot) {
    result <- data_plot[binding_type == "binding site",
                        sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE) /
                          sum(1/.SD[[2]]^2, na.rm = TRUE),
                        .SDcols = c("mean", "sigma")]
    return(result)
  }
  
  # Calculate thresholds for different assay groups
  threshold_darpins <- mean(sapply(c("K13", "K19"), function(a) calculate_threshold(assay_data[[a]])))
  threshold_effectors <- mean(sapply(setdiff(assay_names, c("K13", "K19")), 
                                     function(a) calculate_threshold(assay_data[[a]])))
  
  cat("Threshold for DARPins (K13/K19):", threshold_darpins, "\n")
  cat("Threshold for effectors:", threshold_effectors, "\n")
  
  # Site type classification function
  classify_site_types <- function(data_plot, assay_name, reg_threshold) {
    data_plot[, binding_type_gtp_included := binding_type]
    
    # Include GTP binding sites
    if ("GXPMG_scHAmin_ligand_RAF1" %in% names(data_plot)) {
      data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, 
                binding_type_gtp_included := "GTP binding site"]
    }
    
    # Initialize site types
    data_plot[, site_type := "other"]
    data_plot[binding_type_gtp_included == "binding site", 
              site_type := "binding_interface"]
    data_plot[binding_type_gtp_included == "GTP binding site", 
              site_type := "gtp_binding_interface"]
    
    # Classify allosteric sites
    distance_column <- paste0("scHAmin_ligand_", assay_name)
    data_plot[binding_type_gtp_included == "GTP binding site" &
                mean > reg_threshold &
                (!is.null(data_plot[[distance_column]]) & get(distance_column) >= 5),
              site_type := "allosteric_gtp_pocket"]
    
    data_plot[binding_type_gtp_included == "allosteric site" & mean > reg_threshold,
              site_type := "major_allosteric"]
    
    return(data_plot)
  }
  
  # Apply classification and print results
  final_results <- list()
  
  for (assay in assay_names) {
    # Select appropriate threshold
    reg_threshold <- if (assay %in% c("K13", "K19")) threshold_darpins else threshold_effectors
    final_results[[assay]] <- classify_site_types(assay_data[[assay]], assay, reg_threshold)
    
    # Print classification results
    cat("\n---", assay, "---\n")
    
    site_types <- c("binding_interface", "allosteric_gtp_pocket", "major_allosteric")
    site_labels <- c("Binding interface", "Allosteric GTP pocket", "Major allosteric")
    
    for (i in seq_along(site_types)) {
      positions <- final_results[[assay]][site_type == site_types[i], Pos]
      if (length(positions) > 0) {
        cat(site_labels[i], ":", paste(positions, collapse = ", "), "\n")
      } else {
        cat(site_labels[i], ": none\n")
      }
    }
  }
  
  return(final_results)
}

# Execute analysis
classification_results <- classify_binding_sites(input_files, assay_names, annotations)


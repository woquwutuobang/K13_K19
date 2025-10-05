# KRAS Mutation Effect Heatmap Analysis
# This script analyzes and visualizes high-effect mutations across multiple KRAS binding partners

library(krasddpcams)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(ggnewscale)

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

# Define annotation files and core surface file paths
anno_file <- "path/to/anno_final_for_8.csv"
core_surface_file <- "path/to/KRAS_WT_166_monomer_get_rasa_20250701_2.csv"
dist_cutoff <- 5

# Read annotation data
anno <- fread(anno_file)
core_surface <- fread(core_surface_file)
anno <- merge(core_surface, anno, by = "Pos", all = TRUE)

# Define SASA thresholds for core and surface residues
core_threshold_sasa <- 0.25
surface_threshold_sasa <- 0.26

# Classify core and surface residues
anno <- anno %>%
  mutate(type = case_when(
    RASA <= core_threshold_sasa ~ "core",
    RASA >= surface_threshold_sasa ~ "surface"
  )) %>%
  filter(!is.na(type))

# Read and process ddG data for all assays
read_and_process_ddG <- function(input_path, assay_name) {
  ddG <- krasddpcams__read_ddG(input_path, assay_name)
  ddG <- ddG[, c(1:4, 19:20, 23:24)]
  colnames(ddG)[5:8] <- paste0(colnames(ddG)[5:8], "_", assay_name)
  return(ddG)
}

# Read data for all assays
ddG_list <- lapply(seq_along(input_files), function(i) {
  read_and_process_ddG(input_files[i], assay_names[i])
})

# Merge data from all assays
data_plot_mutation <- Reduce(function(x, y) {
  merge(x, y, by = c("Pos_real", "id", "wt_codon", "mt_codon"), all = TRUE)
}, ddG_list)

# Add annotation information
anno[, Pos_real := Pos]
data_plot_mutation <- merge(data_plot_mutation, anno, by = "Pos_real", all = TRUE)

# Prepare data and calculate thresholds
prepare_data_plot_mutation <- function(data_plot_mutation, dist_cutoff = 5) {
  dp <- copy(data_plot_mutation)
  
  # Define assay groups
  assay_group1 <- c("K13", "K19")
  assay_group2 <- c("RAF1", "RALGDS", "K55", "K27", "PI3KCG", "SOS1")
  
  thresholds <- list()
  
  # Create binding site classification for each assay
  for(assay in c(assay_group1, assay_group2)) {
    # Create binding site type
    dp[, paste0("binding_type_", assay) := "allosteric site"]
    dp[get(paste0("scHAmin_ligand_", assay)) < dist_cutoff, 
       paste0("binding_type_", assay) := "binding site"]
    
    # Create classification including GTP binding sites
    dp[, paste0("binding_type_gtp_included_", assay) := get(paste0("binding_type_", assay))]
    dp[get(paste0("GXPMG_scHAmin_ligand_RAF1")) < dist_cutoff, 
       paste0("binding_type_gtp_included_", assay) := "GTP binding site"]
    
    # Create final site type classification
    dp[, paste0("site_type_", assay) := "Reminder"]
    dp[get(paste0("binding_type_gtp_included_", assay)) == "binding site", 
       paste0("site_type_", assay) := "Binding interface site"]
    dp[get(paste0("binding_type_gtp_included_", assay)) == "GTP binding site", 
       paste0("site_type_", assay) := "GTP binding interface site"]
  }
  
  # Calculate thresholds for group 1 assays
  threshold1_data <- dp[binding_type_K13 == "binding site" | binding_type_K19 == "binding site"]
  if (nrow(threshold1_data) > 0) {
    mean_values <- c(threshold1_data$mean_K13[!is.na(threshold1_data$mean_K13)], 
                     threshold1_data$mean_K19[!is.na(threshold1_data$mean_K19)])
    sigma_values <- c(threshold1_data$std_K13[!is.na(threshold1_data$std_K13)], 
                      threshold1_data$std_K19[!is.na(threshold1_data$std_K19)])
    
    if (length(mean_values) > 0 && length(sigma_values) > 0) {
      thresholds[["threshold_group1"]] <- sum(abs(mean_values) / sigma_values^2, na.rm = TRUE) / 
        sum(1 / sigma_values^2, na.rm = TRUE)
    } else {
      thresholds[["threshold_group1"]] <- NA
    }
  } else {
    thresholds[["threshold_group1"]] <- NA
  }
  
  # Calculate thresholds for group 2 assays
  threshold2_conditions <- paste0("binding_type_", assay_group2, " == 'binding site'", collapse = " | ")
  threshold2_data <- dp[eval(parse(text = threshold2_conditions))]
  
  if (nrow(threshold2_data) > 0) {
    mean_values <- numeric(0)
    sigma_values <- numeric(0)
    
    for(assay in assay_group2) {
      mean_col <- paste0("mean_", assay)
      std_col <- paste0("std_", assay)
      
      if (mean_col %in% names(threshold2_data) && std_col %in% names(threshold2_data)) {
        valid_indices <- !is.na(threshold2_data[[mean_col]]) & !is.na(threshold2_data[[std_col]])
        mean_values <- c(mean_values, threshold2_data[[mean_col]][valid_indices])
        sigma_values <- c(sigma_values, threshold2_data[[std_col]][valid_indices])
      }
    }
    
    if (length(mean_values) > 0 && length(sigma_values) > 0) {
      thresholds[["threshold_group2"]] <- sum(abs(mean_values) / sigma_values^2, na.rm = TRUE) / 
        sum(1 / sigma_values^2, na.rm = TRUE)
    } else {
      thresholds[["threshold_group2"]] <- NA
    }
  } else {
    thresholds[["threshold_group2"]] <- NA
  }
  
  list(data_plot = dp, thresholds = thresholds)
}

# Process data and get thresholds
result <- prepare_data_plot_mutation(data_plot_mutation, dist_cutoff)
data_plot_mutation_processed <- result$data_plot
threshold_group1 <- result$thresholds$threshold_group1
threshold_group2 <- result$thresholds$threshold_group2

# Count high-effect mutations function
count_high_effect <- function(data, assays, thresholds) {
  result <- data.table(Pos_real = unique(data$Pos_real))
  
  # Add type information
  type_info <- unique(data[, .(Pos_real, type, RASA)])
  result <- merge(result, type_info, by = "Pos_real", all.x = TRUE)
  
  for(i in seq_along(assays)) {
    assay <- assays[i]
    mean_col <- paste0("mean_", assay)
    
    if(!mean_col %in% names(data)) {
      warning(paste("Column", mean_col, "not found in data"))
      result[, paste0("count_", assay) := NA]
      next
    }
    
    threshold <- ifelse(assay %in% c("K13", "K19"), thresholds$group1, thresholds$group2)
    
    # Count mutations exceeding threshold for each position
    count_data <- data[, .(count = sum(get(mean_col) >= threshold, na.rm = TRUE)), by = Pos_real]
    result[, paste0("count_", assay) := count_data$count[match(Pos_real, count_data$Pos_real)]]
  }
  
  return(result)
}

assays_group1 <- c("K13", "K19")
assays_group2 <- c("RAF1", "RALGDS", "K55", "K27", "PI3KCG", "SOS1")
all_assays <- c(assays_group1, assays_group2)

# Create thresholds list
thresholds_list <- list(group1 = threshold_group1, group2 = threshold_group2)

# Run function
heatmap_data <- count_high_effect(data_plot_mutation_processed, all_assays, thresholds_list)

# Ensure type information is correct
heatmap_data[, type := ifelse(RASA <= core_threshold_sasa, "core", 
                              ifelse(RASA >= surface_threshold_sasa, "surface", NA))]
heatmap_data <- heatmap_data[!is.na(type)]

# Statistical analysis function
count_high_effect_with_stats <- function(data, assays, thresholds) {
  assay_stats <- list()
  
  for(assay in assays) {
    # Determine threshold
    threshold <- ifelse(assay %in% c("K13", "K19"), thresholds$group1, thresholds$group2)
    
    # Get relevant columns
    mean_col <- paste0("mean_", assay)
    site_type_col <- paste0("site_type_", assay)
    
    if(!mean_col %in% names(data) || !site_type_col %in% names(data)) {
      warning(paste("Required columns for", assay, "not found"))
      next
    }
    
    # Filter valid data
    assay_data <- data[!is.na(get(mean_col))]
    
    # Count interface mutations
    interface_data <- assay_data[get(site_type_col) %in% c("Binding interface site", "GTP binding interface site")]
    n_interface <- sum(interface_data[[mean_col]] >= threshold, na.rm = TRUE)
    
    # Count non-interface mutations
    non_interface_data <- assay_data[!get(site_type_col) %in% c("Binding interface site", "GTP binding interface site")]
    m_non_interface <- sum(non_interface_data[[mean_col]] >= threshold, na.rm = TRUE)
    
    # Further classify non-interface mutations into core and surface
    non_interface_core <- non_interface_data[type == "core"]
    m_core <- sum(non_interface_core[[mean_col]] >= threshold, na.rm = TRUE)
    
    non_interface_surface <- non_interface_data[type == "surface"]
    m_surface <- sum(non_interface_surface[[mean_col]] >= threshold, na.rm = TRUE)
    
    # Background distribution: all mutations
    total_mutations <- nrow(assay_data)
    background_interface <- nrow(assay_data[get(site_type_col) %in% c("Binding interface site", "GTP binding interface site")])
    background_non_interface <- total_mutations - background_interface
    
    # Significance test - interface mutations vs background
    if(background_interface > 0 && total_mutations > 0) {
      n_p_value <- binom.test(n_interface, total_mutations, 
                              p = background_interface/total_mutations)$p.value
    } else {
      n_p_value <- NA
    }
    
    # Significance test - non-interface mutations vs background
    if(background_non_interface > 0 && total_mutations > 0) {
      m_p_value <- binom.test(m_non_interface, total_mutations, 
                              p = background_non_interface/total_mutations)$p.value
    } else {
      m_p_value <- NA
    }
    
    # Fisher's test for core vs surface in non-interface mutations
    if(m_core > 0 || m_surface > 0) {
      # Background distribution: core and surface proportions in non-interface regions
      background_core_non_interface <- nrow(assay_data[!get(site_type_col) %in% c("Binding interface site", "GTP binding interface site") & type == "core"])
      background_surface_non_interface <- nrow(assay_data[!get(site_type_col) %in% c("Binding interface site", "GTP binding interface site") & type == "surface"])
      
      if(background_core_non_interface > 0 && background_surface_non_interface > 0) {
        contingency_core_surface <- matrix(c(m_core, m_surface, 
                                             background_core_non_interface, background_surface_non_interface), nrow = 2)
        core_surface_fisher <- fisher.test(contingency_core_surface)
        core_surface_p <- core_surface_fisher$p.value
        core_surface_or <- core_surface_fisher$estimate
      } else {
        core_surface_p <- NA
        core_surface_or <- NA
      }
    } else {
      core_surface_p <- NA
      core_surface_or <- NA
    }
    
    # Save statistical results
    assay_stats[[assay]] <- list(
      n_interface = n_interface,
      m_non_interface = m_non_interface,
      m_core = m_core,
      m_surface = m_surface,
      n_p_value = n_p_value,
      m_p_value = m_p_value,
      core_surface_p = core_surface_p,
      core_surface_or = core_surface_or,
      threshold = threshold
    )
  }
  
  return(assay_stats)
}

# Convert p-values to stars
p_to_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

# Run statistical analysis
assay_stats <- count_high_effect_with_stats(data_plot_mutation_processed, all_assays, thresholds_list)

# Print detailed statistical results
cat("Detailed statistical results for all binders:\n")
cat("============================================\n\n")

for(assay in all_assays) {
  if(!is.null(assay_stats[[assay]])) {
    stats <- assay_stats[[assay]]
    cat("Binder:", assay, "\n")
    cat("Interface mutations (n):", stats$n_interface, "\n")
    cat("  p-value:", format(stats$n_p_value, scientific = TRUE, digits = 3), "\n")
    cat("  significance:", p_to_stars(stats$n_p_value), "\n")
    cat("Non-interface mutations (m):", stats$m_non_interface, "\n")
    cat("  p-value:", format(stats$m_p_value, scientific = TRUE, digits = 3), "\n")
    cat("  significance:", p_to_stars(stats$m_p_value), "\n")
    cat("  Core mutations:", stats$m_core, "\n")
    cat("  Surface mutations:", stats$m_surface, "\n")
    cat("  Core/Surface Fisher p-value:", format(stats$core_surface_p, scientific = TRUE, digits = 3), "\n")
    cat("  Core/Surface OR:", format(stats$core_surface_or, digits = 3), "\n")
    cat("  Threshold used:", format(stats$threshold, digits = 3), "\n")
    cat("----------------------------------------\n\n")
  }
}

# Combined heatmap plotting function
plot_combined_heatmap_final <- function(count_data, data_plot_mutation, assay_stats) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(tidyr)
  library(ggnewscale)
  library(patchwork)
  
  # Define order
  pos_order <- sort(unique(count_data$Pos_real))
  assay_order <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K13", "K19")
  assay_levels <- rev(assay_order)   # RAF1 at top
  
  # Ensure count_data has correct type information
  if(!"type" %in% names(count_data)) {
    type_info <- unique(data_plot_mutation[, .(Pos_real, type)])
    count_data <- merge(count_data, type_info, by = "Pos_real", all.x = TRUE)
  }
  
  # Melt data - only include cells with data
  count_cols <- paste0("count_", assay_order)
  melt_data <- melt(
    as.data.frame(count_data[, c("Pos_real", "type", count_cols), with = FALSE]),
    id.vars = c("Pos_real", "type"),
    variable.name = "Assay",
    value.name = "Count"
  )
  melt_data$Assay <- gsub("count_", "", melt_data$Assay)
  melt_data$Pos_real <- factor(melt_data$Pos_real, levels = pos_order)
  melt_data$Assay <- factor(melt_data$Assay, levels = assay_levels)
  
  # Filter: only keep Count > 0 and non-NA data
  melt_data_filtered <- melt_data[!is.na(melt_data$Count) & melt_data$Count > 0, ]
  
  # Core/surface background - extend height
  background_data <- unique(count_data[, .(Pos_real, type)])
  background_data <- background_data[order(Pos_real)]
  
  # Create extended background data
  extended_assays <- c("EXTENDED_TOP", assay_levels, "EXTENDED_BOTTOM")
  
  # Create complete background grid
  rects_data <- expand.grid(
    Pos_real = factor(pos_order, levels = pos_order),
    Assay = factor(extended_assays, levels = extended_assays)
  )
  # Merge type information
  rects_data <- merge(rects_data, background_data, by = "Pos_real", all.x = TRUE)
  
  # Binding interface sites
  annotation_list <- list()
  for(assay in assay_order) {
    site_type_col <- paste0("site_type_", assay)
    if(site_type_col %in% names(data_plot_mutation)) {
      site_info <- unique(data_plot_mutation[, .(Pos_real, site_type = get(site_type_col))])
      binding_sites <- site_info[site_type == "Binding interface site", Pos_real]
      if(length(binding_sites) > 0) {
        cat("Binding sites for", assay, ":", binding_sites, "\n")
        annotation_list[[assay]] <- data.frame(Pos_real = binding_sites, Assay = assay)
      }
    }
  }
  annotation_df <- if(length(annotation_list) > 0) do.call(rbind, annotation_list) else NULL
  
  # Debug: check annotation_df
  if(!is.null(annotation_df)) {
    cat("Annotation data frame created with", nrow(annotation_df), "rows\n")
    print(head(annotation_df))
    annotation_df$Pos_real <- factor(annotation_df$Pos_real, levels = pos_order)
    annotation_df$Assay <- factor(annotation_df$Assay, levels = assay_levels)
  } else {
    cat("No binding interface sites found!\n")
    # Check column names
    cat("Available site_type columns:\n")
    print(grep("site_type_", names(data_plot_mutation), value = TRUE))
  }
  
  # Binder annotation information
  annotation_labels <- data.frame(
    Assay = assay_order,
    Label = sapply(assay_order, function(a) {
      if(!is.null(assay_stats[[a]])) {
        stats <- assay_stats[[a]]
        n_star <- p_to_stars(stats$n_p_value)
        m_star <- p_to_stars(stats$m_p_value)
        paste0("n=", stats$n_interface, n_star,
               " m=", stats$m_non_interface, m_star,
               "\n→ m(OR)=", ifelse(!is.na(stats$core_surface_or),
                                    sprintf("%.2f", stats$core_surface_or), "NA"))
      } else "n=NA m=NA\n→ m(OR)=NA"
    })
  )
  
  # Main heatmap
  p_main <- ggplot() +
    # Background core/surface - use extended background
    geom_raster(
      data = rects_data,
      aes(x = Pos_real, y = Assay, fill = type),
      alpha = 0.5
    ) +
    scale_fill_manual(
      values = c("core" = "#3A7D44", "surface" = "#DDF6D2"),
      name = "Residue Type",
      na.value = "white"
    ) +
    # Correctly use new_scale_fill()
    new_scale_fill() +
    # Heatmap cells - only plot actual data regions
    geom_raster(
      data = melt_data_filtered,
      aes(x = Pos_real, y = Assay, fill = Count)
    ) +
    scale_fill_gradientn(
      colors = c("#FFF5F5", "#FFD8D8", "#F4270C"),
      limits = c(0, 20),
      breaks = seq(0, 20, 5),
      name = "Number of high-effect mutations"
    ) +
    # Binding interface sites - use more prominent markers
    {if(!is.null(annotation_df) && nrow(annotation_df) > 0)
      geom_point(
        data = annotation_df,
        aes(x = Pos_real, y = Assay),
        shape = 21, 
        fill = "#0C226F", 
        color = "#0C226F",
        size = 1.5,
        stroke = 1,
        alpha = 0.8
      )
    } +
    # Binder statistical information on right side
    annotate(
      "text",
      x = length(pos_order) + 2,
      y = seq_along(assay_levels) + 1,
      label = annotation_labels$Label,
      hjust = 0, size = 2.6, lineheight = 0.8
    ) +
    labs(x = "Position", y = NULL) +
    # Hide extended labels, only show actual assays
    scale_y_discrete(
      limits = extended_assays,
      labels = function(x) ifelse(x %in% c("EXTENDED_TOP", "EXTENDED_BOTTOM"), "", x)
    ) +
    # Control display range and proportion
    coord_fixed(
      ratio = 3,
      xlim = c(1, length(pos_order)),
      ylim = c(1.5, length(assay_levels) + 0.5),
      clip = "off",
      expand = FALSE
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 9),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(20, 100, 40, 10),
      legend.position = "bottom",
      panel.border = element_blank(),
      panel.spacing = unit(0, "lines")
    ) +
    # Add explanatory text
    annotate("text", x = length(pos_order) + 2, y = length(assay_levels) + 2,
             label = "n: interface, m: non-interface\n→: core/surface OR",
             hjust = 0, size = 2.2, color = "gray40", lineheight = 0.8) +
    # Add binding interface legend
    annotate("point", x = length(pos_order) + 2, y = length(assay_levels) + 3,
             shape = 21, fill = "#0C226F", color = "#0C226F", size = 1.5, stroke = 1) +
    annotate("text", x = length(pos_order) + 2.5, y = length(assay_levels) + 3,
             label = "Binding interface site", hjust = 0, size = 2.2, vjust = 0.5)
  
  # Use patchwork to combine
  p_combined <- p_main + 
    plot_annotation(
      title = "High-effect mutations with core/surface background",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", 
                                  margin = margin(b = 20)),
        plot.margin = margin(20, 20, 40, 20)
      )
    ) &
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = 20))
  
  return(p_combined)
}

# Generate combined heatmap
combined_heatmap <- plot_combined_heatmap_final(heatmap_data, data_plot_mutation_processed, assay_stats)
print(combined_heatmap)

# Save plot
ggsave("combined_mutation_heatmap_final.pdf", 
       combined_heatmap, device = cairo_pdf, width = 18, height = 9)
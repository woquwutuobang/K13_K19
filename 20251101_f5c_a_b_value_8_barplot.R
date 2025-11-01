# Energy-Distance Decay Analysis with Exponential Fitting
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===============================
# Function: Calculate Exponential Fit Parameters Only
# ===============================
calculate_exp_fit_parameters <- function(input, assay_sele, anno_file) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  # Define secondary structures
  rects_sheet <- data.frame(xstart = c(3, 38, 51, 77, 109, 139),
                            xend = c(9, 44, 57, 84, 115, 143),
                            col = c("β1", "β2", "β3", "β4", "β5", "β6"))
  
  rects_helix <- data.frame(xstart = c(15, 67, 87, 127, 148),
                            xend = c(24, 73, 104, 136, 166),
                            col = c("α1", "α2", "α3", "α4", "α5"))
  
  # Prepare results data frame
  results <- data.frame()
  
  # Function to calculate exponential fit
  calculate_exp_fit <- function(xvector, yvector, region_name) {
    # Remove NA values
    valid_idx <- complete.cases(data.frame(x = xvector, y = yvector))
    xvector <- xvector[valid_idx]
    yvector <- abs(yvector[valid_idx])  # Take absolute value of energy
    
    if (length(xvector) < 5) {  # Need enough points for fitting
      return(data.frame(
        region = region_name,
        a_value = NA,
        b_value = NA,
        p_value = NA,
        n_points = length(xvector),
        status = "insufficient_data"
      ))
    }
    
    # Initial parameter optimization
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess <- c(a = 1, b = -0.1)
    opt_params <- tryCatch(
      optim(initial_guess, residual_sum_of_squares, x = xvector, y = yvector)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    # Fit exponential model
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), 
          data = data.frame(x = xvector, y = yvector), 
          start = list(a = opt_params[1], b = opt_params[2]),
          control = nls.control(warnOnly = TRUE, maxiter = 100)),
      error = function(e) NULL
    )
    
    if (!is.null(fit_model) && !is.na(coef(fit_model)["b"])) {
      fit_summary <- summary(fit_model)
      coefs <- fit_summary$coefficients
      
      return(data.frame(
        region = region_name,
        a_value = coefs["a", "Estimate"],
        b_value = coefs["b", "Estimate"],
        p_value = coefs["b", "Pr(>|t|)"],
        n_points = length(xvector),
        status = "success"
      ))
    } else {
      return(data.frame(
        region = region_name,
        a_value = NA,
        b_value = NA,
        p_value = NA,
        n_points = length(xvector),
        status = "fit_failed"
      ))
    }
  }
  
  # Get all distance columns for different assays
  distance_assays <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  distance_cols <- paste0("scHAmin_ligand_", distance_assays)
  
  # Calculate for each distance assay
  for (dist_assay in distance_assays) {
    x_col <- paste0("scHAmin_ligand_", dist_assay)
    
    if (!x_col %in% colnames(anno_final)) next
    
    # Prepare main data
    xvector <- anno_final[[x_col]]
    yvector <- anno_final[[y_col]]
    pos_real <- anno_final$Pos_real
    
    main_df <- data.frame(
      x = xvector,
      y = yvector,
      pos = pos_real
    )
    main_df <- main_df[complete.cases(main_df), ]
    
    if (nrow(main_df) < 5) next
    
    # Calculate for overall data
    overall_result <- calculate_exp_fit(main_df$x, main_df$y, "overall")
    overall_result$assay <- assay_sele
    overall_result$distance_assay <- dist_assay
    results <- rbind(results, overall_result)
    
    # Calculate for sheet regions
    sheet_positions <- unique(unlist(lapply(1:nrow(rects_sheet), function(i) {
      rects_sheet$xstart[i]:rects_sheet$xend[i]
    })))
    
    sheet_df <- main_df[main_df$pos %in% sheet_positions, ]
    if (nrow(sheet_df) >= 5) {
      sheet_result <- calculate_exp_fit(sheet_df$x, sheet_df$y, "sheet")
      sheet_result$assay <- assay_sele
      sheet_result$distance_assay <- dist_assay
      results <- rbind(results, sheet_result)
    }
    
    # Calculate for helix regions
    helix_positions <- unique(unlist(lapply(1:nrow(rects_helix), function(i) {
      rects_helix$xstart[i]:rects_helix$xend[i]
    })))
    
    helix_df <- main_df[main_df$pos %in% helix_positions, ]
    if (nrow(helix_df) >= 5) {
      helix_result <- calculate_exp_fit(helix_df$x, helix_df$y, "helix")
      helix_result$assay <- assay_sele
      helix_result$distance_assay <- dist_assay
      results <- rbind(results, helix_result)
    }
    
    # Calculate for other regions (neither sheet nor helix)
    other_positions <- setdiff(main_df$pos, c(sheet_positions, helix_positions))
    other_df <- main_df[main_df$pos %in% other_positions, ]
    if (nrow(other_df) >= 5) {
      other_result <- calculate_exp_fit(other_df$x, other_df$y, "other")
      other_result$assay <- assay_sele
      other_result$distance_assay <- dist_assay
      results <- rbind(results, other_result)
    }
  }
  
  return(results)
}

# ===============================
# Main Analysis for All Assays
# ===============================
analyze_all_assays_cross <- function(input_template, anno_file, output_file = NULL) {
  
  assays <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  
  all_results <- data.frame()
  
  for (assay in assays) {
    input_file <- gsub("ASSAY", assay, input_template)
    cat("Processing", assay, "...\n")
    
    result <- calculate_exp_fit_parameters(
      input = input_file,
      assay_sele = assay,
      anno_file = anno_file
    )
    
    # Add group information
    result$group <- ifelse(assay %in% c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27"), "BI1", "BI2")
    
    all_results <- rbind(all_results, result)
  }
  
  # Save results if output file specified
  if (!is.null(output_file)) {
    fwrite(all_results, output_file)
    cat("Results saved to:", output_file, "\n")
  }
  
  return(all_results)
}

# ===============================
# Helper Functions
# ===============================

# Helper function to convert p-value to significance stars
get_significance_stars <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# ===============================
# CODE 1: Overall Structure Only (One Plot)
# ===============================

# Function to create overall structure plots
create_overall_plot <- function(results_data, save_path = NULL) {
  
  # Filter data - ONLY KEEP DIAGONAL (assay == distance_assay) and overall region
  plot_data <- results_data %>%
    filter(region == "overall" & !is.na(a_value) & !is.na(b_value) & assay == distance_assay) %>%
    mutate(
      assay = factor(assay, levels = c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")),
      group = factor(group, levels = c("BI1", "BI2")),
      # Add significance stars
      significance_a = sapply(p_value, get_significance_stars),
      significance_b = sapply(p_value, get_significance_stars)
    )
  
  if (nrow(plot_data) == 0) {
    cat("No data available for overall structure\n")
    return(NULL)
  }
  
  # Calculate y positions for significance annotations (above bars)
  plot_data <- plot_data %>%
    mutate(
      y_position_a = a_value * 0.8,
      # 让b值的星号更往下方，使用更大的乘数
      y_position_b = ifelse(b_value >= 0, b_value * 0.8, b_value * 0.8)
    )
  
  # 统一y轴范围
  a_y_limits <- c(0, 1.8)
  b_y_limits <- c(-0.15, 0.05)
  
  # Create a-value plot with colored stars
  p_a <- ggplot(plot_data, aes(x = assay, y = a_value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_a, label = significance_a), 
              position = position_dodge(0.8), 
              vjust = -0.2, size = 4, fontface = "bold", 
              color = "#F4AD0C", show.legend = FALSE) +  
    scale_fill_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    scale_y_continuous(limits = a_y_limits) +  # 统一a值y轴范围
    labs(
      title = "Initial energy (a) Values - Overall Structure\n(Own Assay Energy vs Own Assay Distance)",
      x = "Assay",
      y = "Initial energy (a) value",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  # Create b-value plot with colored stars
  p_b <- ggplot(plot_data, aes(x = assay, y = b_value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_b, label = significance_b), 
              position = position_dodge(0.8), 
              vjust = 1.0,  # 统一使用较大的vjust值让星号更靠下
              size = 4, fontface = "bold", 
              color = "#F4AD0C", show.legend = FALSE) + 
    scale_fill_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    scale_y_continuous(limits = b_y_limits) +  # 统一b值y轴范围
    labs(
      title = "Decay Rate (b) Values - Overall Structure\n(Own Assay Energy vs Own Assay Distance)",
      x = "Assay",
      y = "Decay Rate (b) value",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  # Combine plots
  final_plot <- p_a + p_b + plot_layout(ncol = 2)
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, width = 14, height = 6, device = cairo_pdf)
    cat("Overall structure plot saved to:", save_path, "\n")
  }
  
  return(final_plot)
}

# ===============================
# CODE 2: Secondary Structure Regions (One Plot)
# ===============================

# Function to create secondary structure plots
create_secondary_structure_plot <- function(results_data, save_path = NULL) {
  
  # Define regions and their titles - 按β/α/other顺序
  regions <- list(
    sheet = "β-Sheet Regions", 
    helix = "α-Helix Regions",
    other = "Other Regions"
  )
  
  # Prepare data for all secondary structure regions
  plot_data_list <- list()
  
  for (region_name in names(regions)) {
    region_data <- results_data %>%
      filter(region == region_name & !is.na(a_value) & !is.na(b_value) & assay == distance_assay) %>%
      mutate(
        assay = factor(assay, levels = c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")),
        group = factor(group, levels = c("BI1", "BI2")),
        region_label = factor(regions[[region_name]], 
                              levels = c("β-Sheet Regions", "α-Helix Regions", "Other Regions")),  # 设置因子水平确保顺序
        # Add significance stars
        significance_a = sapply(p_value, get_significance_stars),
        significance_b = sapply(p_value, get_significance_stars)
      )
    
    if (nrow(region_data) > 0) {
      plot_data_list[[region_name]] <- region_data
    }
  }
  
  if (length(plot_data_list) == 0) {
    cat("No data available for secondary structure regions\n")
    return(NULL)
  }
  
  # Combine all region data
  all_plot_data <- bind_rows(plot_data_list)
  
  # 统一y轴范围
  a_y_limits <- c(0, 1.8)
  b_y_limits <- c(-0.15, 0.05)
  
  # Calculate y positions for significance annotations
  all_plot_data <- all_plot_data %>%
    group_by(region) %>%
    mutate(
      y_position_a = a_value * 0.8,
      # 让b值的星号更往下方
      y_position_b = ifelse(b_value >= 0, b_value * 0.8, b_value * 0.8)
    ) %>%
    ungroup()
  
  # Create a-value plot for all secondary structures
  p_a <- ggplot(all_plot_data, aes(x = assay, y = a_value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_a, label = significance_a), 
              position = position_dodge(0.8), 
              vjust = -0.2, size = 4, fontface = "bold", 
              color = "#F4AD0C", show.legend = FALSE) +  
    facet_wrap(~ region_label, ncol = 3, scales = "fixed") +  # 使用fixed确保统一y轴范围
    scale_fill_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    scale_y_continuous(limits = a_y_limits) +  # 统一a值y轴范围
    labs(
      title = "Initial energy (a) Values by Secondary Structure\n(Own Assay Energy vs Own Assay Distance)",
      x = "Assay",
      y = "Initial energy (a) value",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 10)
    )
  
  # Create b-value plot for all secondary structures
  p_b <- ggplot(all_plot_data, aes(x = assay, y = b_value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_b, label = significance_b), 
              position = position_dodge(0.8), 
              vjust = 1.0,  # 统一使用较大的vjust值让星号更靠下
              size = 4, fontface = "bold", 
              color = "#F4AD0C", show.legend = FALSE) +  
    facet_wrap(~ region_label, ncol = 3, scales = "fixed") +  # 使用fixed确保统一y轴范围
    scale_fill_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    scale_y_continuous(limits = b_y_limits) +  # 统一b值y轴范围
    labs(
      title = "Decay Rate (b) Values by Secondary Structure\n(Own Assay Energy vs Own Assay Distance)",
      x = "Assay",
      y = "Decay Rate (b) value",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 10)
    )
  
  # Combine plots
  final_plot <- p_a / p_b + plot_layout(heights = c(1, 1))
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, width = 14, height = 10, device = cairo_pdf)
    cat("Secondary structure plot saved to:", save_path, "\n")
  }
  
  return(final_plot)
}


# ===============================
# Summary Statistics Function
# ===============================
print_summary_statistics <- function(results_data) {
  cat("=== Summary Statistics (Diagonal Only) ===\n")
  
  # Filter for diagonal only
  diagonal_data <- results_data %>% filter(assay == distance_assay)
  
  summary_stats <- diagonal_data %>%
    group_by(region, group) %>%
    summarise(
      n_success_a = sum(!is.na(a_value)),
      n_success_b = sum(!is.na(b_value)),
      mean_a = mean(a_value, na.rm = TRUE),
      mean_b = mean(b_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  cat("\n=== Success Rates by Region (Diagonal Only) ===\n")
  success_rates <- diagonal_data %>%
    group_by(region) %>%
    summarise(
      success_rate_a = mean(!is.na(a_value)) * 100,
      success_rate_b = mean(!is.na(b_value)) * 100,
      .groups = "drop"
    )
  
  print(success_rates)
  
  cat("\n=== P-value Summary (Diagonal Only) ===\n")
  pvalue_summary <- diagonal_data %>%
    filter(!is.na(p_value)) %>%
    group_by(region) %>%
    summarise(
      n_significant = sum(p_value < 0.05),
      n_total = n(),
      significance_rate = n_significant / n_total * 100,
      .groups = "drop"
    )
  
  print(pvalue_summary)
}

# ===============================
# Execute Analysis
# ===============================

# Run analysis
results <- analyze_all_assays_cross(
  input_template = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_ASSAY.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",  
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f5/20251101/exponential_fit_results_diagonal_all_info.csv"
)

# CODE 1: Create overall structure plot
overall_plot <- create_overall_plot(
  results, 
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f5/20251101/overall_structure_8barplot3.pdf"
)

# CODE 2: Create secondary structure plot
secondary_plot <- create_secondary_structure_plot(
  results, 
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f5/20251101/secondary_structure_8barplot3.pdf"
)

# Print summary statistics
#print_summary_statistics(results)
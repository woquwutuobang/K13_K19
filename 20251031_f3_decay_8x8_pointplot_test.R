# Energy-Distance Decay Analysis with Exponential Fitting
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

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
  distance_assays <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")
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
  
  assays <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")
  
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
    result$group <- ifelse(assay %in% c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG"), "BI1", "BI2")
    
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
# New Scatter Plot Functions
# ===============================

# Function to create 8x8 scatter plot matrix
create_scatter_matrix <- function(results_data, region_filter = "overall", 
                                  parameter = c("a_value", "b_value"), 
                                  title_suffix = "") {
  
  # Filter data
  plot_data <- results_data %>%
    filter(region == region_filter & !is.na(.data[[parameter]])) %>%
    mutate(
      assay = factor(assay, levels = c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")),
      distance_assay = factor(distance_assay, levels = c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")),
      group = factor(group, levels = c("BI1", "BI2")),
      # Create significance labels
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 6) +
             theme_void())
  }
  
  # Create the scatter plot matrix
  p <- ggplot(plot_data, aes(x = assay, y = distance_assay)) +
    geom_point(aes(size = abs(.data[[parameter]]), color = group, fill = group), 
               shape = 21, alpha = 0.8) +
    # Add parameter values as text
    geom_text(aes(label = sprintf("%.3f", .data[[parameter]])), 
              size = 3.5, fontface = "bold", vjust = -1.5) +
    # Add significance stars
    geom_text(aes(label = significance), 
              size = 4, fontface = "bold", vjust = 1.5, color = "red") +
    scale_size_continuous(
      name = paste0("|", parameter, "|"),
      range = c(2, 8),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_color_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    scale_fill_manual(values = c("BI1" = "#1B38A6", "BI2" = "#F4270C")) +
    labs(
      title = paste0(toupper(parameter), " Values - ", region_filter, title_suffix),
      x = "Assay (Energy Data)",
      y = "Distance Assay",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    ) +
    # Ensure square aspect ratio for matrix appearance
    coord_fixed(ratio = 1)
  
  return(p)
}

# Master plotting function for scatter matrices
plot_scatter_matrices <- function(results_data, save_path = NULL) {
  
  # Define regions and their titles
  regions <- list(
    overall = "Overall Structure",
    sheet = "β-Sheet Regions", 
    helix = "α-Helix Regions",
    other = "Other Regions"
  )
  
  plots_list <- list()
  
  # Create scatter matrices for each region and parameter
  for (region_name in names(regions)) {
    for (parameter in c("a_value", "b_value")) {
      
      plot_title <- paste0(regions[[region_name]], " - ", 
                           ifelse(parameter == "a_value", "Amplitude (a)", "Decay Rate (b)"))
      
      p <- create_scatter_matrix(
        results_data, 
        region_filter = region_name,
        parameter = parameter,
        title_suffix = paste0(" - ", regions[[region_name]])
      )
      
      plots_list[[paste0(region_name, "_", parameter)]] <- p
    }
  }
  
  # Combine all plots in a 2x2 grid (one for each region, with a and b as subplots)
  final_plot <- wrap_plots(plots_list, ncol = 2, byrow = TRUE) + 
    plot_annotation(
      title = "Exponential Fit Parameters: 8x8 Scatter Matrix",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    ) &
    theme(plot.tag = element_text(size = 12, face = "bold"))
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, width = 20, height = 24, device = cairo_pdf)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(final_plot)
}

# ===============================
# Summary Statistics Function
# ===============================
print_summary_statistics <- function(results_data) {
  cat("=== Summary Statistics ===\n")
  
  summary_stats <- results_data %>%
    group_by(region, group) %>%
    summarise(
      n_success_a = sum(!is.na(a_value)),
      n_success_b = sum(!is.na(b_value)),
      mean_a = mean(a_value, na.rm = TRUE),
      mean_b = mean(b_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  cat("\n=== Success Rates by Region ===\n")
  success_rates <- results_data %>%
    group_by(region) %>%
    summarise(
      success_rate_a = mean(!is.na(a_value)) * 100,
      success_rate_b = mean(!is.na(b_value)) * 100,
      .groups = "drop"
    )
  
  print(success_rates)
}

# ===============================
# Execute Analysis
# ===============================
results <- analyze_all_assays_cross(
  input_template = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_ASSAY.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",  
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251030/exponential_fit_results_scatter.csv"
)

# 创建散点图矩阵
scatter_plots <- plot_scatter_matrices(
  results, 
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure3/20251030/scatter_matrix_comparison.pdf"
)

# 打印统计摘要
print_summary_statistics(results)

# 单独显示每个区域的图形（可选）
# 如果您想单独查看某个区域的图形，可以使用以下代码：
cat("\n=== Generating Individual Plots ===\n")
for (region in c("overall", "sheet", "helix", "other")) {
  for (param in c("a_value", "b_value")) {
    individual_plot <- create_scatter_matrix(results, region, param)
    print(individual_plot)
    cat("Generated plot for", region, "-", param, "\n")
  }
}
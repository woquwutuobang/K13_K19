# Energy-Distance Decay Analysis with Exponential Fitting
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
# ===============================
# Function: Calculate and Plot Exponential Fit
# ===============================
calculate_and_plot_exp_fit <- function(input, assay_sele, anno_file) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  # Distance columns for different assays
  distance_assays <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")
  
  plots_list <- list()
  
  # Loop through all distance assays
  for (dist_assay in distance_assays) {
    x_col <- paste0("scHAmin_ligand_", dist_assay)
    
    if (!x_col %in% colnames(anno_final)) next
    
    # Prepare main data
    xvector <- anno_final[[x_col]]
    yvector <- anno_final[[y_col]]
    pos_real <- anno_final$Pos_real
    
    main_df <- data.frame(
      x = xvector,
      y = abs(yvector),
      pos = pos_real
    )
    main_df <- main_df[complete.cases(main_df), ]
    
    if (nrow(main_df) < 5) next
    
    # Fit exponential model
    fit_model <- tryCatch(
      {
        residual_sum_of_squares <- function(params, x, y) {
          a <- params[1]; b <- params[2]
          predicted <- a * exp(b * x)
          sum((y - predicted)^2, na.rm = TRUE)
        }
        
        initial_guess <- c(a = 1, b = -0.1)
        opt_params <- optim(initial_guess, residual_sum_of_squares, x = main_df$x, y = main_df$y)$par
        
        nls(y ~ a * exp(b * x), 
            data = main_df, 
            start = list(a = opt_params[1], b = opt_params[2]),
            control = nls.control(warnOnly = TRUE, maxiter = 100))
      },
      error = function(e) NULL
    )
    
    # Base plot
    p <- ggplot(main_df, aes(x = x, y = y)) +
      geom_point(alpha = 0.2, size = 1.5, color = "#75C2F6") +
      labs(
        title = paste0(assay_sele, " vs ", dist_assay),
        x = paste0("Distance to ", dist_assay, " (Å)"),
        y = "|ΔΔG| (kcal/mol)",
        subtitle = ifelse(!is.null(fit_model), 
                          paste0("y = ", sprintf("%.3f", coef(fit_model)["a"]), 
                                 " * exp(", sprintf("%.3f", coef(fit_model)["b"]), " * x)"),
                          "Fit failed")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        panel.grid.minor = element_blank()
      ) +
      coord_cartesian(xlim = c(0, 35), ylim = c(0, 2.5))   
    
    # Add fitted curve and annotation
    if (!is.null(fit_model)) {
      x_range <- seq(0, 35, length.out = 100)
      fit_data <- data.frame(x = x_range,
                             y = predict(fit_model, newdata = data.frame(x = x_range)))
      p <- p + geom_line(data = fit_data, aes(x = x, y = y), 
                         color = "#F4270C", linewidth = 1, alpha = 0.8)
      
      fit_summary <- summary(fit_model)
      p_value <- fit_summary$coefficients["b", "Pr(>|t|)"]
      
      p <- p + annotate("text", 
                        x = 20, y = 2.0,   
                        label = paste0("a = ", sprintf("%.3f", coef(fit_model)["a"]),
                                       "\nb = ", sprintf("%.3f", coef(fit_model)["b"]),
                                       "\np = ", sprintf("%.3f", p_value)),
                        size = 3.8,  
                        hjust = 0, 
                        color = "#A31300")
    }
    
    plots_list[[paste0(assay_sele, "_", dist_assay)]] <- p
  }
  
  return(plots_list)
}

# ===============================
# Main Analysis and Plotting Function
# ===============================
plot_all_decay_curves <- function(input_template, anno_file, save_path = NULL) {
  
  assays <- c("RAF1", "SOS1", "K55", "K27", "RALGDS", "PI3KCG", "K19", "K13")
  
  all_plots <- list()
  
  for (assay in assays) {
    input_file <- gsub("ASSAY", assay, input_template)
    cat("Processing", assay, "...\n")
    
    assay_plots <- calculate_and_plot_exp_fit(
      input = input_file,
      assay_sele = assay,
      anno_file = anno_file
    )
    
    all_plots <- c(all_plots, assay_plots)
  }
  
  # Combine all plots
  final_plot <- wrap_plots(all_plots, ncol = 8, nrow = 8) + 
    plot_annotation(
      title = "Energy-Distance Decay Analysis: Exponential Fitting",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, width = 20, height = 20, device = cairo_pdf)
    cat("Plot saved to:", save_path, "\n")
  }
  
  return(final_plot)
}

# ===============================
# Execute
# ===============================
decay_plots <- plot_all_decay_curves(
  input_template = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_ASSAY.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",  
  save_path = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f5/20251101/energy_distance_decay_curves_8x8 scatterplot_9.pdf"
)

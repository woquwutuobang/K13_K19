# Energy-Distance Decay Analysis with Exponential Fitting
# This script analyzes the relationship between binding energy changes and distance from binding partners

library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# ===============================
# General Function: Distance Effect Fitting (Exponential Model)
# ===============================
plot_energy_distance_decay_expfit <- function(input, assay_sele, anno_file,
                                              x_cols = c("scHAmin_ligand_K13", "scHAmin_ligand_RAF1"),
                                              titles = c("Distance to K13", "Distance to RAF1"),
                                              x_range = c(0, 35), y_range = c(0, 3),
                                              save_path = NULL) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, c(3, 20:22)]
  colnames(data)[2:4] <- paste0(colnames(data)[2:4], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno_final <- merge(anno, data, by = "Pos", all = TRUE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  plot_list <- list()
  
  # 2. Iterate through all x_cols
  for (i in seq_along(x_cols)) {
    xvector <- anno_final[[x_cols[i]]]
    yvector <- abs(anno_final[[y_col]])  # Take absolute value of energy
    
    df <- data.frame(x = xvector, y = yvector)
    df <- df[complete.cases(df), ]
    
    # Initial parameter optimization
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    initial_guess <- c(a = 1, b = -0.1)
    opt_params <- tryCatch(
      optim(initial_guess, residual_sum_of_squares, x = df$x, y = df$y)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    # Fit exponential model
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), data = df, start = list(a = opt_params[1], b = opt_params[2])),
      error = function(e) NULL
    )
    
    # Fit results and annotation
    fit_df <- data.frame()
    annotation_text <- NULL
    
    if (!is.null(fit_model)) {
      x_seq <- seq(min(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE), length.out = 200)
      y_fit <- predict(fit_model, newdata = data.frame(x = x_seq))
      fit_df <- data.frame(x = x_seq, y = y_fit)
      
      fit_summary <- summary(fit_model)
      coefs <- fit_summary$coefficients
      a_val <- round(coefs["a", "Estimate"], 3)
      b_val <- round(coefs["b", "Estimate"], 3)
      p_val_b <- coefs["b", "Pr(>|t|)"]
      
      p_text <- if (is.na(p_val_b)) {
        "p = NA"
      } else if (p_val_b < 0.001) {
        "p < 0.001"
      } else if (p_val_b < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b, 3))
      }
      
      annotation_text <- paste0("a = ", a_val, ", b = ", b_val, "\n", p_text)
    }
    
    # Calculate density data
    density_data <- density(df$x, na.rm = TRUE)
    density_df <- data.frame(
      x = density_data$x,
      density = density_data$y
    )
    
    # Main plot
    p_main <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.1, size = 0.7, color = "#48B3AF") +
      coord_cartesian(xlim = x_range, ylim = y_range) +
      theme_classic(base_size = 13) +
      labs(
        x = titles[i],
        y = paste0("ΔΔG (", assay_sele, " binding)")
      )
    
    if (nrow(fit_df) > 0) {
      p_main <- p_main + geom_line(data = fit_df, aes(x = x, y = y),
                                   inherit.aes = FALSE, color = "#9A3F3F", linewidth = 0.8)
    }
    if (!is.null(annotation_text)) {
      p_main <- p_main + annotate("text",
                                  x = Inf, y = Inf,
                                  hjust = 1.1, vjust = 1.5,
                                  label = annotation_text,
                                  size = 4, color = "black")
    }
    
    # Density plot - Add axis information
    p_density <- ggplot(density_df, aes(x = x, y = density)) +
      geom_area(fill = "#48B3AF", alpha = 0.4) +
      coord_cartesian(xlim = x_range) +
      theme_classic(base_size = 10) +  # Use theme_classic instead of theme_void
      theme(
        plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5),  # Adjust margins
        axis.title.x = element_blank(),  # Remove x-axis title (main plot already has it)
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.title.y = element_text(size = 8),  # Reduce y-axis title font size
        axis.text.y = element_text(size = 7),   # Reduce y-axis tick label size
        panel.grid = element_blank()     # Ensure no grid lines
      ) +
      labs(
        x = NULL,
        y = "Density"
      )
    
    # Combine main plot and density plot
    p_combined <- p_main / p_density + 
      plot_layout(heights = c(3, 1))
    
    plot_list[[i]] <- p_combined
  }
  
  # 3. Combine plots
  final_plot <- wrap_plots(plot_list, ncol = 2)
  
  # 4. Save
  if (!is.null(save_path)) {
    ggsave(save_path, final_plot, device = cairo_pdf, height = 6, width = 8)
  }
  
  return(final_plot)
}

# ===============================
# Analysis for Different Binding Partners
# ===============================

### RAF1
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/RAF1_expfit.pdf"
)
print(plot_result)

### SOS1
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_SOS.txt",
  assay_sele = "SOS1",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/SOS1_expfit.pdf"
)
print(plot_result)

### K55
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_K55.txt",
  assay_sele = "K55",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/K55_expfit.pdf"
)
print(plot_result)

### K27
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_K27.txt",
  assay_sele = "K27",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/K27_expfit.pdf"
)
print(plot_result)

### RALGDS
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_RAL.txt",
  assay_sele = "RALGDS",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/RALGDS_expfit.pdf"
)
print(plot_result)

### PIK3CG
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_PI3.txt",
  assay_sele = "PIK3CG",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/PIK3CG_expfit.pdf"
)
print(plot_result)

### K13
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/K13_expfit.pdf"
)
print(plot_result)

### K19
plot_result <- plot_energy_distance_decay_expfit(
  input = "path/to/weights_Binding_K19.txt",
  assay_sele = "K19",
  anno_file = "path/to/annotations.csv",
  save_path = "path/to/K19_expfit.pdf"
)
print(plot_result)
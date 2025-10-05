# KRAS DARPins Binding Correlation Analysis
# This script analyzes the correlation between mutation effects on K13 and K19 binding

library(data.table)
library(krasddpcams)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

# Read ddG data for K13 and K19 binding
ddG_k13 <- krasddpcams__read_ddG(
  ddG = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13"
)

ddG_k19 <- krasddpcams__read_ddG(
  ddG = "path/to/weights_Binding_K19.txt", 
  assay_sele = "K19"
)

# Define binding interface residues (same for both DARPins)
binding_interface_residues <- c(98, 107, 101, 102, 99, 136, 68, 95, 137, 94, 133, 90, 129, 87, 91, 88)

# Filter for interface residues and select relevant columns
ddG_k13_interface <- ddG_k13[Pos %in% binding_interface_residues, c(1:3, 23:27)]
ddG_k19_interface <- ddG_k19[Pos %in% binding_interface_residues, c(1:3, 23:27)]

# Reshape data to wide format
ddG_k13_wide <- ddG_k13_interface %>%
  spread(key = assay, value = `mean_kcal/mol`)

ddG_k19_wide <- ddG_k19_interface %>%
  spread(key = assay, value = `mean_kcal/mol`)

# Merge K13 and K19 data
ddG_combined <- merge(ddG_k13_wide, ddG_k19_wide, by = c("mt", "Pos_real"), all = TRUE)

# Calculate correlation coefficient
correlation_coef <- cor(ddG_combined$K13, ddG_combined$K19, use = "complete.obs")

# Identify outliers based on linear regression residuals
fit <- lm(K19 ~ K13, data = ddG_combined)
ddG_combined$abs_residual <- abs(resid(fit))

# Mark top N outliers
num_outliers <- 10
ddG_combined <- ddG_combined %>%
  mutate(
    outlier_rank = rank(-abs_residual, ties.method = "first"),
    is_outlier = outlier_rank <= num_outliers
  )

# Get outliers for labeling
outliers_label <- ddG_combined %>%
  filter(is_outlier == TRUE) %>%
  arrange(desc(abs_residual))

# Create correlation plot with error bars
correlation_plot <- ggplot(ddG_combined, aes(x = K13, y = K19)) +
  # Regression line
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8, alpha = 0.5) +
  
  # Error bars for all points
  geom_errorbar(aes(ymin = K19 - `std_kcal/mol.y`, ymax = K19 + `std_kcal/mol.y`),
                width = 0.05, color = "grey70", alpha = 0.3) +
  geom_errorbarh(aes(xmin = K13 - `std_kcal/mol.x`, xmax = K13 + `std_kcal/mol.x`),
                 height = 0.05, color = "grey70", alpha = 0.3) +
  
  # All points
  geom_point(color = "#0C226F", alpha = 0.6, size = 2) +
  
  # Highlighted error bars for outliers
  geom_errorbar(data = outliers_label,
                aes(ymin = K19 - `std_kcal/mol.y`, ymax = K19 + `std_kcal/mol.y`),
                width = 0.08, color = "#F4270C", alpha = 0.6, size = 0.8) +
  geom_errorbarh(data = outliers_label,
                 aes(xmin = K13 - `std_kcal/mol.x`, xmax = K13 + `std_kcal/mol.x`),
                 height = 0.08, color = "#F4270C", alpha = 0.6, size = 0.8) +
  
  # Highlight outlier points
  geom_point(data = outliers_label, color = "#F4270C", size = 3) +
  
  # Label outliers
  geom_text_repel(
    data = outliers_label,
    aes(label = paste0(mt, " (Δ=", round(abs_residual, 2), ")")),
    color = "#F4270C",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.7,
    point.padding = 0.3,
    segment.color = "#F4270C",
    segment.alpha = 0.3
  ) +
  
  # Add statistics annotation
  annotate("text",
           x = min(ddG_combined$K13, na.rm = TRUE),
           y = max(ddG_combined$K19, na.rm = TRUE),
           label = paste(
             "Pearson r =", round(correlation_coef, 3), "\n",
             "Top", num_outliers, "deviating points labeled\n",
             "Error bars show ±1 SD"
           ),
           hjust = 0, vjust = 1, size = 4
  ) +
  
  labs(
    title = "K13 vs K19 Binding ΔΔG Correlation",
    x = "ΔΔG K13 Binding (kcal/mol)",
    y = "ΔΔG K19 Binding (kcal/mol)",
    caption = "Error bars represent standard deviation. Red points show largest deviations from linear relationship."
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.caption = element_text(size = 9, color = "grey50", hjust = 0)
  )

# Display plot
print(correlation_plot)

# Save plot
ggsave("K13_K19_binding_correlation.pdf",
       plot = correlation_plot,
       device = cairo_pdf,
       height = 6,
       width = 6,
       dpi = 300)

# Generate outlier analysis report
cat("=== Top", num_outliers, "Most Deviating Mutations ===\n")
outlier_analysis <- outliers_label %>%
  mutate(
    K13_error_range = 2 * `std_kcal/mol.x`,
    K19_error_range = 2 * `std_kcal/mol.y`,
    likely_measurement_noise = abs_residual < (K13_error_range + K19_error_range) / 2
  ) %>%
  select(mt, Pos_real, K13, `std_kcal/mol.x`, K19, `std_kcal/mol.y`, 
         abs_residual, K13_error_range, K19_error_range, likely_measurement_noise) %>%
  arrange(desc(abs_residual))

print(outlier_analysis)

# Summary statistics
noise_count <- sum(outlier_analysis$likely_measurement_noise, na.rm = TRUE)
cat("\n=== Summary ===\n")
cat("Correlation coefficient:", round(correlation_coef, 3), "\n")
cat("Outliers potentially due to measurement noise:", noise_count, "/", num_outliers, "\n")
cat("Outliers likely representing biological differences:", num_outliers - noise_count, "/", num_outliers, "\n")
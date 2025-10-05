library(krasddpcams)
library(data.table)
library(dplyr)
library(beeswarm)
library(wesanderson)
library(ggplot2)
library(scales)

# Function to define amino acid chemical types
define_chemotype <- function(aa) {
  aromatic <- c("F", "W", "Y")
  aliphatic <- c("A", "V", "I", "L", "M")
  polar_uncharged <- c("S", "T", "N", "Q", "C")
  positive <- c("K", "R", "H")
  negative <- c("D", "E")
  special <- c("G", "P")
  
  if (aa %in% aromatic) return("Aromatic")
  if (aa %in% aliphatic) return("Aliphatic")
  if (aa %in% polar_uncharged) return("Polar uncharged")
  if (aa %in% positive) return("Positive")
  if (aa %in% negative) return("Negative")
  if (aa %in% special) return("Special")
  return(NA)
}

# Main function to create beeswarm plots of ΔΔG values for binding interface residues
plot_ddG_beeswarm <- function(ddG_file, assay_sele, residues, output_file, 
                              width = 14, height = 8) {
  
  ## ===== 1. Read and process data =====
  ddG <- krasddpcams__read_ddG(ddG_file, assay_sele)
  ddG <- ddG[, c(1:3, 23, 26, 27)]
  
  # Add chemical type classification
  ddG$Chemotype <- sapply(ddG$mt_codon, define_chemotype)
  ddG$sites <- paste(ddG$wt_codon, ddG$Pos_real, sep = "")
  
  ## ===== 2. Filter for interface residues =====
  df_plot <- ddG[ddG$sites %in% residues, ]
  
  # Define color scheme for chemical types
  chemotype.cols <- c(
    "Aromatic"       = wes_palette("Darjeeling2", 6, type = "continuous")[1],
    "Aliphatic"      = wes_palette("Darjeeling2", 6, type = "continuous")[2],
    "Polar uncharged"= wes_palette("Darjeeling2", 6, type = "continuous")[3],
    "Positive"       = wes_palette("Darjeeling2", 6, type = "continuous")[4],
    "Negative"       = wes_palette("Darjeeling2", 6, type = "continuous")[5],
    "Special"        = wes_palette("Darjeeling2", 6, type = "continuous")[6]
  )
  
  ## ===== 3. Create beeswarm plot =====
  # Use cairo_pdf for high-quality output
  cairo_pdf(output_file, width = width, height = height)
  par(mar = c(8, 6, 2, 2))
  
  # Set residue order (reversed for plotting)
  residues_ordered <- rev(residues)
  
  # Prepare beeswarm data
  beeswarm.out <- split(df_plot$`mean_kcal/mol`, df_plot$sites)
  
  # Create empty boxplot to set up coordinate system
  boxplot(beeswarm.out, col = "white", border = "white", outline = F, horizontal = F,
          ylim = c(min(df_plot$`mean_kcal/mol`, na.rm = T), 
                   max(df_plot$`mean_kcal/mol`, na.rm = T)),
          frame = F, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  
  # Reference line at ΔΔG = 0
  abline(h = 0, col = "black", lty = "dotted", lwd = 2)
  
  # ===== Add mean ΔΔG bars for each site =====
  site_means <- tapply(df_plot$`mean_kcal/mol`, df_plot$sites, mean, na.rm = TRUE)
  site_means <- site_means[residues_ordered]  # Order by reversed residue list
  
  for (i in seq_along(residues_ordered)) {
    segments(x0 = i - 0.3, x1 = i + 0.3, 
             y0 = site_means[i],
             y1 = site_means[i], 
             col = alpha("grey40", 0.7), lwd = 4)
  }
  
  # ===== Create beeswarm plot =====
  beeswarm_result <- beeswarm(df_plot$`mean_kcal/mol` ~ factor(df_plot$sites, levels = residues_ordered),
                              method = "swarm", 
                              do.plot = FALSE)
  
  # Match data with chemical types
  plot_data <- data.frame(
    x = beeswarm_result$x,
    y = beeswarm_result$y,
    site = beeswarm_result$x.orig,
    value = beeswarm_result$y.orig
  )
  
  plot_data <- merge(plot_data, df_plot[, c("sites", "mean_kcal/mol", "mt_codon", "Chemotype")],
                     by.x = c("site", "value"), 
                     by.y = c("sites", "mean_kcal/mol"),
                     all.x = TRUE)
  
  # Assign colors based on chemical type
  plot_data$color <- chemotype.cols[plot_data$Chemotype]
  
  # Plot amino acid labels
  text(x = plot_data$x, 
       y = plot_data$y, 
       labels = plot_data$mt_codon, 
       col = plot_data$color, 
       cex = 2)
  
  # X-axis labels
  axis(1, at = 1:length(residues_ordered), labels = residues_ordered, las = 2, cex.axis = 1.2)
  
  # Y-axis label (ΔΔG)
  mtext(side = 2, line = 3, cex = 1.5,
        expression(Delta * Delta * "G (kcal/mol)"))
  
  # Y-axis ticks
  axis(2, las = 2, cex.axis = 1.2)
  
  # Legend for chemical types
  legend("bottom", horiz = TRUE,
         legend = names(chemotype.cols),
         pch = 21, pt.bg = chemotype.cols,
         bty = "n", cex = 1.2, inset = c(0, -0.2), xpd = TRUE)
  
  dev.off()
  
  cat("Plot saved to:", output_file, "\n")
  return(invisible())
}

# ===== Usage Examples =====

# Define binding interface residue list
residues_list <- c("K88", "E91", "T87", "Q129", "F90", "L133", "H94", 
                   "Y137", "H95", "R68", "S136", "Q99", "R102", "K101", 
                   "E107", "E98")

# Generate beeswarm plot for K13 binding interface
plot_ddG_beeswarm(
  ddG_file = "path/to/weights_Binding_K13.txt",
  assay_sele = "K13",
  residues = residues_list,
  output_file = "path/to/binding_interface_ddG_beeswarm_K13.pdf"
)

# Generate beeswarm plot for K19 binding interface
plot_ddG_beeswarm(
  ddG_file = "path/to/weights_Binding_K19.txt", 
  assay_sele = "K19",
  residues = residues_list,
  output_file = "path/to/binding_interface_ddG_beeswarm_K19.pdf"
)

# Summary of analysis
cat("\nBinding Interface ΔΔG Analysis Summary:\n")
cat("========================================\n")
cat("• Analyzed", length(residues_list), "binding interface residues\n")
cat("• Chemical types classified for all amino acid substitutions\n")
cat("• Beeswarm plots show distribution of ΔΔG values\n")
cat("• Gray bars indicate mean ΔΔG for each residue position\n")
cat("• Colors represent amino acid chemical properties\n")
cat("• Dotted line at ΔΔG = 0 indicates no change from wild-type\n")
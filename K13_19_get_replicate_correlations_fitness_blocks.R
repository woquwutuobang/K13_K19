# Fitness Correlation Analysis Across Replicates
# This script creates correlation plots to assess reproducibility between experimental replicates

library(data.table)
library(krasddpcams)
library(ggplot2)
library(GGally)
library(ggpubr)

# Define color scheme
colour_scheme <- list(
  "blue" = "#1B38A6",
  "red" = "#F4270C", 
  "orange" = "#F4AD0C",
  "green" = "#09B636",
  "yellow" = "#F1DD10",
  "purple" = "#6D17A0",
  "pink" = "#FFB0A5",
  "light orange" = "#FFE4A5",
  "light blue" = "#9DACE3",
  "light green" = "#97E9AD",
  "light red" = "#FF6A56",
  "dark red" = "#A31300",
  "dark blue" = "#0C226F",
  "dark green" = "#007A20"
)

#' Create fitness correlation plot across experimental blocks
#'
#' @param input Normalized fitness data
#' @param assay_name Name of the assay for plot title
#' @param colour_scheme Color scheme for different blocks
#' @return GGally plot object showing correlations between replicates

krasddpcams__plot_fitness_correlation_blocks <- function(
    input,
    assay_name,
    colour_scheme
){
  d <- GGally::ggpairs(
    input,
    columns = c(18,19,20), ##K13/K19
    #columns = c(18, 20, 22),##abundance/RAF1
    columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
    mapping = ggplot2::aes(color = block),
    lower = list(continuous = function(data, mapping, ...){
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_bin2d(bins = 100, alpha = 0.2) +
        ggplot2::scale_fill_gradient(low = "white", high = "black")
    }),
    upper = list(continuous = GGally::wrap("cor", 
                                           mapping = ggplot2::aes(color = block), 
                                           size = 8 * 0.35)),
    diag = list(continuous = "blankDiag")
  ) +
    ggplot2::scale_color_manual(values = c(
      colour_scheme[["red"]], 
      colour_scheme[["green"]], 
      colour_scheme[["blue"]]
    )) +
    ggplot2::ggtitle(assay_name) +
    ggpubr::theme_classic2() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 8),
      strip.text.x = ggplot2::element_text(size = 8),
      strip.text.y = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8)
    )
  
  return(d)
}


# =====================
# Load and Process Data
# =====================

# Load K19 binding data
K19_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/K19_block1.RData",
  block2_dimsum_df = "path/to/K19_block2.RData",
  block3_dimsum_df = "path/to/K19_block3.RData"
)

# Load K13 binding data
K13_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/K13_block1.RData",
  block2_dimsum_df = "path/to/K13_block2.RData",
  block3_dimsum_df = "path/to/K13_block3.RData"
)

# Load protein stability data
stability_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/abundance_block1.RData",
  block2_dimsum_df = "path/to/abundance_block2.RData",
  block3_dimsum_df = "path/to/abundance_block3.RData"
)

# Load RAF1 binding data
RAF1_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "path/to/RAF1_block1.RData",
  block2_dimsum_df = "path/to/RAF1_block2.RData",
  block3_dimsum_df = "path/to/RAF1_block3.RData"
)

# =====================
# Generate Correlation Plots
# =====================

# K13 correlation plot
d_K13 <- krasddpcams__plot_fitness_correlation_blocks(K13_nor_df, "K13", colour_scheme)
print(d_K13)
ggplot2::ggsave("path/to/K13_fitness_correlation.pdf", d_K13, device = cairo_pdf, height = 4, width = 4)

# K19 correlation plot
d_K19 <- krasddpcams__plot_fitness_correlation_blocks(K19_nor_df, "K19", colour_scheme)
print(d_K19)
ggplot2::ggsave("path/to/K19_fitness_correlation.pdf", d_K19, device = cairo_pdf, height = 4, width = 4)

# Stability correlation plot
d_stability <- krasddpcams__plot_fitness_correlation_blocks(stability_nor_df, "stability", colour_scheme)
print(d_stability)
ggplot2::ggsave("path/to/stability_fitness_correlation.pdf", d_stability, device = cairo_pdf, height = 4, width = 4)

# RAF1 correlation plot
d_RAF1 <- krasddpcams__plot_fitness_correlation_blocks(RAF1_nor_df, "RAF1", colour_scheme)
print(d_RAF1)

ggplot2::ggsave("path/to/RAF1_fitness_correlation.pdf", d_RAF1, device = cairo_pdf, height = 4, width = 4)

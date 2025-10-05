# Fitness Heatmap Visualization for KRAS Mutations
# This script creates heatmaps to visualize fitness effects of single amino acid mutations

library(wlab.block)
library(data.table)
library(ggplot2)

#' Create fitness heatmap for single amino acid mutations
#'
#' @param input Normalized fitness data for single mutations
#' @param wt_aa Wild-type amino acid sequence
#' @param title Plot title
#' @param legend_limits Limits for color legend
#' @return ggplot object showing fitness heatmap

fitness_heatmap_optimization <- function(input, wt_aa, title = "fitness", legend_limits = NULL) {
  # Define amino acid order
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  num <- nchar(wt_aa) + 1
  
  # Prepare single mutation data
  input_single <- input
  input_single[, position := AA_Pos1]
  input_single[, WT_AA := wtcodon1]
  
  # Create heatmap template
  heatmap_tool_fitness <- data.table(
    wtcodon1 = rep(unlist(strsplit(wt_aa, "")), each = 21),
    position = rep(2:num, each = 21),
    codon1 = c(unlist(aa_list), "*")
  )
  
  # Merge with fitness data
  heatmap_tool_fitness_anno_single <- merge(
    input_single, heatmap_tool_fitness,
    by = c("wtcodon1", "position", "codon1"), all = TRUE
  )
  
  # Set factor levels for amino acids
  heatmap_tool_fitness_anno_single <- within(
    heatmap_tool_fitness_anno_single,
    codon1 <- factor(codon1, levels = c(
      "*", "D", "E", "R", "H", "K", "S", "T", "N", "Q",
      "C", "G", "P", "A", "V", "I", "L", "M", "F", "W", "Y"
    ))
  )
  
  # Set wild-type fitness to zero
  heatmap_tool_fitness_anno_single[wtcodon1 == codon1, nor_fitness_nooverlap := 0]
  
  # Create heatmap plot
  ggplot() +
    theme_classic() +
    geom_tile(
      data = heatmap_tool_fitness_anno_single[position > 1, ],
      aes(x = position, y = codon1, fill = nor_fitness_nooverlap)
    ) +
    scale_x_discrete(limits = c(2:num), labels = c(2:num)) +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 0.5, hjust = 0.5,
        color = c(NA, NA, NA, rep(c("black", NA, NA, NA, NA), nchar(wt_aa) %/% 5))
      )
    ) +
    scale_fill_gradient2(
      limits = legend_limits,
      low = "#F4270C", mid = "gray", high = "#1B38A6",
      name = "Fitness",
      na.value = "white",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ylab("Mutant aa") +
    ggtitle(title) +
    labs(fill = NULL) +
    geom_text(
      data = heatmap_tool_fitness_anno_single[position > 1 & wtcodon1 == codon1, ],
      aes(x = position, y = codon1),
      label = "-", size = 3
    ) +
    theme(
      text = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = c(1, 1.38),
      title = element_text(size = 8, face = "bold"),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.text = element_text(size = 8),
      axis.title.x = element_text(size = 8, face = "plain"),
      axis.title.y = element_text(size = 8, face = "plain"),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(
        family = "Courier", angle = 90, size = 9.5,
        vjust = 0.5, hjust = 0.5,
        margin = margin(0, -0.5, 0, 0, "mm")
      ),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(4, "mm"),
      legend.key.size = unit(1, "mm"),
      plot.margin = margin(0, -0, 0, 0)
    ) +
    coord_fixed()
}

# =====================
# Analysis for Different Assays
# =====================

# Define wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Abundance fitness heatmap
nor_fit <- nor_fitness(
  block1 = "path/to/abundance_block1.RData",
  block2 = "path/to/abundance_block2.RData",
  block3 = "path/to/abundance_block3.RData"
)
nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
nor_fit_single <- pos_id(nor_fit_single, wt_aa)
p_abundance <- fitness_heatmap_optimization(nor_fit_single, wt_aa, title = "KRAS-Abundance")
ggsave("path/to/KRAS_Abundance_heatmap.pdf", p_abundance, height = 6, width = 20)

# RAF1 binding fitness heatmap
nor_fit <- nor_fitness(
  block1 = "path/to/RAF1_block1.RData",
  block2 = "path/to/RAF1_block2.RData",
  block3 = "path/to/RAF1_block3.RData"
)
nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
nor_fit_single <- pos_id(nor_fit_single, wt_aa)
p_raf1 <- fitness_heatmap_optimization(nor_fit_single, wt_aa, title = "KRAS-RAF1")
ggsave("path/to/KRAS_RAF1_heatmap.pdf", p_raf1, height = 6, width = 20)

# K13 DARPin binding fitness heatmap
nor_fit <- nor_fitness(
  block1 = "path/to/K13_block1.RData",
  block2 = "path/to/K13_block2.RData",
  block3 = "path/to/K13_block3.RData"
)
nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
nor_fit_single <- pos_id(nor_fit_single, wt_aa)
p_k13 <- fitness_heatmap_optimization(nor_fit_single, wt_aa, title = "KRAS-DARPin K13")
ggsave("path/to/KRAS_K13_heatmap.pdf", p_k13, height = 6, width = 20)

# K19 DARPin binding fitness heatmap
nor_fit <- nor_fitness(
  block1 = "path/to/K19_block1.RData",
  block2 = "path/to/K19_block2.RData",
  block3 = "path/to/K19_block3.RData"
)
nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
nor_fit_single <- pos_id(nor_fit_single, wt_aa)
p_k19 <- fitness_heatmap_optimization(nor_fit_single, wt_aa, title = "KRAS-DARPin K19")
ggsave("path/to/KRAS_K19_heatmap.pdf", p_k19, height = 6, width = 20)
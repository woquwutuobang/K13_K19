library(ggplot2)
library(data.table)
library(dplyr)
library(wlab.block)
library(krasddpcams)
library(ggpubr)
library(rlang) 

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

colour_scheme <- list(
  "blue" = "#1B38A6",        # rgb(27, 56, 166)
  "red" = "#F4270C",         # rgb(244, 39, 12)
  "orange" = "#F4AD0C",      # rgb(244, 173, 12)
  "green" = "#09B636",       # rgb(9, 182, 54)
  "yellow" = "#F1DD10",      # rgb(241, 221, 16)
  "purple" = "#6D17A0",      # rgb(109, 23, 160)
  "pink" = "#FFB0A5",        # rgb(255, 176, 165)
  "light orange" = "#FFE4A5", # rgb(255, 228, 165)
  "light blue" = "#9DACE3",   # rgb(157, 172, 227)
  "light green" = "#97E9AD",  # rgb(151, 233, 173)
  "light red" = "#FF6A56",    # rgb(255, 106, 86)
  "dark red" = "#A31300",     # rgb(163, 19, 0)
  "dark blue" = "#0C226F",    # rgb(12, 34, 111)
  "dark green" = "#007A20"    # rgb(0, 122, 32)
)

# ================================
# Function: Scatter plot (ΔΔGf vs ΔΔGb)
# ================================
krasddpcams__plot_scatter_ddGb_ddGf <- function(ddG1, assay1, ddG2, assay2,
                                                anno, binder = c("K13", "K19"),
                                                colour_scheme,
                                                xlim = c(-1, 3),
                                                ylim = c(-2, 3)) {
  binder <- match.arg(binder)
  
  # ---- 1. Read and merge ddG data ----
  ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1)
  ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2)
  all_ddG <- rbind(ddG1, ddG2)
  all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ assay, value.var = "mean_kcal/mol")
  all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", by.y = "Pos")
  
  # ---- 2. Define binding interface ----
  dist_col <- paste0("scHAmin_ligand_", binder)
  all_ddG_dc_anno[, `:=`(binding_type, "others")]
  all_ddG_dc_anno[get(dist_col) <= 5, `:=`(binding_type, "binding interface")]
  
  # ---- 3. Create plot ----
  p <- ggplot() +
    geom_point(
      data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "others", ],
      aes(x = !!sym(assay1), y = !!sym(binder)),
      color = "black", alpha = 0.6, size = 0.1
    ) +
    geom_point(
      data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "binding interface"],
      aes(x = !!sym(assay1), y = !!sym(binder)),
      color = colour_scheme[["red"]], alpha = 0.6, size = 0.5
    ) +
    
    # ---- Set coordinate limits ----
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    
    theme_classic() +
    labs(
      x = "Folding ΔΔG (kcal/mol)",
      y = paste0(binder, " Binding ΔΔG (kcal/mol)")
    ) +
    theme(
      text = element_text(size = 8),
      legend.position = "right",
      legend.text = element_text(size = 8),
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = 8, colour = "black"
      ),
      axis.text.y = element_text(
        size = 8, colour = "black",
        vjust = 0.5, hjust = 0.5,
        margin = margin(0, 0, 0, 0, "mm")
      ),
      axis.title = element_text(size = 8),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(3.1, "mm"),
      legend.key.size = unit(1, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    coord_fixed(ratio = 1, expand = TRUE, clip = "on")
  
  return(p)
}


anno<-fread("path/to/anno_final_for_5.csv")



krasddpcams__plot_scatter_ddGb_ddGf_K13(ddG1="path/to/weights_Folding.txt",
                                        assay1="folding",
                                        ddG2="path/to/weights_Binding_K13.txt",
                                        assay2="K13",
                                        anno=anno, colour_scheme)

ggplot2::ggsave("figure2a_scatter_ddGb_ddGf_K13.pdf", device = cairo_pdf,height = 4,width=4)


krasddpcams__plot_scatter_ddGb_ddGf_K19(ddG1="path/to/weights_Folding.txt",
                                        assay1="folding",
                                        ddG2="path/to/weights_Binding_K19.txt",
                                        assay2="K19",
                                        anno=anno, colour_scheme)

ggplot2::ggsave("figure2a_scatter_ddGb_ddGf_K19.pdf", device = cairo_pdf,height = 4,width=4)


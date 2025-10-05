library(data.table)
library(wlab.block)
library(scales)  
library(grid)
library(cowplot)
library(ggplot2)



ddG_heatmap_optimization<-function (input, wt_aa, title = "folding free energy change", 
                                    legend_limits = NULL) 
{
  ddG <- fread(input)
  num <- nchar(wt_aa) + 1
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:num, each = 20), mt_codon = unlist(aa_list))
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", 
                                         "mt_codon"), all = T)
  input_heatmap <- within(ddG, mt_codon <- factor(mt_codon, 
                                                  levels = c("D", "E", "R", "H", "K", "S", "T", "N", "Q", 
                                                             "C", "G", "P", "A", "V", "I", "L", "M", "F", "W", 
                                                             "Y")))
  input_heatmap[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 
                                           0)]
  ggplot2::ggplot() + ggpubr::theme_classic2() + 
    ggplot2::geom_tile(data = input_heatmap[Pos_real >1, ], ggplot2::aes(x = Pos_real, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::scale_x_discrete(limits = c(2:num), labels = c(2:num)) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits, 
                                  low = "#1B38A6", mid = "gray", high = "#F4270C", name = expression(Delta * Delta * "G (kcal/mol)" ),
                                  na.value = "white",
                                  guide = ggplot2::guide_colorbar(
                                    title.position = "top",   
                                    title.hjust = 1           
                                  )) + 
    ggplot2::ggtitle(title) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5,vjust = 0.5, hjust = 0.5, color = c(NA, NA, NA, rep(c("black",NA, NA, NA, NA), nchar(wt_aa)%/%5)))) + 
    ggplot2::geom_text(data = input_heatmap[Pos_real > 1 & wt_codon == mt_codon, ], ggplot2::aes(x = Pos_real,y = mt_codon), label = "-", size = 3) + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::ylab("Mutant AA") +
    ggplot2::xlab("Position") +
    ggplot2::theme(text = ggplot2::element_text(size = 8), 
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.ticks.y = ggplot2::element_blank(), 
                   legend.position = c(1,1.40),title = ggplot2::element_text(size = 8,face = "bold"),
                   legend.justification = c(1, 1),
                   legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(size = 8), 
                   axis.title.x = ggplot2::element_text(size = 8, face = "plain"), 
                   axis.title.y = ggplot2::element_text(size = 8, face = "plain"), 
                   axis.text.x = ggplot2::element_text(size = 8, angle = 90, 
                                                       vjust = 0.5, hjust = 1), 
                   axis.text.y = ggplot2::element_text(family = "Courier", angle = 90, size = 9.5, vjust = 0.5, hjust = 0.5, 
                                                       margin = ggplot2::margin(0, -0.5, 0, 0, "mm")), 
                   legend.key.height = ggplot2::unit(3.1, "mm"), 
                   legend.key.width = ggplot2::unit(4, "mm"), 
                   legend.key.size = ggplot2::unit(1, "mm"), plot.margin = ggplot2::margin(0,-0, 0, 0)) + 
    ggplot2::coord_fixed()
}




wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

## folding
ddG_heatmap_optimization(
  input="path/to/weights_Folding.txt",
  wt_aa=wt_aa,
  title="KRAS-Folding free energy changes"
)

ggplot2::ggsave("KRAS-Folding free energy changes.pdf", device = cairo_pdf,height = 6,width=20)



# RAF1
ddG_heatmap_optimization(
  input="path/to/weights_RAF1.txt",
  wt_aa=wt_aa,
  title="KRAS-RAF1 binding free energy changes"
)

ggplot2::ggsave("KRAS-RAF1 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# SOS1
ddG_heatmap_optimization(
  input="path/to/weights_SOS.txt",
  wt_aa=wt_aa,
  title="KRAS-SOS1 binding free energy changes"
)

ggplot2::ggsave("KRAS-SOS1 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# K55
ddG_heatmap_optimization(
  input="path/to/weights_K55.txt",
  wt_aa=wt_aa,
  title="KRAS-K55 binding free energy changes"
)

ggplot2::ggsave("KRAS-K55 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# K27
ddG_heatmap_optimization(
  input="path/to/weights_K27.txt",
  wt_aa=wt_aa,
  title="KRAS-K27 binding free energy changes"
)

ggplot2::ggsave("KRAS-K27 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# RALGDS
ddG_heatmap_optimization(
  input="path/to/weights_RAL.txt",
  wt_aa=wt_aa,
  title="KRAS-RALGDS binding free energy changes"
)

ggplot2::ggsave("KRAS-RALGDS binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# PIK3CG
ddG_heatmap_optimization(
  input="path/to/weights_PI3.txt",
  wt_aa=wt_aa,
  title="KRAS-PIK3CG binding free energy changes"
)

ggplot2::ggsave("KRAS-PIK3CG binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# DARPin K13
ddG_heatmap_optimization(
  input="path/to/weights_K13.txt",
  wt_aa=wt_aa,
  title="KRAS-DARPin_K13 binding free energy changes"
)

ggplot2::ggsave("KRAS-DARPin_K13 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)

# DARPin K19
ddG_heatmap_optimization(
  input="path/to/weights_K19.txt",
  wt_aa=wt_aa,
  title="KRAS-DARPin_K19 binding free energy changes"
)

ggplot2::ggsave("KRAS-DARPin_K19 binding free energy changes.pdf", device = cairo_pdf, height = 6, width = 20)


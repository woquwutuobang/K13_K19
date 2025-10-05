library(ggplot2)
library(data.table)
library(dplyr)
library(wlab.block)
library(krasddpcams)
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

colour_scheme<- list(
  "blue"="#1B38A6",#rgb(27, 56, 166)
  "red"="#F4270C",#rgb(244, 39, 12)
  "orange"="#F4AD0C",#rgb(244, 173, 12)
  "green"="#09B636",#rgb(9, 182, 54)
  "yellow"="#F1DD10",#rgb(241, 221, 16)
  "purple"="#6D17A0",#rgb(109, 23, 160)
  "pink"="#FFB0A5",#rgb(255, 176, 165)
  "light orange"="#FFE4A5",#rgb(255, 228, 165)
  "light blue"="#9DACE3",#rgb(157, 172, 227)
  "light green"="#97E9AD",#rgb(151, 233, 173)
  "light red"="#FF6A56",#rgb(255, 106, 86)
  "dark red"="#A31300",#rgb(163, 19, 0)
  "dark blue"="#0C226F",#rgb(12, 34, 111)
  "dark green"="#007A20"#rgb(0, 122, 32)
)



####krasddpcams__plot_scatter_ddGb_ddGf_K13
krasddpcams__plot_scatter_ddGb_ddGf_K13<-function (ddG1 = ddG1, assay1_sele = assay1_sele, ddG2 = ddG2, 
                                                   assay2_sele = assay2_sele, anno = anno, colour_scheme) 
{
  ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1_sele)
  ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2_sele)
  all_ddG <- rbind(ddG1, ddG2)
  all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ 
                        assay, value.var = "mean_kcal/mol")
  all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", 
                           by.y = "Pos")
  all_ddG_dc_anno[, `:=`(K13_type_bs, "others")]
  all_ddG_dc_anno[scHAmin_ligand_K13 <= 5, `:=`(K13_type_bs, 
                                                "binding interface")]
  ggplot2::ggplot() + ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                                   1 & K13_type_bs == "others", ], ggplot2::aes(x = folding, 
                                                                                                                y = K13), color = "black", alpha = 0.6, size = 0.1) + 
    ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                 1 & K13_type_bs == "binding interface"], ggplot2::aes(x = folding, 
                                                                                                       y = K13), color = colour_scheme[["red"]], alpha = 0.6, 
                        size = 0.5) + ggpubr::theme_classic2() + ggplot2::labs(color = NULL) + 
    ggplot2::xlab("Folding ddG") + ggplot2::ylab("Binding ddG") + 
    ggplot2::theme(text = ggplot2::element_text(size = 8), 
                   legend.position = "right", legend.text = ggplot2::element_text(size = 8), 
                   axis.text.x = ggplot2::element_text(size = 8, vjust = 0.5, 
                                                       hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8, 
                                                                                                         vjust = 0.5, hjust = 0.5, margin = ggplot2::margin(0, 
                                                                                                                                                            0, 0, 0, "mm")), legend.key.height = ggplot2::unit(3.1, 
                                                                                                                                                                                                               "mm"), legend.key.width = ggplot2::unit(3.1, 
                                                                                                                                                                                                                                                       "mm"), legend.key.size = ggplot2::unit(1, "mm"),
                   
                   plot.margin = ggplot2::margin(0, 0, 0, 0)) + ggplot2::coord_fixed(ratio = 1, 
                                                                                     xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}


####krasddpcams__plot_scatter_ddGb_ddGf_K13
krasddpcams__plot_scatter_ddGb_ddGf_K13<-function (ddG1 = ddG1, assay1_sele = assay1_sele, ddG2 = ddG2, 
                                                   assay2_sele = assay2_sele, anno = anno, colour_scheme) 
{
  ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1_sele)
  ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2_sele)
  all_ddG <- rbind(ddG1, ddG2)
  all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ 
                        assay, value.var = "mean_kcal/mol")
  all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", 
                           by.y = "Pos")
  all_ddG_dc_anno[, `:=`(K13_type_bs, "others")]
  all_ddG_dc_anno[scHAmin_ligand_K13 <= 5, `:=`(K13_type_bs, 
                                                "binding interface")]
  ggplot2::ggplot() + ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                                   1 & K13_type_bs == "others", ], ggplot2::aes(x = folding, 
                                                                                                                y = K13), color = "black", alpha = 0.6, size = 0.1) + 
    ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                 1 & K13_type_bs == "binding interface"], ggplot2::aes(x = folding, 
                                                                                                       y = K13), color = colour_scheme[["red"]], alpha = 0.6, 
                        size = 0.5) + ggpubr::theme_classic2() + ggplot2::labs(color = NULL) + 
    ggplot2::xlab("Folding ddG") + ggplot2::ylab("K13 Binding ddG") + 
    ggplot2::theme(text = ggplot2::element_text(size = 8), 
                   legend.position = "right", legend.text = ggplot2::element_text(size = 8), 
                   axis.text.x = ggplot2::element_text(size = 8, vjust = 0.5, 
                                                       hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8, 
                                                                                                         vjust = 0.5, hjust = 0.5, margin = ggplot2::margin(0, 
                                                                                                                                                            0, 0, 0, "mm")), legend.key.height = ggplot2::unit(3.1, 
                                                                                                                                                                                                               "mm"), legend.key.width = ggplot2::unit(3.1, 
                                                                                                                                                                                                                                                       "mm"), legend.key.size = ggplot2::unit(1, "mm"),
                   
                   plot.margin = ggplot2::margin(0, 0, 0, 0)) + ggplot2::coord_fixed(ratio = 1, 
                                                                                     xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}




####krasddpcams__plot_scatter_ddGb_ddGf_K19
krasddpcams__plot_scatter_ddGb_ddGf_K19<-function (ddG1 = ddG1, assay1_sele = assay1_sele, ddG2 = ddG2, 
                                                   assay2_sele = assay2_sele, anno = anno, colour_scheme) 
{
  ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1_sele)
  ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2_sele)
  all_ddG <- rbind(ddG1, ddG2)
  all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ 
                        assay, value.var = "mean_kcal/mol")
  all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", 
                           by.y = "Pos")
  all_ddG_dc_anno[, `:=`(K19_type_bs, "others")]
  all_ddG_dc_anno[scHAmin_ligand_K19 <= 5, `:=`(K19_type_bs, 
                                                "binding interface")]
  ggplot2::ggplot() + ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                                   1 & K19_type_bs == "others", ], ggplot2::aes(x = folding, 
                                                                                                                y = K19), color = "black", alpha = 0.6, size = 0.1) + 
    ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
                                                 1 & K19_type_bs == "binding interface"], ggplot2::aes(x = folding, 
                                                                                                       y = K19), color = colour_scheme[["red"]], alpha = 0.6, 
                        size = 0.5) + ggpubr::theme_classic2() + ggplot2::labs(color = NULL) + 
    ggplot2::xlab("Folding ddG") + ggplot2::ylab("K19 Binding ddG") + 
    ggplot2::theme(text = ggplot2::element_text(size = 8), 
                   legend.position = "right", legend.text = ggplot2::element_text(size = 8), 
                   axis.text.x = ggplot2::element_text(size = 8, vjust = 0.5, 
                                                       hjust = 0.5), axis.text.y = ggplot2::element_text(size = 8, 
                                                                                                         vjust = 0.5, hjust = 0.5, margin = ggplot2::margin(0, 
                                                                                                                                                            0, 0, 0, "mm")), legend.key.height = ggplot2::unit(3.1, 
                                                                                                                                                                                                               "mm"), legend.key.width = ggplot2::unit(3.1, 
                                                                                                                                                                                                                                                       "mm"), legend.key.size = ggplot2::unit(1, "mm"),
                   
                   plot.margin = ggplot2::margin(0, 0, 0, 0)) + ggplot2::coord_fixed(ratio = 1, 
                                                                                     xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
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

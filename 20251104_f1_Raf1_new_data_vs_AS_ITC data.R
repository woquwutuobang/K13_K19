library(krasddpcams)

colour_scheme<-list(
  "blue"="#1B38A6",#rgb(27, 56, 166)
  "red"="#F4270C",#rgb(244, 39, 12)
  "orange"="#F4AD0C",#rgb(244, 173, 12)
  "green"="#09B636",#rgb(9, 182, 54)
  "yellow"="#F1DD10",#rgb(241, 221, 16)
  "purple"="#C68EFD",#rgb(198, 142, 253)
  "hot pink"="#FF0066",#rgb(255, 0, 102)
  "light blue"="#75C2F6",#rgb(117, 194, 246)
  "light red"="#FF6A56",#rgb(255, 106, 86)       # The red ones are unified with this, but the heatmap ones remain unchanged.
  "dark red"="#A31300",#rgb(163, 19, 0)
  "dark green"="#007A20",#rgb(0, 122, 32)
  "pink"="#FFB0A5" #rgb(255, 176, 165)
)


krasddpcams__plot_RAF_invitro_cor(input="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
                                  assay_name = "RAF1", colour_scheme)


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f1/20251104/figure1k_RAF1_new_invitro_ddG_cor.pdf", device = cairo_pdf,height = 50,width=50,units = "mm")

library(wlab.block)
library(data.table)
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

nor_fit<-nor_fitness(block1="../../../final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                     block2="../../../final_fitness_for_plot/20250525_RDatas/Abundance_block2_Q20_fitness_replicates_fullseq.RData",
                     block3="../../../final_fitness_for_plot/20250525_RDatas/Abundance_block3_Q20_rbg_filter_fitness_replicates_fullseq_20250710.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap_optimization(nor_fit_single,wt_aa,title = "KRAS-Abundance")
ggplot2::ggsave("./result/20250721/KRAS-Abundance.pdf",height = 6,width=20)







fitness_heatmap_optimization<-function (input, wt_aa, title = "fitness", legend_limits = NULL) 
{
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  num <- nchar(wt_aa) + 1
  input_single <- input
  input_single[, `:=`(position, AA_Pos1)]
  input_single[, `:=`(WT_AA, wtcodon1)]
  heatmap_tool_fitness <- data.table(wtcodon1 = rep(unlist(strsplit(wt_aa, 
                                                                    "")), each = 21), position = rep(2:num, each = 21), codon1 = c(unlist(aa_list), 
                                                                                                                                   "*"))
  heatmap_tool_fitness_anno_single <- merge(input_single, heatmap_tool_fitness, 
                                            by = c("wtcodon1", "position", "codon1"), all = T)
  heatmap_tool_fitness_anno_single <- within(heatmap_tool_fitness_anno_single, 
                                             codon1 <- factor(codon1, levels = c("*", "D", "E", "R", 
                                                                                 "H", "K", "S", "T", "N", "Q", "C", "G", "P", "A", 
                                                                                 "V", "I", "L", "M", "F", "W", "Y")))
  heatmap_tool_fitness_anno_single[wtcodon1 == codon1, `:=`(nor_fitness_nooverlap, 
                                                            0)]
  ggplot2::ggplot() + ggplot2::theme_classic() + ggplot2::geom_tile(data = heatmap_tool_fitness_anno_single[position > 
                                                                                                              1, ], ggplot2::aes(x = position, y = codon1, fill = nor_fitness_nooverlap)) + 
    ggplot2::scale_x_discrete(limits = c(2:num), labels = c(2:num)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, 
                                                       vjust = 0.5, hjust = 0.5, color = c(NA, NA, NA, rep(c("black", 
                                                                                                             NA, NA, NA, NA), nchar(wt_aa)%/%5)))) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits, 
                                  low = "#F4270C", mid = "gray", high = "#1B38A6", name = "Fitness",
                                  na.value = "white",
                                  guide = ggplot2::guide_colorbar(
                                    title.position = "top",   # 标题在上方
                                    title.hjust = 0.5           # 标题右对齐（0左，0.5中，1右）
                                  )) + 
    ggplot2::ylab("Mutant aa") + 
    ggplot2::ggtitle(title) + 
    ggplot2::labs(fill = NULL) + ggplot2::geom_text(data = heatmap_tool_fitness_anno_single[position > 
                                                                                              1 & wtcodon1 == codon1, ], ggplot2::aes(x = position, 
                                                                                                                                      y = codon1), label = "-", size = 3) + 
    ggplot2::theme(text = ggplot2::element_text(size = 8),
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = c(1,1.38),title = ggplot2::element_text(size = 8,face = "bold"),
                   legend.justification = c(1, 1),
                   legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(size = 8), 
                   axis.title.x = ggplot2::element_text(size = 8, face = "plain"), 
                   axis.title.y = ggplot2::element_text(size = 8, face = "plain"), 
                   axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1), 
                   axis.text.y = ggplot2::element_text(family = "Courier", angle = 90, size = 9.5, vjust = 0.5, hjust = 0.5, margin = ggplot2::margin(0, -0.5, 0, 0, "mm")), 
                   legend.key.height = ggplot2::unit(3.1, "mm"), legend.key.width = ggplot2::unit(4, "mm"), 
                   legend.key.size = ggplot2::unit(1, "mm"), plot.margin = ggplot2::margin(0, -0, 0, 0)) + 
    ggplot2::coord_fixed()
}

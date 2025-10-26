library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"


# Secondary structure information from the paper
rects_sheet <- data.frame(xstart = c(3,38,51,77,109,139),
                          xend = c(9,44,57,84,115,143),
                          col = c("β1","β2","β3","β4","β5","β6"))

rects_helix <- data.frame(xstart = c(15,67,87,127,148),
                          xend = c(24,73,104,136,166),
                          col = c("α1","α2","α3","α4","α5"))

# Function for plotting one assay (either K13 or K19)
Manhatta_plot_single_assay_for_BI1_binder <- function(input, assay_sele, anno, rects_sheet, rects_alpha, wt_aa) {
  # Read and process ddG data
  ddG <- fread(input)
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  anno <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  # Define binding site types
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # Merge mutation data with site type information
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # Classify mutation types
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Orthosteric site huge differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "Orthosteric site small differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "GTP binding allosteric mutation"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "GTP binding other mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == TRUE, mutation_type := "Allosteric mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == FALSE, mutation_type := "Other mutation"]
  
  # Set factor levels for mutation types
  data_plot_mutation <- within(data_plot_mutation,
                               mutation_type <- factor(mutation_type,
                                                       levels = c("Orthosteric site huge differences",
                                                                  "Orthosteric site small differences",
                                                                  "GTP binding allosteric mutation",
                                                                  "GTP binding other mutation",
                                                                  "Allosteric mutation",
                                                                  "Other mutation")))
  
  # Create Manhattan plot
  p <- ggplot() +
    # Secondary structure background
    geom_rect(data = rects_sheet, aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5),
              fill = "#75C2F6", alpha = 0.06) +
    geom_rect(data = rects_alpha, aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5),
              fill = "#C68EFD", alpha = 0.06) +
    # Mutation points
    geom_point(data = data_plot_mutation,
               aes(x = Pos_real, y = `mean_kcal/mol`, color = mutation_type), size = 1) +
    # Color scheme with transparency
    scale_color_manual(
      values = c(
        alpha("#F4270C", 1),    # 100% opacity   BI1 color
        alpha("#FFB0A5", 1),  # 10% opacity
        alpha("#007A20", 1),    # 100% opacity
        alpha("#09B636", 0.5),  # 60% opacity
        alpha("grey", 0.8),    # 100% opacity
        alpha("gray", 0.8)      # 80% opacity
      )
    ) +
    # Reference line
    geom_hline(yintercept = 0, linetype = 2) +
    # Axis settings
    scale_x_continuous(expand = c(1 / 188, 11 / 188)) +
    ylab(paste0("Binding ΔΔGb (", assay_sele, ") (kcal/mol)")) +
    xlab("Amino acid position") +
    labs(color = NULL) +
    # Secondary structure labels (rotated 90 degrees)
    annotate("text", x = (rects_sheet$xstart[1] + rects_sheet$xend[1]) / 2, y = 3.1, 
             label = "strand", size = 2.5, vjust = 0.5, hjust = 1, angle = 90) +
    annotate("text", x = (rects_helix$xstart[1] + rects_helix$xend[1]) / 2, y = 3.1, 
             label = "helix", size = 2.5, vjust = 0.5, hjust = 1, angle = 90) +
    # Theme settings
    theme_classic2() +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          text = element_text(size = 8),
          legend.position = "none",
          legend.text = element_text(size = 8),
          strip.background = element_rect(colour = "black", fill = "white")) +
    coord_fixed(ratio = 10, xlim = c(-0.5, 190), ylim = c(-1.5, 3.5))
  
  return(p)
}



## control
p_RAF<-Manhatta_plot_single_assay_for_BI1_binder(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_RAF.txt",
  assay_sele = "RAF",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final3.csv",
  rects_sheet = rects_sheet,
  rects_alpha = rects_helix,
  wt_aa = wt_aa
)
p_RAF


# BI 1 RAF1
p_RAF1<-Manhatta_plot_single_assay_for_BI1_binder(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  rects_sheet = rects_sheet,
  rects_alpha = rects_helix,
  wt_aa = wt_aa
)

p_RAF1

# Save the plot
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/figure2b_new_RAF1_Manhattan_plot.pdf", p_RAF1,
                device = cairo_pdf, height = 4.5, width = 12)





# BI 1 K27
p_K27<-Manhatta_plot_single_assay_for_BI1_binder(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  assay_sele = "K27",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  rects_sheet = rects_sheet,
  rects_alpha = rects_helix,
  wt_aa = wt_aa
)

p_K27

# Save the plot
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure2/20251026/figure2b_new_K27_Manhattan_plot.pdf", p_K27,
                device = cairo_pdf, height = 4.5, width = 12)

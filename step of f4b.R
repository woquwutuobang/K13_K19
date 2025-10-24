
## new_f4c
##20251024



## first panel  versions 1
## show each allosteric sites

## 位点信息---all allosteric sites
K13_allosteric_site <- c(15, 16, 17, 35, 145, 10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)


K13_RAF1_common_allosteric_sites<-c(10,15,16,17,35,55,77,79,145,159,163)

K13_unique_allosteric_sites<-c(19,21,24,53,82,93,151)
RAF1_unique_allosteric_sites<-c(20,28,32,34,54,57,58,60,112,113,114,134,144,146,156)


## color setting:K13 unique----"#09B636";RAF1 unique-----"#F4AD0C";Both:----"#C68EFD"

## open "C:\\Users\\36146\\OneDrive - USTC\\DryLab\\base_information_for_K13_K19_project\\6vjj.pdb"
## set bgColor #ffffff00
## open "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h46.pdb"
## matchmaker #2 to #1;hide solvent;color #1/B black;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;color #1/A lightgrey
## sel #1/A:20,28,32,34,54,57,58,60,112,113,114,134,144,146,156;color sel #F4AD0C;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #1/A:10,15,16,17,35,55,77,79,145,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## select :GNP;surface sel;color sel #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/RAF1_allosites.cxs"

## color #2/B black;color #2/A lightgrey
## sel #2/A:19,21,24,53,82,93,151;color sel #09B636;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #2/A:10,15,16,17,35,55,77,79,145,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## select :GDP;surface sel;color sel #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/K13_allosites.cxs"

## 两个结构merge版本
## hide #2/A cartoon;hide #2/A atoms
## sel #1/A:19,21,24,53,82,93,151;color sel #09B636;show sel atoms;style sel sphere;size sel atomRadius 2;
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/K13_RAF1_merge_major_allosites.cxs"






## first panel  versions 2
## show each allosteric sites

## 位点信息---major allosteric sites
K13_allosteric_site <- c(10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)


K13_RAF1_common_allosteric_sites<-c(10,55,77,79,159,163)

K13_unique_allosteric_sites<-c(19,21,24,53,82,93,151)
RAF1_unique_allosteric_sites<-c(20,54,58,112,113,114,134,144,156)

## color setting:K13 unique----"#09B636";RAF1 unique-----"#F4AD0C";Both:----"#C68EFD"

## open "C:\\Users\\36146\\OneDrive - USTC\\DryLab\\base_information_for_K13_K19_project\\6vjj.pdb"
## set bgColor #ffffff00
## open "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6h46.pdb"
## matchmaker #2 to #1;hide solvent;color #1/B black;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;color #1/A lightgrey
## sel #1/A:20,54,58,112,113,114,134,144,156;color sel #F4AD0C;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #1/A:10,55,77,79,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## select :GNP;surface sel;color sel #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/RAF1_major_allosites.cxs"

## color #2/B black;color #2/A lightgrey
## sel #2/A:19,21,24,53,82,93,151;color sel #09B636;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #2/A:10,15,16,17,35,55,77,79,145,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## select :GDP;surface sel;color sel #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/K13_major_allosites.cxs"

## 两个结构merge版本
## hide #2/A cartoon;hide #2/A atoms
## sel #1/A:19,21,24,53,82,93,151;color sel #09B636;show sel atoms;style sel sphere;size sel atomRadius 2;
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/K13_RAF1_merge_allosites.cxs"






### median ddGb of K13/RAF1 to structure

### second panel& third panel
## RAF1&K13
## open "C:\\Users\\36146\\OneDrive - USTC\\DryLab\\Data_analysis_scripts\\median ddG endow different type structure\\results\\20250904_median for structure\\median\\6vjj_RAF1_median_ddG_use_0901_data_2.pdb"
## set bgColor white;
## open "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/median ddG endow different type structure/results/20250904_median for structure/median/6h46_K13_median_ddG_use_0901_data.pdb"
## matchmaker #2 to #1;hide solvent;color #1/B black;color #2/B black;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;
## color bfactor #1/A range -1,1;color bfactor #2/A range -1,1;
## select :GNP;surface sel;color sel #1B38A6;select :GDP;surface sel;color sel #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/RAF1_median_to_structure.cxs"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure4/20251024/K13_median_to_structure.cxs"





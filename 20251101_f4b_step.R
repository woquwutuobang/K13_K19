## 20251101
## f4b

## show each allosteric sites

## 位点信息---all allosteric sites
K13_allosteric_site <- c(15, 16, 17, 35, 145, 10, 19, 21, 24, 53, 55, 77, 79, 82, 93, 151, 159, 163)
RAF1_allosteric_site <- c(15, 16, 17, 28, 32, 34, 35, 57, 60, 145, 146, 10, 20, 54, 55, 58, 77, 79, 112, 113, 114, 134, 144, 156, 159, 163)

##########===========
## common allosteric sites
K13_RAF1_common_allosteric_sites<-c(10,15,16,17,35,55,77,79,145,159,163)

## unique allosteric sites
K13_unique_allosteric_sites<-c(19,21,24,53,82,93,151)
RAF1_unique_allosteric_sites<-c(20,28,32,34,54,57,58,60,112,113,114,134,144,146,156)


## color setting:K13 unique----"#FFB0A5";RAF1 unique-----"#F4AD0C";Both:----"#C68EFD"

## K13
## open "C:\\Users\\36146\\OneDrive - USTC\\DryLab\\base_information_for_K13_K19_project\\6h46.pdb"
## set bgColor #ffffff00
## open "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/6vjj.pdb"
## matchmaker #2 to #1;hide solvent;color #1/B black;lighting soft;graphics silhouettes true color black width 1;cartoon style width 2.5 thickness 0.8;color #1/A lightgrey;hide #1/B atoms;
## select :GDP;surface sel;color sel #09B636
## sel #1/A:19,21,24,53,82,93,151;color sel #FFB0A5;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #1/A:10,15,16,17,35,55,77,79,145,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## color #1/B #1B38A6
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/K13_allosteric_sites_show_with_cartoon.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K13_allosteric_sites_show_with_cartoon.mp4




## RAF1
## color #2/B black;color #2/A lightgrey;hide #2/B atoms
## select :GNP;surface sel;color sel #09B636
## sel #2/A:20,28,32,34,54,57,58,60,112,113,114,134,144,146,156;color sel #F4AD0C;show sel atoms;style sel sphere;size sel atomRadius 2;
## sel #2/A:10,15,16,17,35,55,77,79,145,159,163;color sel #C68EFD;show sel atoms;style sel sphere;size sel atomRadius 2;
## color #2/B #F4270C
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f4/20251101/RAF1_allosteric_sites_show_with_cartoon.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\RAF1_allosteric_sites_show_with_cartoon.mp4






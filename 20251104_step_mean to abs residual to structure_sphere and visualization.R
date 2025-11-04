## 20251104
## sphere structure and visualization 8 BPs
## 用位点mean绝对值进行拟合，然后mean减去拟合值得到残差，|残差|赋值结构


## BI2
## K13(front view)
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\sf8\\20251104\\K13_residuals.pdb"
## set bgColor white;hide #1/A cartoon;show #1/A atoms;style #1/A sphere;color bfactor #1/A range -1,1;color #1/B black;lighting full;graphics silhouettes true color black width 1;sel #1/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #1/A atomRadius 2.2;hide #1/B atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K13_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K13_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View" 
## set bgColor white;hide #1/A cartoon;show #1/A atoms;style #1/A sphere;color bfactor #1/A range -1,1;color #1/B black;lighting full;graphics silhouettes true color black width 1;sel #1/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #1/A atomRadius 2.2;hide #1/B atoms
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K13_mean_to_residual_sphere_visuliazation.mp4





## K19
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\K13_mean_to_residual_sphere.cxs"
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\sf8\\20251104\\K19_residuals.pdb"
## matchmaker #!2 to #1;
## set bgColor white;hide #2/A cartoon;show #2/A atoms;style #2/A sphere;color bfactor #2/A range -1,1;color #2/B black;lighting full;graphics silhouettes true color black width 1;sel #2/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/A atomRadius 2.2;hide #2/B atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K19_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K19_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K19_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K19_mean_to_residual_sphere_visuliazation.mp4





## BI1
## RAF1(front view) 
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\sf8\\20251104\\RAF1_residuals.pdb"
## set bgColor white;hide #1/A cartoon;show #1/A atoms;style #1/A sphere;color bfactor #1/A range -1,1;color #1/B black;lighting full;graphics silhouettes true color black width 1;sel #1/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #1/A atomRadius 2.2;hide #1/B atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/RAF1_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\RAF1_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/RAF1_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\RAF1_mean_to_residual_sphere_visuliazation.mp4






## RALGDS
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\RAF1_mean_to_residual_sphere.cxs"
## open "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/RALGDS_residuals.pdb"
## matchmaker #2/B to #1/A pairing ss;hide #2/C,D cartoon;hide #2/C,D atoms;
## set bgColor white;hide #2/B cartoon;show #2/B atoms;style #2/B sphere;color bfactor #2/B range -1,1;color #2/A black;lighting full;graphics silhouettes true color black width 1;sel #2/A;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/B atomRadius 2.2;hide #2/A atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/RALGDS_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\RALGDS_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/RALGDS_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\RALGDS_mean_to_residual_sphere_visuliazation.mp4




## PI3KCG
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\RAF1_mean_to_residual_sphere.cxs"
## open "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/PI3KCG_residuals.pdb"
## matchmaker #2/B to #1/A pairing ss;
## set bgColor white;hide #2/B cartoon;show #2/B atoms;style #2/B sphere;color bfactor #2/B range -1,1;color #2/A black;lighting full;graphics silhouettes true color black width 1;sel #2/A;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/B atomRadius 2.2;hide #2/A atoms
## select #2/A:144-202;hide sel cartoons;select #2/A:313-1084;hide sel cartoons
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/PI3KCG_RBD_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\PI3KCG_RBD_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/PI3KCG_RBD_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\PI3KCG_RBD_mean_to_residual_sphere_visuliazation.mp4




## SOS1
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\RAF1_mean_to_residual_sphere.cxs"
## open "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/SOS1_Q_residuals.pdb"
## matchmaker #2/Q to #1/A pairing ss;hide #2/R cartoon;hide #2/S atoms;hide #2/R atoms
## set bgColor white;hide #2/Q cartoon;show #2/Q atoms;style #2/Q sphere;color bfactor #2/Q range -1,1;color #2/S black;lighting full;graphics silhouettes true color black width 1;sel #2/S;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/S atomRadius 2.2;hide #2/S atoms;hide #2/R atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/SOS1_Q_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\SOS1_Q_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/SOS1_Q_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\SOS1_Q_mean_to_residual_sphere_visuliazation.mp4







## K55
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\RAF1_mean_to_residual_sphere.cxs"
## open "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K55_residuals.pdb"
## matchmaker #!2 to #1
## set bgColor white;hide #2/A cartoon;show #2/A atoms;style #2/A sphere;color bfactor #2/A range -1,1;color #2/B black;lighting full;graphics silhouettes true color black width 1;sel #2/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/A atomRadius 2.2;hide #2/B atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K55_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K55_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K55_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K55_mean_to_residual_sphere_visuliazation.mp4




## K27
## open "C:\\Users\\36146\\OneDrive - USTC\\Manuscripts\\K13_K19\\figures\\20251031_version_all_figure\\f6\\20251104\\RAF1_mean_to_residual_sphere.cxs"
## open "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/sf8/20251104/K27_residuals.pdb"
## matchmaker #!2 to #1;hide #2/C,E,D,F,H,G atoms;show #2/B cartoon;hide #2/B atoms
## set bgColor white;hide #2/A cartoon;show #2/A atoms;style #2/A sphere;color bfactor #2/A range -1,1;color #2/B black;lighting full;graphics silhouettes true color black width 1;sel #2/B;cartoon style width 2.5 thickness 0.8;hide solvent;size #2/A atomRadius 2.2;hide #2/B atoms
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K27_mean_to_residual_sphere.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K27_mean_to_residual_sphere.mp4
## turn y 270;ui tool show "Side View"
## save "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f6/20251104/K27_mean_to_residual_sphere_visulization.cxs"
## movie record;turn y 3 360;wait 200;movie stop
## movie encode C:\Users\36146\OneDrive\K27_mean_to_residual_sphere_visuliazation.mp4












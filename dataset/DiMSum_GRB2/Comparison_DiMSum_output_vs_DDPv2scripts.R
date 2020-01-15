### Comparison v2 vs v3 fitness scores DDP ###
# JDE October 2019

#modified 2019/11/28 by Joern adding also comparison of uncertainties

#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer


require(data.table)
require(ggplot2)
theme_set(theme_classic(base_size=9))
require(gridExtra)


### load GRB2 single mutation data - Old version of fitness scores
# use  source data with read count threshold 20 and non-logged fitness values 
singles_binding = fread("dataset/GRB2_CYC_singles_readT20.txt")
singles_stab = fread("dataset/GRB2_GPD_singles_readT20.txt")
# merge both datasets (binding and stability assays)
singles = merge(singles_binding[Mut!='*',.(Pos,WT_AA,Mut, sigma_binding = sigma,binding = fitness)],singles_stab[Mut!='*',.(Pos,WT_AA,Mut,sigma_stab= sigma,stability= fitness)],by=c("Pos","WT_AA","Mut"),all=T)
singles[,diff_binding_stability := (binding-stability)]
singles[,diff_stability_binding := (stability-binding)]
singles[,sigma_diff := sqrt((sigma_binding^2) + (sigma_stab^2))]
singles[,Mut := factor(Mut,levels=strsplit('RHKDESTNQCGPAVILMFYW','')[[1]])]
#few singles with infinite binding values, kick them out
singles[is.infinite(binding),binding := NA]
singles = singles[!is.na(binding)]
singles = singles[!is.na(stability)] 

### load GRB2 single mutatnt data - New version of scores (directly from DiMSum output)
singles_binding_2 = fread("dataset/DiMSum_last_version/fitness_singles_CYC.txt")
singles_stab_2 = fread("dataset/DiMSum_last_version/fitness_singles_GPD.txt")

singles_binding_2 = singles_binding_2[is.reads0 == T & mean_count > 20 & Mut !="*" & WT_AA != "*",]
singles_stab_2 = singles_stab_2[is.reads0 == T & mean_count > 20 & Mut !="*" & WT_AA != "*",]


### merge old vs new
bind <- merge(singles_binding, singles_binding_2, by = c("WT_AA", "Pos", "Mut"), suffixes = c("_old", "_new") )
stab <- merge(singles_stab, singles_stab_2, by = c("WT_AA", "Pos", "Mut"), suffixes = c("_old", "_new") )

conversion_binding_new2old = bind[,nls(formula = fitness_old ~ a*fitness_new + b,start = list(a=0.3,b=1))]
conversion_stability_new2old = stab[,nls(formula = fitness_old ~ a*fitness_new + b,start = list(a=0.3,b=1))]

p1a <- ggplot(bind, aes(fitness_new, fitness_old)) + geom_point() +
  ggtitle(paste0("deepPCA, spearman rho = ", round(cor(bind$fitness_old, bind$fitness_new, method = "spearman"), 2))) 
p1b <- ggplot(stab, aes(fitness_new, fitness_old)) + geom_point() +
  ggtitle(paste0("stabilityPCA, spearman rho = ", round(cor(stab$fitness_old, stab$fitness_new, method = "spearman"), 2)))
p1c <- ggplot(bind, aes(coef(conversion_binding_new2old)[1]*sigma_new, sigma_old)) + geom_point() +
  ggtitle(paste0("deepPCA, spearman rho = ", round(cor(bind$sigma_old, bind$sigma_new, method = "spearman"), 2))) + 
  scale_x_log10() + scale_y_log10() + geom_abline(color="red")
p1d <- ggplot(stab, aes(coef(conversion_stability_new2old)[1]*sigma_new, sigma_old)) + geom_point() +
  ggtitle(paste0("stabilityPCA, spearman rho = ", round(cor(stab$sigma_old, stab$sigma_new, method = "spearman"), 2)))  + 
  scale_x_log10() + scale_y_log10() + geom_abline(color="red")

grid.arrange(p1a, p1b,p1c,p1d)


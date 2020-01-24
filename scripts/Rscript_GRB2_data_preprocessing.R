#Rscript for preprocessing GRB2 data to use in all other scripts
# this is the first stuff we did in the figure script, but now for the new DiMSum data with proper error model


#### version: 1.0
#### created: 2019/11/28 by Joern
#### last modified by: Julia

#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")

require(data.table)
require(ggplot2)
require(gridExtra)
require(cowplot)
library(ggExtra)
require(PRROC)
theme_set(theme_classic(base_size=9))

#########################################################
################# LOAD AND PROCESS DATA ################# 
#########################################################

### load GRB2 single mutation data
# use  source data with read count threshold 20 and non-logged fitness values
singles_binding = fread("dataset/DiMSum_last_version/GRB2_epPCR_CYC_fitness_singles.txt")
singles_stab = fread("dataset/DiMSum_last_version/GRB2_epPCR_GPD_fitness_singles.txt")

# merge both datasets (binding and stability assays)
singles = merge(singles_binding[Mut!='*' & mean_count > 20,.(Pos,WT_AA,Mut,mean_count,sigma_binding = sigma,binding = fitness)],
                singles_stab[Mut!='*' & mean_count > 20,.(Pos,WT_AA,Mut,mean_count,sigma_stab= sigma,stability= fitness)],by=c("Pos","WT_AA","Mut"),all=T)
singles[,diff_binding_stability := (binding-stability)]
singles[,diff_stability_binding := (stability-binding)]
singles[,sigma_diff := sqrt((sigma_binding^2) + (sigma_stab^2))]
singles[,Mut := factor(Mut,levels=strsplit('RHKDESTNQCGPAVILMFYW','')[[1]])]

# few singles with infinite binding values, kick them out
singles[is.infinite(binding),binding := NA]
singles = singles[!is.na(binding)]
singles = singles[!is.na(stability)] 


### load distances between GRB2 and ligand
GRB2_GAB2_distances = fread("dataset/PDB_contactmap_2vwf_AB.txt")
GRB2_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
#add minimal side-chain heavy atom distance to ligand to singles data.table
singles = merge(singles,GRB2_GAB2_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1],by.x="Pos",by.y="Pos1")


### load RSA values
RSA_chainA = fread("dataset/2vwf_A.rsa",skip=8,nrows = 56)
names(RSA_chainA) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
RSA_chainAB = fread("dataset/2vwf_AB.rsa",skip=8,nrows = 56+14)
names(RSA_chainAB) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
# add to singles data.table
singles = merge(singles,RSA_chainA[chain == "A",.(Pos,RSA_unbound = RSA_all_rel)],by="Pos")
singles = merge(singles,RSA_chainAB[chain == "A",.(Pos,RSA_bound = RSA_all_rel)],by="Pos")
# >>> all residues that loose RSA when bound to ligand are on interface with ligand
# ggplot(singles,aes(RSA_unbound,RSA_bound,color=HAmin_ligand < 5)) + geom_point(size=3)

## split data into core, surface and binding
# ggplot(singles,aes(RSA_unbound,stability)) + geom_point() +geom_smooth()
# >>> set core threshodl as RSA <= 10
threshold_RSA = 10
threshold_ligand_dist = 5 
singles[RSA_unbound <= threshold_RSA,type := "core"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "ligand_binding"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist,type := "surface"]
singles[,type := factor(type,levels=c("core","surface","ligand_binding"))]
singles[,.(Npos = length(unique(Pos))),type]

### write table of processed single mutants
write.table(x=singles,file = "dataset/DiMSum_last_version/GRB2_singles_dataset.txt",col.names = T,row.names=F,quote=F)


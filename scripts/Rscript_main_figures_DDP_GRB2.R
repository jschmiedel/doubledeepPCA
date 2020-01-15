##### double deep main figure plots ######

#### version: 4.0
#### last modified: 2019/11/11
#### created by: Joern
#### last modified by: Julia
# all main figures in a single script (splited old fig2 into 3 figures)
# fixed bug in p-value calculations
# add new panel A in fig4


#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")

require(data.table)
require(ggplot2)
theme_set(theme_classic(base_size=9))
require(gridExtra)
require(cowplot)
library(ggExtra)
require(PRROC)


## Colors
col_purple = "#9161A8"
col_blue =  "#0066CC"
col_orange = "#F7941E"
col_red = "#EF4136"



#########################################################
################# LOAD AND PROCESS DATA ################# 
#########################################################

### load GRB2 single mutation data
# use  source data with read count threshold 20 and non-logged fitness values
singles_binding = fread("dataset/GRB2_CYC_singles_readT20.txt")
singles_stab = fread("dataset/GRB2_GPD_singles_readT20.txt")

# merge both datasets (binding and stability assays)
singles = merge(singles_binding[Mut!='*',.(Pos,WT_AA,Mut, sigma_binding = sigma,binding = fitness)],singles_stab[Mut!='*',.(Pos,WT_AA,Mut,sigma_stab= sigma,stability= fitness)],by=c("Pos","WT_AA","Mut"),all=T)
singles[,diff_binding_stability := (binding-stability)]
singles[,diff_stability_binding := (stability-binding)]
singles[,sigma_diff := sqrt((sigma_binding^2) + (sigma_stab^2))]
singles[,Mut := factor(Mut,levels=strsplit('RHKDESTNQCGPAVILMFYW','')[[1]])]

# few singles with infinite binding values, kick them out
singles[is.infinite(binding),binding := NA]
singles = singles[!is.na(binding)]
singles = singles[!is.na(stability)] 

# add p-values using t-test for the 3 replicates
ttest <- function(av, se, df = 2, mu = 1) {
  tstat <- (av-mu)/se
  # Using T-dist
  pval = 2*pt(abs(tstat), df, lower=FALSE)
  return(pval)
}
singles[,pval_bind := p.adjust(ttest(binding, sigma_binding), method = "fdr")]
singles[,pval_stab := p.adjust(ttest(stability, sigma_stab), method = "fdr")]
singles[,pval_diff := p.adjust(ttest(diff_binding_stability, sigma_diff, mu = 0), method = "fdr")]

pval_sig_thrs = 0.1


# mutations that affect binding (growth phenotype)
singles[, affects_binding := pval_bind < pval_sig_thrs & binding < 1]
# mutations that affect stability (both biochem and growth)
singles[, affects_stability := pval_stab < pval_sig_thrs & stability < 1]
# mutations effects binding (biochemical phenotype)
singles[, affects_binding_biochem := pval_bind < pval_sig_thrs & pval_diff < pval_sig_thrs & diff_binding_stability < 0]


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
ggplot(singles,aes(RSA_unbound,RSA_bound,color=HAmin_ligand < 5)) + geom_point(size=3)

## split data into core, surface and binding
ggplot(singles,aes(RSA_unbound,stability)) + geom_point() +geom_smooth()
# >>> set core threshodl as RSA <= 10
threshold_RSA = 10
threshold_ligand_dist = 5 
singles[RSA_unbound <= threshold_RSA,type := "core"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "ligand_binding"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist,type := "surface"]
singles[,type := factor(type,levels=c("core","surface","ligand_binding"))]
singles[,.(Npos = length(unique(Pos))),type]


# additional functions and thresdholds
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
# detrimental fitness threshold
threshold_fitness = 0.5

### write table of processed single mutants
write.table(x=singles,file = "dataset/GRB2_singles_dataset.txt",col.names = T,row.names=F,quote=F)



##################################################################
############################ Figure 2 ############################
##################################################################

### panel A: fitness distributions and enrichment of core/surface/binding positions in lethal mutations
panel2A_1 = ggplot(singles,aes(stability)) +
  geom_density(adjust=.5, fill="grey95")  +
  geom_vline(xintercept = threshold_fitness,lty=2) +
  coord_cartesian(xlim = range(c(singles$binding,singles$stability),na.rm=T)) +
  scale_x_continuous(breaks = seq(0,1,0.25),expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0.1)) +
  geom_text(inherit.aes = F,data=data.table(x=threshold_fitness,y=1.5,label="detrimental"),aes(x,y,label=label),hjust=1.1) +
  geom_segment(aes(x = threshold_fitness, y = 1.4, xend = threshold_fitness-0.2, yend = 1.4, colour = "segment"),arrow =  arrow(length = unit(0.03, "npc")),color="black") +
  labs(x = "fitness stabilityPCA") +
  ggtitle("stabilityPCA")
panel2A_2 = ggplot(singles,aes(binding)) +
  geom_density(adjust=.5, fill="grey95")  +
  geom_vline(xintercept = threshold_fitness,lty=2) +
  coord_cartesian(xlim = range(c(singles$binding,singles$stability),na.rm=T)) +
  scale_x_continuous(breaks = seq(0,1,0.25),expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0.1),labels = fmt_dcimals(1)) +
  geom_text(inherit.aes = F,data=data.table(x=threshold_fitness,y=2,label="detrimental"),aes(x,y,label=label),hjust=1.1) +
  geom_segment(aes(x = threshold_fitness, y = 1.8, xend = threshold_fitness-0.2, yend = 1.8, colour = "segment"),arrow =  arrow(length = unit(0.03, "npc")),color="black") +
  labs(x = "fitness bindingPCA") +
  ggtitle("bindingPCA")
panel2A_3 = ggplot(singles[,.(lethal = sum(stability < threshold_fitness,na.rm=T)/sum(!is.na(stability)) / singles[,sum(stability < threshold_fitness,na.rm=T)/sum(!is.na(stability))]),type],
                  aes(x=type,y=lethal,fill=type)) +
  geom_bar(stat="identity", width = 0.8, fill="grey95", col="black") +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0.25,0.5,1,1.5,2,3)) +
  theme(legend.position='none') +
  coord_cartesian(ylim = c(0.4,2.25)) +
  labs(x="",y="enrichment detrimental mut.",fill="") +
  ggtitle("")
panel2A_4 = ggplot(singles[,.(lethal = sum(binding < threshold_fitness,na.rm=T)/sum(!is.na(binding)) / singles[,sum(binding < threshold_fitness,na.rm=T)/sum(!is.na(binding))]),type],
                  aes(x=type,y=lethal,fill=type)) +
  geom_bar(stat="identity", width = 0.8, fill="grey95", col="black") +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0.25,0.5,1,1.5,2,3)) +
  theme(legend.position=c(0.5,0.15), legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0.4,2.25)) +
  labs(x="",y="enrichment detrimental mut.",fill="") +
  ggtitle("")

p2A = grid.arrange(grobs = list(panel2A_1,panel2A_2,panel2A_3,panel2A_4),nrow=2,
               layout_matrix = rbind(c(1,1,1,1,1,3,3),
                                     c(2,2,2,2,2,4,4)))


### panel B: comparison stability and binding fitness colouring by distance to the ligand residues

p2B = ggplot(singles,
                aes(binding, stability,color=log10(HAmin_ligand))) +
  geom_point(size=4) +
  geom_abline(linetype=2) +
  scale_color_gradient2("min dist to ligand [??]", low = col_orange, high = col_purple, mid = "grey90", midpoint = 0.9, breaks = log10(c(5,10,15,20)), labels = c(5,10,15,20)) +
  scale_x_continuous(breaks = seq(0,1,0.25),expand=c(0,0.0), limits = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0,1,0.25),expand=c(0,0.0),  limits = c(0, 1.2)) +
  labs(color="min.dist[A]",x="fitness bindingPCA",y="fitness stabilityPCA") +
  theme(legend.direction = "horizontal", legend.position = c(0.8,0.075))
p2B

### figure 1 assembled
fig2 <- grid.arrange(p2A, p2B, layout_matrix = rbind(c(1,1,1,2,2),
                                             c(1,1,1,2,2)))
ggsave("figures/fig2_assembled.pdf", fig2, width = 9, height = 3.5)


##################################################################
############################ Figure 3 ############################
##################################################################

### panel A: heatmaps of mutation effects
setkey(singles,Mut,Pos)

panel3A_1 = ggplot(singles,aes(Pos,Mut,fill=stability)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[,.(mean(stability,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,56,5)) +
  coord_cartesian(ylim=c(-0.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 1,limits=c(0,1),breaks = seq(0,1,0.5)) +
  labs(x="",y="aa mutation") +
  ggtitle("stabilityPCA")

panel3A_2 = ggplot(singles,aes(Pos,Mut,fill=binding)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[is.finite(binding),.(mean(binding,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,56,5)) +
  coord_cartesian(ylim=c(-0.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 1,limits=c(0,1),breaks = seq(0,1,0.5)) +
  labs(x="",y="aa mutation") +
  ggtitle("bindingPCA")

panel3A_3 = ggplot(singles,aes(Pos,Mut,fill=diff_binding_stability)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[is.finite(diff_binding_stability),.(mean(diff_binding_stability,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  geom_text(inherit.aes = F,data=singles[HAmin_ligand < threshold_ligand_dist,.(Pos,GAB2_AA)],aes(x=Pos,y=-2,label=GAB2_AA), size = 2) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,56,5)) +
  coord_cartesian(ylim=c(-2.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 0,limits=c(-1,1),breaks = seq(-1,1,0.5)) +
  labs(x="position",y="aa mutation",fill="bin-sta") +
  ggtitle("Difference")
## position 37 is close to ligand but doesn't affect binding more than stability; >> is only touching the C-terminal kink of the ligand

p3A = grid.arrange(grobs = list(panel3A_1,panel3A_2,panel3A_3),nrow=3)
ggsave(plot=p3A,file = "figures/fig3_panelA.pdf",height=unit(9,"cm"),width=unit(7,"cm"))


### panel B: GRB2-GAB2 crystal structure residues coloured by fitness or fitness difference
# create txt file to run in pymol

transCM = fread("dataset/PDB_contactmap_2vwf_AB.txt")                                                                              
singlesCM = merge(singles,transCM[,.(trans_dist = min(scHAmin)),Pos1],by.x="Pos",by.y="Pos1")
singlesCM_avg = singlesCM[,.(diff_binding_stability = mean(binding-stability),
                           binding = mean(binding),
                           trans_dist = unique(trans_dist), 
                           stability = mean(stability)),Pos]

script_file = "figures/pymol/Fig3_panelB.txt"
preferred_view = "set_view (  0.900416851, 0.410176218,-0.144926280,-0.081169695,-0.168878928,-0.982286036,-0.427388251, 0.896231771,-0.118769340,0.000000000, 0.000000000, -102.566505432,-2.309910297,13.775758743, 0.139800072,80.864212036,  124.268798828,  -20.000000000 )"
col_grey = "[0.8, 0.8, 0.8]"

# content in the pymol script
pymol_script = "reinitialize"
pymol_script[length(pymol_script)+1] = "fetch 2vwf, async=0"

# reference GRB2-GAB2
pymol_script[length(pymol_script)+1] = "hide everything"
pymol_script[length(pymol_script)+1] = "show cartoon, chain A"
pymol_script[length(pymol_script)+1] = paste("set cartoon_color,", col_grey)
pymol_script[length(pymol_script)+1] = "show sticks, chain B"
pymol_script[length(pymol_script)+1] = "set stick_color, orange, chain B"
pymol_script[length(pymol_script)+1] = "remove resn hoh"
pymol_script[length(pymol_script)+1] = "set stick_radius, 0.4"
pymol_script[length(pymol_script)+1] = "set ray_opaque_background, 0"
pymol_script[length(pymol_script)+1] = "set ray_shadow, 0"
pymol_script[length(pymol_script)+1] = "set ray_trace_fog, 0"
pymol_script[length(pymol_script)+1] = "set antialias, 1"
pymol_script[length(pymol_script)+1] = "bg_color white"
pymol_script[length(pymol_script)+1] = preferred_view
pymol_script[length(pymol_script)+1] = "zoom center, 25"
pymol_script[length(pymol_script)+1] = "rotate y, 45"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_reference_0.png, dpi=600")

# color GRB2 by the difference of bindingPCA-stabilityPCA fitness
for (i in 1:nrow(singlesCM_avg)) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=",singlesCM_avg[i,diff_binding_stability])
}
pymol_script[length(pymol_script)+1] = "hide labels"
pymol_script[length(pymol_script)+1] = "show sticks, chain A"
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white_blue, chain A, minimum=-0.5, maximum=0.5'
pymol_script[length(pymol_script)+1] = "set ray_opaque_background, 0"
pymol_script[length(pymol_script)+1] = "set ray_shadow, 0"
pymol_script[length(pymol_script)+1] = "set ray_trace_fog, 0"
pymol_script[length(pymol_script)+1] = "set antialias, 1"
pymol_script[length(pymol_script)+1] = "bg_color white"
pymol_script[length(pymol_script)+1] = "hide cartoon, chain A"
pymol_script[length(pymol_script)+1] = "show spheres, chain A"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_diff_spheres_0.png, dpi=600")

# color GRB2 by the bindingPCA fitness
for (i in 1:nrow(singlesCM_avg)) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=",singlesCM_avg[i,binding])
}
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white, chain A, minimum=0.5, maximum=1'
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_binding_spheres_0.png, dpi=600")



# color GRB2 by the stabilityPCA fitness
for (i in 1:nrow(singlesCM_avg)) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=",singlesCM_avg[i,stability])
}
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white, chain A, minimum=0.5, maximum=1'
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_stability_spheres_0.png, dpi=600")

# write pymol script in a .txt
write(x = pymol_script,file = script_file)


##################################################################
############################ Figure 4 ############################
##################################################################

## Aggregate fitness and delta fitness values per position. Use min value or average.
# min fitness or delta fintess values
singles_min = singles[is.finite(binding),.(diff = min(diff_binding_stability,na.rm=T),
                                           mindist = min(HAmin_ligand),
                                           stab = min(stability,na.rm=T),
                                           bind = min(binding,na.rm=T)),.(Pos,RSA_unbound)]

# average fitness or delta fitness
singles_avg = singles[,.(meandist = mean(HAmin_ligand),
                         diff = mean(diff_binding_stability,na.rm=T),
                         mindist = min(HAmin_ligand),
                         stab = mean(stability,na.rm=T),
                         dist5A = min(HAmin_ligand) < threshold_ligand_dist,
                         bind = mean(binding,na.rm=T), 
                         pval_bind = p.adjust(t.test(binding, mu=1, alternative="less")$p.val, method = "fdr"),
                         pval_stab = p.adjust(t.test(stability, mu=1, alternative="less")$p.val, method = "fdr"),
                         pval_diff = p.adjust(t.test(diff_binding_stability, mu=0, alternative="less")$p.val, method = "fdr")) ,.(Pos,RSA_unbound,RSA_bound)]


### Panel A: distribution of fitness effects vs the distance to ligand, with p-value as dot size
Y0 = melt(singles_avg,id.vars=c("Pos","mindist","RSA_unbound","RSA_bound","dist5A"))
Y_pval = Y0[grepl("pval",variable)]
Y = Y0[!grepl("pval",variable)]
Y = merge(Y,Y_pval[,.(var = strsplit(as.character(variable),"_")[[1]][2],pval=value),.(Pos,variable)][,.(Pos,variable=var,pval)],by=c("Pos","variable"))
Y[,variable := factor(variable,levels=c("stab","bind","diff"))]
levels(Y$variable) = c("stabilityPCA","bindingPCA","Difference")

panel_4A = ggplot(Y) +
  geom_point(aes(value, mindist,color=mindist < 5,shape=pval < 0.05), size=3) + 
  geom_density(aes(value,..density..,fill=mindist < 5),adjust=1,alpha=0.3) +
  # geom_hline(yintercept = -log10(pval_sig_thrs), alpha=0.5, linetype=3) +
  facet_wrap(variable ~., scales = "free_x") + 
  scale_x_continuous(expand = c(0,0.05)) +
  # scale_color_gradient2("-log10(pval)", mid = "grey90", midpoint = -log10(0.05)) +
  scale_color_manual("FDR < 0.05",values = c(col_purple,col_orange)) +
  scale_shape_manual("FDR < 0.05", values = c(21,24)) +
  scale_fill_manual(values = c(col_purple,col_orange)) +
  labs(x="average fitness per residue",y="minimal distance to GAB2",color="FDR < 0.05",fill="distance < 5Å")
panel_4A
ggsave("figures/fig4_panelA.pdf", panel_4A, width = 8.5, height = 3)

# correlation between -log10(pval) and distance
cor(-log10(singles_avg$pval_stab), log10(singles_avg$mindist))
cor(-log10(singles_avg$pval_bind), log10(singles_avg$mindist))
cor(-log10(singles_avg$pval_diff), log10(singles_avg$mindist))

# correlation between fitness and distance
cor((singles_avg$stab), log10(singles_avg$mindist))
cor((singles_avg$bind), log10(singles_avg$mindist))
cor((singles_avg$diff), log10(singles_avg$mindist))



### panel B:  ROC or precision recall curves
# predicting interaction surface residues based on the average or the minimum fitness value per position

# using min fitness or delta fitness per position
pr_stab = pr.curve(singles_min[mindist < threshold_ligand_dist,-stab],
                   singles_min[mindist > threshold_ligand_dist,-stab],curve=T)
pr_bind = pr.curve(singles_min[mindist < threshold_ligand_dist,-bind],
                   singles_min[mindist > threshold_ligand_dist,-bind],curve=T)
pr_diff_bs = pr.curve(singles_min[mindist < threshold_ligand_dist,-diff],
                      singles_min[mindist > threshold_ligand_dist,-diff],curve=T)
# using average
pr_stab_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-stab],
                       singles_avg[mindist > threshold_ligand_dist,-stab],curve=T)
pr_bind_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-bind],
                       singles_avg[mindist > threshold_ligand_dist,-bind],curve=T)
pr_diff_bs_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-diff],
                          singles_avg[mindist > threshold_ligand_dist,-diff],curve=T)

# merge both metrics to show in a single plot
pr_all = rbind(data.table(metric = "min", fitness = "stability",recall = rev(pr_stab$curve[,1]),precision = rev(pr_stab$curve[,2])),
               data.table(metric = "min", fitness = "binding",recall = rev(pr_bind$curve[,1]),precision = rev(pr_bind$curve[,2])),
               data.table(metric = "min", fitness = "difference",recall = rev(pr_diff_bs$curve[,1]),precision = rev(pr_diff_bs$curve[,2])),
               data.table(metric = "avg", fitness = "stability",recall = rev(pr_stab_avg$curve[,1]),precision = rev(pr_stab_avg$curve[,2])),
               data.table(metric = "avg", fitness = "binding",recall = rev(pr_bind_avg$curve[,1]),precision = rev(pr_bind_avg$curve[,2])),
               data.table(metric = "avg", fitness = "difference",recall = rev(pr_diff_bs_avg$curve[,1]),precision = rev(pr_diff_bs_avg$curve[,2])))
pr_all[,fitness := factor(fitness,levels=c("stability","binding","difference"))]
pr_all[,metric := factor(metric, levels = c("min", "avg"))]

panel_4B_1 = ggplot(pr_all,aes(recall,precision,color=fitness)) +
  geom_line(size=1, aes(linetype=metric), show.legend = F) + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0.01)) +
  scale_color_manual(values = c(col_purple, col_orange , "grey50"))
panel_4B_1

pr_all_AUC = rbind(data.table(metric = "min", fitness = "stability", AUC=round(pr_stab$auc.integral,2)),
               data.table(metric = "min", fitness = "binding", AUC=round(pr_bind$auc.integral,2)),
               data.table(metric = "min", fitness = "difference", AUC=round(pr_diff_bs$auc.integral,2)),
               data.table(metric = "avg", fitness = "stability", AUC=round(pr_stab_avg$auc.integral,2)),
               data.table(metric = "avg", fitness = "binding", AUC=round(pr_bind_avg$auc.integral,2)),
               data.table(metric = "avg", fitness = "difference", AUC=round(pr_diff_bs_avg$auc.integral,2)))
pr_all_AUC[,fitness := factor(fitness,levels=c("stability","binding","difference"))]
pr_all_AUC[,metric := factor(metric, levels = c("avg", "min"))]
panel_4B_2 = ggplot(pr_all_AUC, aes(x=fitness, y=metric, fill=AUC)) +
  geom_tile(color="black", size=0.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylab("metric per position used") + xlab("fitness score used") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0,1)) +
  geom_text(aes(label=AUC), color="white", size=5)
panel_4B_2 

# assemble panel B
panel_4B = grid.arrange(panel_4B_1, panel_4B_2, layout_matrix = rbind(c(1,1),
                                                                   c(1,1),
                                                                   c(2,2)))
ggsave("figures/fig4_panelB.pdf", panel_4B, width = 3.5, height = 4.5)





##################################################################
########################## Text numbers ##########################
##################################################################
# Number of mutations altering stability, binding only, both or none.

## Per mutation
singles[is.finite(binding) & is.finite(stability),.N] #total of 504 mutations
# bindingPCA vs stabilityPCA
singles[is.finite(binding) & is.finite(stability),.(total=.N,freq = .N/singles[is.finite(binding) & is.finite(stability),.N]),.(affects_binding,affects_stability)]

# Binding affinity only using delta fitness score
singles[is.finite(binding) & is.finite(stability),.(total=.N,freq = .N/singles[is.finite(binding) & is.finite(stability),.N]),.(affects_binding_biochem,affects_stability)]


## Per position
# mutations that affect binding (growth phenotype)
singles_avg[, affects_binding := pval_bind < pval_sig_thrs & bind < 1]
# mutations that affect stability (both biochem and growth)
singles_avg[, affects_stability := pval_stab < pval_sig_thrs & stab < 1]
# mutations effects binding (biochemical phenotype)
singles_avg[, affects_binding_biochem := pval_bind < pval_sig_thrs & pval_diff < pval_sig_thrs & diff < 0]

# Number of positions
singles_avg[is.finite(bind) & is.finite(stab),.N] #total of 56 positions
# bindingPCA vs stabilityPCA
singles_avg[is.finite(bind) & is.finite(stab),.(total=.N,freq = .N/singles_avg[is.finite(bind) & is.finite(stab),.N]),.(affects_binding,affects_stability)]

# Binding affinity only using average delta fitness score
singles_avg[is.finite(bind) & is.finite(stab),.(total=.N,freq = .N/singles_avg[is.finite(bind) & is.finite(stab),.N]),.(affects_binding_biochem,affects_stability)]




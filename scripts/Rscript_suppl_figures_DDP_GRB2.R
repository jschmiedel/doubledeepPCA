##### double deep supplementary figure plots ######

#### version: 1.0
#### last modified: 2019/11/11
#### created by: Julia
#### last modified by: Julia
# put all supplementary figs together

library(gdata)
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_classic(base_size=9))
require(data.table)
require(GGally)
require(PRROC)

#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")


## Colors
col_purple = "#9161A8"
col_blue =  "#0066CC" # "#4990C9"
col_orange = "#F7941E"
col_red = "#EF4136"


### load GRB2 single mutation data
# use  source data with read count threshold 20 and non-logged fitness values
singles_binding = fread("dataset/GRB2_CYC_singles_readT20.txt")
singles_binding[, input_mean_reads := rowMeans(singles_binding[,.(input1_e1_s0_bNA_count_agg, input2_e2_s0_bNA_count_agg, input3_e3_s0_bNA_count_agg)])]
singles_stab = fread("dataset/GRB2_GPD_singles_readT20.txt")
singles_stab[, input_mean_reads := rowMeans(singles_stab[,.(input1_e1_s0_bNA_count_agg, input2_e2_s0_bNA_count_agg, input3_e3_s0_bNA_count_agg)])]


# merge both datasets (binding and stability assays)
singles = merge(singles_binding[Mut!='*',.(Pos,WT_AA,Mut, sigma_binding = sigma,binding = fitness, input_mean_reads_binding = input_mean_reads)],
                singles_stab[Mut!='*',.(Pos,WT_AA,Mut,sigma_stab= sigma, stability= fitness, input_mean_reads_stab = input_mean_reads)], by=c("Pos","WT_AA","Mut"),all=T)

singles[,diff_binding_stability := (binding-stability)]
singles[,diff_stability_binding := (stability-binding)]
singles[, input_mean_reads := rowMeans(singles[, .(input_mean_reads_binding, input_mean_reads_stab)])]
singles[,sigma_diff := sqrt((sigma_binding^2) + (sigma_stab^2))]
singles[,Mut := factor(Mut,levels=strsplit('RHKDESTNQCGPAVILMFYW','')[[1]])]
#few singles with infinite binding values, kick them out
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
#add to singles data.table
singles = merge(singles,RSA_chainA[chain == "A",.(Pos,RSA_unbound = RSA_all_rel)],by="Pos")
singles = merge(singles,RSA_chainAB[chain == "A",.(Pos,RSA_bound = RSA_all_rel)],by="Pos")

threshold_RSA = 10
threshold_ligand_dist = 5 
singles[RSA_unbound <= threshold_RSA,type := "core"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "ligand_binding"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist,type := "surface"]
singles[,type := factor(type,levels=c("core","surface","ligand_binding"))]
singles[,.(Npos = length(unique(Pos))),type]

################################################################################
############################ Supplementary Figure 1 ############################
################################################################################
# Growth rate of WTs, controls and individual single mutant validations 

### load data
# OD600 measurements with time
a1 = read.delim("dataset/GRB2_mutations_MTXassay1.txt") # bindingPCA data
a1 <- a1[grepl("[BCD]", a1$row),] # Selecting rows which competition media has 200 ug/ml MTX
a2 = read.delim("dataset/GRB2_mutations_MTXassay2.txt") # stabilityPCA data
# mapping from plasmid id to construct information
pl2id_1 = read.xls("dataset/plasmid2id_MTXassay1.xlsx",sheet = 1, header = TRUE)
pl2id_1 <- pl2id_1[, 1:4]
pl2id_2 = read.xls("dataset/plasmid2id_MTXassay2.xlsx",sheet = 1, header = TRUE)
pl2id_2 <- pl2id_2[, 1:3]

## Calculate real ODs substracting the blank for each of the two plate assays separatetly
dtp_1 = 5
od1 = as.matrix(a1[,-(1:dtp_1)])
real_od1 = do.call("rbind",apply(od1, 1, function(x){
  x - a1[a1$plasmid == "blank", -(1:dtp_1)]
}))

dtp_2 = 4
od2 = as.matrix(a2[,-(1:dtp_2)])
blanks <- colMeans(a2[a2$plasmid == "blank", -(1:dtp_2)])
real_od2 = as.data.frame(t(apply(od2, 1, function(x){
  x - blanks
})))

# Create file with real ODs and plasmid name
ds1 = data.frame(real_od1)
ds1$plasmid = a1$plasmid
ds1$well = paste(a1$row, a1$column, sep = "")
ds1 = ds1[!(ds1$plasmid == "blank"),]
DS1 = merge(ds1, pl2id_1, by.x = "plasmid", by.y="plasmid_short")
DS_melt1 = melt(DS1, id.vars = c("plasmid_long","domain", "assay", "well", "plasmid"))
DS_melt1$time = sapply(as.character(DS_melt1$variable), function(x){as.numeric(substr(x,2,nchar(x)))/3600})
DS_melt1 = DS_melt1[, c("plasmid_long","domain", "assay", "well","variable", "value", "time")]
colnames(DS_melt1) <- c("plasmid","domain", "assay", "well","variable", "value", "time")

ds2 = data.frame(real_od2)
ds2$plasmid = a2$plasmid
ds2$well = paste(a2$row, a2$column, sep = "")
ds2 = ds2[!(ds2$plasmid == "blank"),]
DS2 = merge(ds2, pl2id_2, by = "plasmid")
DS_melt2 = melt(DS2, id.vars = c("plasmid","domain", "assay", "well"))
DS_melt2$time = sapply(as.character(DS_melt2$variable), function(x){as.numeric(substr(x,2,nchar(x)))/3600})

## Calculate Growth rate (slope of linear fit on exponential phase)
gr1 <- do.call("rbind", lapply(unique(DS_melt1$well), function(x){
  ds = DS_melt1[DS_melt1$well == x,]
  mid_od = max(ds$value)/2
  time_od50 = max(ds$time[ds$value <= mid_od])
  lmds = data.frame(t = ds$time[ds$time <= time_od50+2 & ds$time >= time_od50-2])
  lmds$ods = log(ds$value[ds$time %in% lmds$t])
  lmfit = lm(formula = ods ~ t, data = lmds)
  data.frame(plasmid = unique(ds$plasmid), well = x, domain = unique(ds$domain), assay = unique(ds$assay), 
             start_od = min(ds$value), mid_od = mid_od, start_time = min(ds$time), time_od50 = time_od50, growth_rate_slope = lmfit$coefficients[2])
}))

gr2 <- do.call("rbind", lapply(unique(DS_melt2$well), function(x){
  ds = DS_melt2[DS_melt2$well == x & DS_melt2$time <= 50,]
  mid_od = max(ds$value)/2
  time_od50 = max(ds$time[ds$value <= mid_od])
  lmds = data.frame(t = ds$time[ds$time <= time_od50+2 & ds$time >= time_od50-2])
  lmds$ods = log(ds$value[ds$time %in% lmds$t])
  lmfit = lm(formula = ods ~ t, data = lmds)
  data.frame(plasmid = unique(ds$plasmid), well = x, domain = unique(ds$domain), assay = unique(ds$assay), 
             start_od = min(ds$value), mid_od = mid_od, start_time = min(ds$time), time_od50 = time_od50, growth_rate_slope = lmfit$coefficients[2])
}))

gr <- rbind(gr1[!gr1$domain == "Neg control 2",], gr2)
gr$domain <- factor(gr$domain, levels = c("Pos control", "GRB2", "M46W", "M46L", "F7H", "Y51R", "P48H", "Y51S", "Y51D", "Neg control"))


### Panel A: plot growth rates of WTs, controls and GRB2 single mutants for both assays all together

dodge <- position_dodge(width = 0.7)
panel_1A <- ggplot(gr, aes(x=domain, y=growth_rate_slope, fill = assay, color=assay)) + 
  geom_boxplot(width = 0.7, position = dodge) + geom_point(aes(fill=assay),position = position_jitterdodge(), alpha = 0.6, show.legend = F) + 
  xlab("variant") + ylab("Growth rate") +
  # scale_fill_manual("", values = c("#c8b0d3", "#fbc98e")) +
  # scale_color_manual("", values = c(col_purple, col_orange)) +
  scale_fill_manual("", values = c("#fbc98e","#c8b0d3")) +
  scale_color_manual("", values = c(col_orange, col_purple)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.direction = "horizontal", legend.position = "top")
panel_1A
ggsave(plot=panel_1A,file = "figures/Supp_fig_1A.pdf",height=unit(4,"cm"),width=unit(7,"cm"))


### panel B: comparison of WT relative growth rates to fitness values from sequencing data
# Get the average growth rate per GRB2 mutant in each assay 
gr_mean <- aggregate(growth_rate_slope ~ domain+assay, data=gr[, c("growth_rate_slope", "domain", "assay")], FUN = mean)
gr_mean_stab <- gr_mean[gr_mean$assay == "stabilityPCA",]
gr_mean_stab$gr_rel2wt <-  gr_mean_stab$growth_rate_slope-gr_mean_stab$growth_rate_slope[gr_mean_stab$domain == "GRB2"]
gr_mean_binding <- gr_mean[gr_mean$assay == "bindingPCA",]
gr_mean_binding$gr_rel2wt <-  gr_mean_binding$growth_rate_slope-gr_mean_binding$growth_rate_slope[gr_mean_binding$domain == "GRB2"]

# merge it with the sequencing data
singles_binding[, id := paste(WT_AA, Pos, Mut, sep="")]
singles_stab[, id := paste(WT_AA, Pos, Mut, sep="")]

GRvsSeq <- rbind(merge(gr_mean_binding, singles_binding, by.x="domain", by.y = "id"), merge(gr_mean_stab, singles_stab, by.x="domain", by.y = "id"))

panel_1B <- ggplot(GRvsSeq, aes(x=gr_rel2wt, y=fitness)) + geom_point(aes(color=assay)) +
  geom_text(aes(label=domain, color=assay, x=gr_rel2wt+0.007, y=fitness-0.04), show.legend = F) +
  theme(legend.position = c(0.15, 0.9)) +
  scale_color_manual(values = c(col_orange, col_purple)) +
  annotate("text", x=-0.11, y=0.8, label=paste0("r = ", round(cor(GRvsSeq$fitness, GRvsSeq$gr_rel2wt, method = "spearman"),2))) +
  xlab("growth rate relative to WT (individual measurements)") + ylab("fitness relative to WT (deep sequencing)")
panel_1B
ggsave(plot=panel_1B,file = "figures/Supp_fig_1B.pdf",height=unit(4.5,"cm"),width=unit(5,"cm"))
# p-value of this last correlation
cor.test(GRvsSeq$fitness, GRvsSeq$gr_rel2wt, method = "spearman") # p-value < 2.2e-16


################################################################################
############################ Supplementary Figure 2 ############################
################################################################################
# Mutational frequency and bias of the error prone PCR library

### Load dataset
# Raw DiMSum dataset (not the already processed ones)
load("dataset/DiMSum_GRB2_CYC_variant_data_merge_Q25.RData")
binding <- variant_data_merge[Ndel_nt == 0 & Nins_nt == 0,]
rm(variant_data_merge)

load("dataset/DiMSum_GRB2_GPD_variant_data_merge_Q25.RData")
stability <- variant_data_merge[Ndel_nt == 0 & Nins_nt == 0,]
rm(variant_data_merge)

stability_bynumnt <- data.frame(num_nt_mut = 0:12,
                                rep1 = aggregate(input1_e1_s0_bNA_count ~ Nmut_nt, data = stability, sum)$input1_e1_s0_bNA_count*100/sum(stability$input1_e1_s0_bNA_count),
                                rep2 = aggregate(input2_e2_s0_bNA_count ~ Nmut_nt, data = stability, sum)$input2_e2_s0_bNA_count*100/sum(stability$input2_e2_s0_bNA_count),
                                rep3 = aggregate(input3_e3_s0_bNA_count ~ Nmut_nt, data = stability, sum)$input3_e3_s0_bNA_count*100/sum(stability$input3_e3_s0_bNA_count))
stability_bynumnt_melt <- melt(stability_bynumnt, id.vars = "num_nt_mut")

binding_bynumnt <- data.frame(num_nt_mut = 0:12,
                              rep1 = aggregate(input1_e1_s0_bNA_count ~ Nmut_nt, data = binding, sum)$input1_e1_s0_bNA_count*100/sum(binding$input1_e1_s0_bNA_count),
                              rep2 = aggregate(input2_e2_s0_bNA_count ~ Nmut_nt, data = binding, sum)$input2_e2_s0_bNA_count*100/sum(binding$input2_e2_s0_bNA_count),
                              rep3 = aggregate(input3_e3_s0_bNA_count ~ Nmut_nt, data = binding, sum)$input3_e3_s0_bNA_count*100/sum(binding$input3_e3_s0_bNA_count))
binding_bynumnt_melt <- melt(binding_bynumnt, id.vars = "num_nt_mut")


### panel A: proportion of reads per nt mutation ######

panel_2A1 <- ggplot(stability_bynumnt_melt[stability_bynumnt_melt$num_nt_mut <=6,], aes(x = num_nt_mut, y=value, color=variable)) + 
  geom_point(size=2) +  geom_line(alpha=0.5) + scale_color_manual("", values = c(col_purple, col_red, col_orange)) +
  xlab("Number of nucleotide mutations") + ylab("Percentage of variants read counts") +
  scale_x_continuous(breaks = seq(0,6,1)) + theme(legend.position = c(0.8, 0.7))+
  ggtitle("stabilityPCA")
panel_2A1

panel_2A2 <- ggplot(binding_bynumnt_melt[binding_bynumnt_melt$num_nt_mut <=6,], aes(x = num_nt_mut, y=value, color=variable)) + 
  geom_point(size=2) +  geom_line(alpha=0.5) + scale_color_manual("", values = c(col_purple, col_red, col_orange)) +
  xlab("Number of nucleotide mutations") + ylab("Percentage of variants read counts") +
  scale_x_continuous(breaks = seq(0,6,1)) + theme(legend.position = c(0.8, 0.7)) +
  ggtitle("bindingPCA")
panel_2A2

panel_2A <- grid.arrange(panel_2A1, panel_2A2, nrow=2)
ggsave(plot=panel_2A,file = "figures/Supp_fig_2A.pdf",height=unit(4,"cm"),width=unit(5,"cm"))


##### panel B: mutational bias of the epPCR ######

wildtype = binding[Nmut_nt == 0,]
wt_ntseq_split = strsplit(wildtype[WT==T,nt_seq],"")[[1]][1:171]
wt_ntseq = paste0(wt_ntseq_split,collapse = "")

binding_1 = binding[Nmut_nt == 1,]
binding_1[,Pos := which(strsplit(nt_seq,"")[[1]][1:171] !=wt_ntseq_split),nt_seq]
binding_1[,Mut := toupper(strsplit(nt_seq,"")[[1]][Pos]),nt_seq]
binding_1[,WT_nt := toupper(wt_ntseq_split[Pos]),nt_seq]
binding_1[, id := paste(WT_nt, Pos, Mut, sep = "")]

stability_1 = stability[Nmut_nt == 1,]
stability_1[,Pos := which(strsplit(nt_seq,"")[[1]][1:171] !=wt_ntseq_split),nt_seq]
stability_1[,Mut := toupper(strsplit(nt_seq,"")[[1]][Pos]),nt_seq]
stability_1[,WT_nt := toupper(wt_ntseq_split[Pos]),nt_seq]
stability_1[, id := paste(WT_nt, Pos, Mut, sep = "")]

# Categorize mutations given mutational transition or transversion
list_regex = c("(A[0-9]+G)|(T[0-9]+C)", "(G[0-9]+A)|(C[0-9]+T)", "(A[0-9]+T)|(T[0-9]+A)", "(A[0-9]+C)|(T[0-9]+G)", "(G[0-9]+C)|(C[0-9]+G)", "(G[0-9]+T)|(C[0-9]+A)")
list_transi_transv = c("A>G", "G>A", "A>T", "A>C", "G>C", "G>T")

binding_1$bias = do.call("c", lapply(binding_1$id, function(x){
  list_transi_transv[do.call("c", lapply(1:length(list_regex), function(y){grepl(list_regex[y], x)}))]
}))
stability_1$bias = do.call("c", lapply(stability_1$id, function(x){
  list_transi_transv[do.call("c", lapply(1:length(list_regex), function(y){grepl(list_regex[y], x)}))]
}))

# plot reads distribution of each mutational bias for each of the replicates
reps = c("input1_e1_s0_bNA_count", "input2_e2_s0_bNA_count", "input3_e3_s0_bNA_count")

pL2B1 <- lapply(1:3, function(x){
  ds = as.data.frame(stability_1)
  DS = setNames(ds[, c(reps[x], "bias")], c("reads", "bias"))
  ggplot(DS, aes(x=reads, color=bias)) + 
    geom_density(aes(fill=bias), alpha=0.3) +
    scale_x_log10() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    xlab(paste0("number of input read counts (replicate ", x, ")")) +
    theme(legend.position = c(0.9, 0.8))
})
panel_2B1 <- do.call("grid.arrange", c(pL2B1, top="stabilityPCA"))

pL2B2 <- lapply(1:3, function(x){
  ds = as.data.frame(binding_1)
  DS = setNames(ds[, c(reps[x], "bias")], c("reads", "bias"))
  ggplot(DS, aes(x=reads, color=bias)) + 
    geom_density(aes(fill=bias), alpha=0.3, show.legend = F) +
    scale_x_log10() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    xlab(paste0("number of input read counts (replicate ", x, ")"))
})
panel_2B2 <- do.call("grid.arrange", c(pL2B2, top="bindingPCA"))

panel_2B = grid.arrange(panel_2B1, panel_2B2, nrow=1)
ggsave(plot=panel_2B,file = "figures/Supp_fig_2B.pdf",height=unit(5.5,"cm"),width=unit(9,"cm"))


#### NUMBERS #####
# Number of average nt mutations in GRB2
mean(do.call("c", lapply(1:nrow(stability_bynumnt_melt), function(x){
  rep(stability_bynumnt_melt$num_nt_mut[x], round(stability_bynumnt_melt$value[x]))
}))) 
# 1.63

mean(do.call("c", lapply(1:nrow(binding_bynumnt_melt), function(x){
  rep(binding_bynumnt_melt$num_nt_mut[x], round(binding_bynumnt_melt$value[x]))
})))
# 1.65



################################################################################
############################ Supplementary Figure 3 ############################
################################################################################
# Correalation between replicates

### Load doubles
doubles_binding = fread("dataset/GRB2_CYC_doubles_readT20.txt")
doubles_stab = fread("dataset/GRB2_GPD_doubles_readT20.txt")

### Functions for correlation plots
lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point(alpha=0.5)+
    # scale_x_continuous(limits = c(-1.81,0.15))+
    # scale_y_continuous(limits = c(-1.81,0.15)) +
    scale_x_continuous(breaks = seq(0,1,.5))+
    scale_y_continuous(breaks = seq(0,1,.5)) +
    geom_density2d(color="grey80", alpha=0.7) +
    geom_abline(linetype=2, color="grey50")
}  
diagfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_density(fill="grey80") +
    # scale_x_continuous(limits = c(-0.1,1.25))
    scale_x_continuous(breaks = seq(0,1,.5))
  
}

lowerfun_d <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    stat_binhex() + scale_fill_gradient(low = "grey80", high = "black") +
    # scale_x_continuous(limits = c(-2.41,0.27))+
    # scale_y_continuous(limits = c(-2.41,0.27)) +
    scale_x_continuous(breaks = seq(0,1,.5))+
    scale_y_continuous(breaks = seq(0,1,.5)) +
    geom_abline(linetype=2, color="grey50")
}  
diagfun_d <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_density(fill="grey80") +
    # scale_x_continuous(limits = c(-2.41,0.27))
    scale_x_continuous(breaks = seq(0,1,.5))
}


##### panel A: Correlations singles stability  
panel_3A <- ggpairs(singles_stab[is.finite(fitness1) & is.finite(fitness2) & is.finite(fitness3), 
                                c("fitness1", "fitness2", "fitness3")], 
                   lower = list(continuous = wrap(lowerfun)), 
                   diag = list(continuous = wrap(diagfun))) +
  theme(legend.position = "none", panel.grid.major = element_blank()) + ggtitle("Singles, stabilityPCA")
ggsave("figures/Supp_fig_3A.pdf", panel_3A, width = 5, height = 4.5)


##### panel B: Correlations singles binding
panel_3B <- ggpairs(singles_binding[is.finite(fitness1) & is.finite(fitness2) & is.finite(fitness3),
                                   c("fitness1", "fitness2", "fitness3")], 
                   lower = list(continuous = wrap(lowerfun)), 
                   diag = list(continuous = wrap(diagfun))) +
  theme(legend.position = "none", panel.grid.major = element_blank())  + ggtitle("Singles, bindingPCA")
ggsave("figures/Supp_fig_3B.pdf", panel_3B, width = 5, height = 4.5)


##### panel C: Correlations doubles stability  
panel_3C <- ggpairs(doubles_stab[is.finite(fitness1_uncorr) & is.finite(fitness2_uncorr) & is.finite(fitness3_uncorr),
                                c("fitness1_uncorr", "fitness2_uncorr", "fitness3_uncorr")], 
                   lower = list(continuous = wrap(lowerfun_d)), 
                   diag = list(continuous = wrap(diagfun_d)))  +
  theme(legend.position = "none", panel.grid.major = element_blank()) + ggtitle("Doubles, stabilityPCA")
ggsave("figures/Supp_fig_3C.pdf", panel_3C, width = 5, height = 4.5)



##### panel D: Correlations doubles binding
panel_3D <- ggpairs(doubles_binding[is.finite(fitness1_uncorr) & is.finite(fitness2_uncorr) & is.finite(fitness3_uncorr), 
                                   c("fitness1_uncorr", "fitness2_uncorr", "fitness3_uncorr")], 
                   lower = list(continuous = wrap(lowerfun_d)), 
                   diag = list(continuous = wrap(diagfun_d)))  +
  theme(legend.position = "none", panel.grid.major = element_blank()) + ggtitle("Doubles, bindingPCA")
ggsave("figures/Supp_fig_3D.pdf", panel_3D, width = 5, height = 4.5)




### NUMBERS ###

# bindingPCA single mutants obtained
singles_binding[is.finite(fitness) & !STOP,.N] #599
# stabilityPCA single mutants obtained
singles_stab[is.finite(fitness) & !STOP,.N] #527
# number of shared singles 
nrow(merge(singles_stab[is.finite(fitness) & !STOP], singles_binding[is.finite(fitness) & !STOP], by=c("Mut", "Pos"), all = F)) #504

# median number of aa substitutions per position 
median(table(singles_binding[is.finite(fitness) & !STOP, Pos])) #6
median(table(singles_stab[is.finite(fitness) & !STOP, Pos])) #6
median(table(singles[, Pos])) #6
# min number of aa substitutions per position 
min(table(singles_binding[is.finite(fitness) & !STOP, Pos])) #4
min(table(singles_stab[is.finite(fitness) & !STOP, Pos])) #4
min(table(singles[, Pos])) #6
# max number of aa substitutions per position 
max(table(singles_binding[is.finite(fitness) & !STOP, Pos])) #7
max(table(singles_stab[is.finite(fitness) & !STOP, Pos])) #7
max(table(singles[, Pos])) #6

# % of detrimental mutations in both assays
threshold_fitness = 0.5
singles_binding[is.finite(fitness) & !STOP,.(below_threshold = sum(fitness < threshold_fitness),total=.N,frac = sum(fitness < threshold_fitness)/.N)]
# 24.8%
singles_stab[is.finite(fitness) & !STOP,.(below_threshold = sum(fitness < threshold_fitness),total=.N,frac = sum(fitness < threshold_fitness)/.N)]#)*100/nrow(singles_binding) #
# 27.13%

## doubles
doubles_binding[!STOP & !is.na(fitness) & is.finite(fitness),.N]
#24952
doubles_stab[!STOP & !is.na(fitness) & is.finite(fitness),.N]
#19985
merge(doubles_binding[!STOP & !is.na(fitness) & is.finite(fitness),],
      doubles_stab[!STOP & !is.na(fitness) & is.finite(fitness),],by=c("Pos1","Pos2","Mut1","Mut2"),all=F)[,.N]
#17464 shared



################################################################################
############################ Supplementary Figure 4 ############################
################################################################################
# Comparison of doubles's fitness values in bindingPCA and stabilityPCA 

# merge doubles
doubles = merge(doubles_binding[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,binding = fitness)], doubles_stab[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,stability = fitness)], by=c("Mut1", "Pos1", "Mut2", "Pos2", "WT_AA1","WT_AA2"), all = F)

### load distances between GRB2 and ligand
GRB2_GAB2_distances = fread("dataset/PDB_contactmap_2vwf_AB.txt")
GRB2_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]

temp = merge(doubles,GRB2_GAB2_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1],by.x="Pos1",by.y="Pos1")
doubles = merge(temp,GRB2_GAB2_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1],by.x="Pos2",by.y="Pos1", suffixes = c("_Pos1", "_Pos2"))
doubles_av = doubles[,.(stab = mean(stability,na.rm=T),
                        bind = mean(binding,na.rm=T)),.(Pos1,Pos2,HAmin_ligand_Pos1, HAmin_ligand_Pos2)]


panel_4A = ggplot(doubles_av[, .(HAmin_ligand_mean = mean(HAmin_ligand_Pos1, HAmin_ligand_Pos2)), .(bind, stab)],
                 aes(bind, stab,color=log10(HAmin_ligand_mean))) +
  geom_point(size=3) +
  geom_abline(linetype=2) +
  scale_color_gradient2("average min dist [??]", low = col_orange, high = col_purple, mid = "grey90", midpoint = 0.9, breaks = log10(c(5,10,15)), labels = c(5,10,15)) +
  labs(color="min.dist[Å]",x="fitness bindingPCA",y="fitness stabilityPCA") +
  theme(legend.direction = "vertical")
panel_4A
ggsave(plot=panel_4A,file = "figures/Supp_fig_4A.pdf",height=unit(5,"cm"),width=unit(7,"cm"))



################################################################################
############################ Supplementary Figure 5 ############################
################################################################################
# Figure rotating the structures of GRB2-GAB2 with each residue coloured by the average fitness

transCM = fread("dataset/PDB_contactmap_2vwf_AB.txt")                                                                              
singlesCM = merge(singles,transCM[,.(trans_dist = min(scHAmin)),Pos1],by.x="Pos",by.y="Pos1")
singlesCM_avg = singlesCM[,.(diff_binding_stability = mean(binding-stability),
                             binding = mean(binding),
                             trans_dist = unique(trans_dist), 
                             stability = mean(stability)),Pos]

script_file = "figures/pymol/Suppl_fig5.txt"
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
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_reference_1.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, 90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_reference_2.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, -90"
pymol_script[length(pymol_script)+1] = "rotate x, -90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_reference_3.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate x, 90"


# color GRB2 by the difference of deepPCA-stabilityPCA fitness
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
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_diff_spheres_1.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, 90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_diff_spheres_2.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, -90"
pymol_script[length(pymol_script)+1] = "rotate x, -90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_diff_spheres_3.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate x, 90"

# color GRB2 by the deepPCA fitness
for (i in 1:nrow(singlesCM_avg)) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=",singlesCM_avg[i,binding])
}
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white, chain A, minimum=0.5, maximum=1'
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_binding_spheres_1.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, 90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_binding_spheres_2.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, -90"
pymol_script[length(pymol_script)+1] = "rotate x, -90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_binding_spheres_3.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate x, 90"


# color GRB2 by the stabilityPCA fitness
for (i in 1:nrow(singlesCM_avg)) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=",singlesCM_avg[i,stability])
}
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white, chain A, minimum=0.5, maximum=1'
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_stability_spheres_1.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, 90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_stability_spheres_2.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate y, -90"
pymol_script[length(pymol_script)+1] = "rotate x, -90"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 001-GRB2-GAB2_stability_spheres_3.png, dpi=600")
pymol_script[length(pymol_script)+1] = "rotate x, 90"

# write pymol script in a .txt
write(x = pymol_script,file = script_file)




################################################################################
############################ Supplementary Figure X ############################
################################################################################

ggplot(singles, aes(x=input_mean_reads_binding, y=input_mean_reads_stab)) + geom_point() + geom_abline(col="red")
ggplot(singles, aes(x=(input_mean_reads))) + geom_histogram(bins = 50)

read_thrs = c(15000,10000,7500,5000,2500,1000,500,250,100,50,25,10)
AUC_thrs_reads = do.call("rbind",lapply(read_thrs, function(x){
  singles_avg = singles[input_mean_reads > x,.(meandist = mean(HAmin_ligand),
                           diff = mean(diff_binding_stability,na.rm=T),
                           diff_min = min(diff_binding_stability, na.rm = T),
                           stab = mean(stability,na.rm=T),
                           stab_min = min(stability,na.rm=T),
                           bind = mean(binding,na.rm=T), 
                           bind_min = min(binding,na.rm=T), 
                           mindist = min(HAmin_ligand),
                           dist5A = min(HAmin_ligand) < threshold_ligand_dist,
                           aa_num = length(unique(Mut))), .(Pos,RSA_unbound,RSA_bound)]
  
  pr_stab = pr.curve(singles_avg[mindist < threshold_ligand_dist,-stab_min],
                     singles_avg[mindist > threshold_ligand_dist,-stab_min],curve=T)
  pr_bind = pr.curve(singles_avg[mindist < threshold_ligand_dist,-bind_min],
                     singles_avg[mindist > threshold_ligand_dist,-bind_min],curve=T)
  pr_diff_bs = pr.curve(singles_avg[mindist < threshold_ligand_dist,-diff_min],
                        singles_avg[mindist > threshold_ligand_dist,-diff_min],curve=T)
  pr_stab_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-stab],
                     singles_avg[mindist > threshold_ligand_dist,-stab],curve=T)
  pr_bind_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-bind],
                     singles_avg[mindist > threshold_ligand_dist,-bind],curve=T)
  pr_diff_bs_avg = pr.curve(singles_avg[mindist < threshold_ligand_dist,-diff],
                        singles_avg[mindist > threshold_ligand_dist,-diff],curve=T)
  
  data.table(AUC = c( round(pr_stab_avg$auc.integral,2), round(pr_bind_avg$auc.integral,2), round(pr_diff_bs_avg$auc.integral,2),
                      round(pr_stab$auc.integral,2), round(pr_bind$auc.integral,2), round(pr_diff_bs$auc.integral,2) ),
             metric = c("avg", "avg", "avg", "min", "min", "min"),
             fitness = c("stability", "binding", "difference","stability", "binding", "difference"),
             read_thrs = rep(x, 6),
             median_aa_num = rep(median(singles_avg$aa_num), 6),
             mean_aa_num = rep(median(singles_avg$aa_num), 6),
             total_num_aa = rep(sum(singles_avg$aa_num)))
}))

AUC_thrs_reads[,metric := factor(metric, levels = c("min", "avg"))]
AUC_thrs_reads[, pct_total_num_aa := round(total_num_aa*100/1064, 1)]

panel_6A <- ggplot(AUC_thrs_reads, aes(x=read_thrs, y=AUC, color=fitness)) + geom_line(aes(linetype=metric)) + 
  geom_point() +
  scale_color_manual(values = c( col_orange , "grey50",col_purple)) +
  scale_x_log10() +
  xlab("Minimum number of read counts") +
  geom_text(data=AUC_thrs_reads[metric == "min" & fitness== "difference",], aes(label=paste0(pct_total_num_aa, "%"), x=read_thrs, y=AUC+0.04), show.legend = F) +
  geom_text(data=AUC_thrs_reads[metric == "avg" & fitness== "difference",], aes(label=paste0(median_aa_num, " aa"), x=read_thrs, y=AUC-0.03), show.legend = F)
panel_6A

ggsave(plot=panel_6A,file = "figures/Supp_fig_6A.pdf",height=unit(5,"cm"),width=unit(8,"cm"))




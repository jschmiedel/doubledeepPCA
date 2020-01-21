### import and process PDZ data 

## source functions
filelist = list.files('functions/')
sapply(paste0('functions/',filelist),source,.GlobalEnv)

# create directory structure
dir.create("processed_data/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/preprocessing/", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
theme_set(theme_classic(base_size=9))

### load the different data files
## wildtype
wildtype_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_wildtype.txt")
wildtype_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_wildtype.txt")
wildtype_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_wildtype.txt")
wildtype_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_wildtype.txt")
#merge
wildtype = merge(merge(merge(wildtype_epPCR_stabilityPCA[,.(Pos = NA,WT,STOP,STOP_readthrough,Nham_nt,Nham_aa,Nmut_codons,
                                                mean_count_sEP = mean_count,s_fitness_EP = 0,
                                                s_sigma_EP = 1/(sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))],
                 wildtype_epPCR_bindingPCA[,.(Pos = NA,WT,STOP,STOP_readthrough,Nham_nt,Nham_aa,Nmut_codons,
                                                mean_count_bEP = mean_count,b_fitness_EP = 0,
                                                b_sigma_EP = 1/(sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))]),
                 wildtype_NM_stabilityPCA[,.(Pos = NA,WT,STOP,STOP_readthrough,Nham_nt,Nham_aa,Nmut_codons,
                                              mean_count_sNM = mean_count,s_fitness_NM = 0,
                                              s_sigma_NM = 1/(sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))]),
                 wildtype_NM_bindingPCA[,.(Pos = NA,WT,STOP,STOP_readthrough,Nham_nt,Nham_aa,Nmut_codons,
                                             mean_count_bNM = mean_count,b_fitness_NM = 0,
                                             b_sigma_NM = 1/(sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))])
                 
## singles
singles_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_singles.txt")
singles_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_singles.txt")
singles_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_singles.txt")
singles_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_singles.txt")
#merge
singles = merge(merge(merge(singles_NM_stabilityPCA[,.(Pos,WT_AA,Mut,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                                       mean_count_sNM = mean_count,
                                                       s_fitness_NM = fitness,s_sigma_NM = sigma)],
                            singles_epPCR_stabilityPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                                          mean_count_sEP = mean_count,
                                                          s_fitness_EP = fitness,s_sigma_EP = sigma)],all=T),
                      singles_NM_bindingPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                               mean_count_bNM = mean_count,
                                               b_fitness_NM = fitness,b_sigma_NM = sigma)],all=T),
                singles_epPCR_bindingPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                            mean_count_bEP = mean_count,
                                            b_fitness_EP = fitness,b_sigma_EP = sigma)],all=T)

## doubles
doubles_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_doubles.txt")
doubles_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_doubles.txt")
doubles_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_doubles.txt")
doubles_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_doubles.txt")
#merge
doubles = merge(merge(merge(doubles_NM_stabilityPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                                       mean_count_sNM = mean_count,
                                                       s_fitness_NM = fitness_uncorr,s_sigma_NM = sigma_uncorr)],
                            doubles_epPCR_stabilityPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                                          mean_count_sEP = mean_count,
                                                          s_fitness_EP = fitness_uncorr,s_sigma_EP = sigma_uncorr)],all=T),
                      doubles_NM_bindingPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                               mean_count_bNM = mean_count,
                                               b_fitness_NM = fitness_uncorr,b_sigma_NM = sigma_uncorr)],all=T),
                doubles_epPCR_bindingPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
                                            mean_count_bEP = mean_count,
                                            b_fitness_EP = fitness_uncorr,b_sigma_EP = sigma_uncorr)],all=T)

### variant statistics
singles[,.N,.(!is.na(s_fitness_NM),!is.na(s_fitness_EP),!is.na(b_fitness_NM),!is.na(b_fitness_EP))]
cutoff=25
singles[,.N,.(mean_count_sNM>cutoff,mean_count_sEP>cutoff,mean_count_bNM>cutoff,mean_count_bEP>cutoff)]
#512 variants in all four datasets
#1000 variants only in NM datasets
# ~150 variants in other combinations

### plot input count distributions between the datasets
ggpairs(singles,aes(color=factor(Nham_nt),alpha=0.1),
        columns = grep("mean_count",names(singles),value=T),
        lower = list(continuous = ggpairs_log10_points),
        diag = list(continuous = ggpairs_log10_diagonal),
        upper = list(continous = ggpairs_log10_cor))
ggsave("results/preprocessing/PDZ_singles_count_distribution_4samples.pdf",width=8,height=8)
# as expected datasets from epPCR or NM are nearly identical
# in NM all codon changes with different # nt mutations are equally abundant
# in epPCR, there's huge differences
# single nt variants have similar abundance in epPCR and NM; all higher nt changes are more abundant in NM
singles[!is.na(s_fitness_EP),.(.N,mean(mean_count_sEP)),Nham_nt][order(Nham_nt)]
singles[!is.na(s_fitness_NM),.(.N,mean(mean_count_sNM)),Nham_nt][order(Nham_nt)]
#well kind of; actually average NM counts are higher

### fitness versus input counts
X=melt(singles,measure.vars = grep("mean_count",names(singles),value=T))[,.(Pos,Mut,mean_count = value,variable)]
Y=melt(singles,measure.vars = grep("fitness",names(singles),value=T))[,.(Pos,Mut,fitness = value,variable)]
X[,id := strsplit(as.character(variable),"_")[[1]][3],variable]
Y[,id := paste0(strsplit(as.character(variable),"_")[[1]][c(1,3)],collapse=""),variable]
ggplot(merge(X[,.(mean_count,id,Pos,Mut)],Y[,.(fitness,id,Pos,Mut)]),aes(mean_count,fitness)) +
  geom_point() +
  scale_x_log10() +
  geom_vline(xintercept = 15) +
  facet_wrap(id~.)
ggsave("results/preprocessing/PDZ_mean_count_threshold.pdf",width=8,height=8)
## 25 mean_count seems to be a good cutoff that covers most of the fitness range


#merge fitness and sigma from NM and epPCR
ggplot(singles,aes(s_fitness_EP,s_fitness_NM,color=mean_count_sNM > 25 & mean_count_sEP > 25)) +
  geom_density2d() + 
  geom_point() +
  geom_smooth(method="lm") +
  geom_abline()
ggsave("results/preprocessing/PDZ_stabilityPCA_fitness_comparision1.pdf",width=8,height=8)
ggplot(melt(singles,measure.vars = grep("s_fitness",names(singles),value=T)),aes(value,color=variable)) +
  geom_density()
ggsave("results/preprocessing/PDZ_stabilityPCA_fitness_comparision2.pdf",width=5,height=5)
#distributions don't quite match
#use scale and shift procedure to make them more similar

# normalize fitness by scaling and shifting
set.seed(1603)
minF = function(p) {
  F_norm = (F_data + matrix(p[(Nsamples+1):(2*Nsamples)],nrow=nrow(F_data),ncol=Nsamples,byrow = T)) * matrix(p[1:Nsamples],nrow=nrow(F_data),ncol=Nsamples,byrow = T)
  F_avg = rowMeans(F_data + matrix(p[(Nsamples+1):(2*Nsamples)],nrow=nrow(F_data),ncol=Nsamples,byrow = T))
  diffF = sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}
F_data = singles[mean_count_sNM > 100 & mean_count_sEP > 100,as.matrix(.SD),.SDcols = grep("s_fitness_(NM|EP)$",names(singles))]
Nsamples = 2
x=nlm(minF,rep(c(1,0),each=Nsamples))
print(x)
p = x$estimate
# p[1:Nsamples] = p[1:Nsamples]/p[1] # used in DiMSum
fitness_norm_model = data.table(t(p))
names(fitness_norm_model) = c("scale_sNM","scale_sEP","shift_sNM","shift_sEP")

#wild-type correction such that mean(wild-type) = 0
wt_corr = wildtype[,rowMeans((.SD + unlist(fitness_norm_model[,.SD,,.SDcols = grep("shift",names(fitness_norm_model))])) * 
                                      unlist(fitness_norm_model[,.SD,,.SDcols = grep("scale",names(fitness_norm_model))])),
                    ,.SDcols = grep("s_fitness",names(wildtype))]
## normalize fitness values
singles[,s_fitness_NM_norm := (s_fitness_NM + fitness_norm_model$shift_sNM)*fitness_norm_model$scale_sNM - wt_corr]
singles[,s_fitness_EP_norm := (s_fitness_EP + fitness_norm_model$shift_sEP)*fitness_norm_model$scale_sEP - wt_corr] #why is this happening? looks like the scaling should be > 1 to increase spread
ggplot(melt(singles[mean_count_sNM > 100 & mean_count_sEP > 100],measure.vars = grep("s_fitness",names(singles),value=T)),aes(exp(value),color=variable)) +
  geom_density() +
  geom_vline(xintercept = 1)
# this is caused by low count variants


# singles[,s_fitness := sum(c(s_fitness_NM/s_sigma_NM^2,s_fitness_EP/s_sigma_EP^2),na.rm=T)/sum(c(s_sigma_NM^-2,s_sigma_EP^-2),na.rm=T),.(Pos,Mut)]
# ggpairs(singles,columns=grep("s_fitness",names(singles),value=T))

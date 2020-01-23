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

## Colors
col_purple = "#9161A8"
col_blue =  "#0066CC"
col_orange = "#F7941E"
col_red = "#EF4136"

### load the different data files
## wildtype
wildtype_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_wildtype.txt")
wildtype_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_wildtype.txt")
wildtype_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_wildtype.txt")
wildtype_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_wildtype.txt")
wt_AA_seq = wildtype_epPCR_stabilityPCA$aa_seq
wt_AA_seq_split = strsplit(wt_AA_seq,"")[[1]]
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
# singles_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_singles.txt")
# singles_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_singles.txt")
# singles_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_singles.txt")
# singles_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_singles.txt")
#merge
# singles = merge(merge(merge(singles_NM_stabilityPCA[,.(Pos,WT_AA,Mut,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                        mean_count_sNM = mean_count,
#                                                        s_fitness_NM = fitness,s_sigma_NM = sigma)],
#                             singles_epPCR_stabilityPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                           mean_count_sEP = mean_count,
#                                                           s_fitness_EP = fitness,s_sigma_EP = sigma)],all=T),
#                       singles_NM_bindingPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                mean_count_bNM = mean_count,
#                                                b_fitness_NM = fitness,b_sigma_NM = sigma)],all=T),
#                 singles_epPCR_bindingPCA[,.(Pos,WT_AA,Mut,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                             mean_count_bEP = mean_count,
#                                             b_fitness_EP = fitness,b_sigma_EP = sigma)],all=T)

# ## doubles
# doubles_epPCR_stabilityPCA = fread("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/fitness_doubles.txt")
# doubles_epPCR_bindingPCA = fread("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/fitness_doubles.txt")
# doubles_NM_stabilityPCA = fread("dataset/PDZ/01c-PDZ_NM_stabilityPCA/fitness_doubles.txt")
# doubles_NM_bindingPCA = fread("dataset/PDZ/01d-PDZ_NM_bindingPCA/fitness_doubles.txt")
# #merge
# doubles = merge(merge(merge(doubles_NM_stabilityPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                        mean_count_sNM = mean_count,
#                                                        s_fitness_NM = fitness_uncorr,s_sigma_NM = sigma_uncorr)],
#                             doubles_epPCR_stabilityPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                           mean_count_sEP = mean_count,
#                                                           s_fitness_EP = fitness_uncorr,s_sigma_EP = sigma_uncorr)],all=T),
#                       doubles_NM_bindingPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                                mean_count_bNM = mean_count,
#                                                b_fitness_NM = fitness_uncorr,b_sigma_NM = sigma_uncorr)],all=T),
#                 doubles_epPCR_bindingPCA[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nham_aa,Nmut_codons,STOP,STOP_readthrough,
#                                             mean_count_bEP = mean_count,
#                                             b_fitness_EP = fitness_uncorr,b_sigma_EP = sigma_uncorr)],all=T)


##############################################
### intermediate data files before merging ###
##############################################
#### NM stabiltiy
load("dataset/PDZ/01c-PDZ_NM_stabilityPCA/01c-PDZ_NM_stabilityPCA_fitness_intermediate.RData")
#normalize by generations
all_variants[,g1 := 4.8]
all_variants[,g2 := 5]
all_variants[,g3 := 5.1]
all_variants[,fitness1 := fitness1_uncorr / g1]
all_variants[,fitness2 := fitness2_uncorr / g2]
all_variants[,fitness3 := fitness3_uncorr / g3]
all_variants[,sigma1 := sigma1_uncorr / g1]
all_variants[,sigma2 := sigma2_uncorr / g2]
all_variants[,sigma3 := sigma3_uncorr / g3]
#merge fitness, error and input counts
all_variants[,fitness := rowSums(.SD[,1:3]/(.SD[,(3+1):(2*3)]^2),na.rm=T) / 
               rowSums(1/(.SD[,(3+1):(2*3)]^2),na.rm=T),,
             .SDcols = c(grep("^fitness[123]$",names(all_variants)),
                         grep("^sigma[123]$",names(all_variants)))]
all_variants[,sigma := sqrt(1/rowSums(1/.SD^2,na.rm=T)),,
             .SDcols = grep("^sigma[123]$",names(all_variants))]
all_variants[,mean_count := rowMeans(.SD),,.SDcols = grep("count_e[123]_s0",names(all_variants))]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1,.(aa_seq,count_nt1 = log10(mean_count),fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2,.(aa_seq,count_nt2 = log10(mean_count),fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3,.(aa_seq,count_nt3 = log10(mean_count),fitness_nt3 = fitness)]

X=merge(merge(nt1,nt2,all=T),nt3,all=T)
ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 2 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)) & (count_nt3 > 1 | is.na(count_nt3))))
ggsave("results/preprocessing/PDZ_Nhamnt_count_fitness_dependencies_stabilityNM.pdf",width=6,height=6)
#if Nham_nt==1, variants below 10^2 read counts are sequencing errors
#if Nham_nt > 1, variants below 10^1 read counts are sequencing errors
#also there's a proportionality between readcounts of the same AA variants that could be taken into account


## merge based on hamming-distance specific read thresholds
all_data_rest = all_variants[Nham_aa == 1 & ((Nham_nt == 1 & mean_count > 100) | (Nham_nt > 1 & mean_count > 10)),
             .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
               sigma = sqrt(1/sum(sigma^-2)),
               mean_count = sum(mean_count),Nsyn = .N,
               Nham_aa = unique(Nham_aa),WT = unique(WT),STOP = unique(STOP),STOP_readthrough=unique(STOP_readthrough)),aa_seq]
all_data = all_variants[Nham_aa == 1 & mean_count > 10,.(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
                           sigma = sqrt(1/sum(sigma^-2)),
                           mean_count = sum(mean_count),Nsyn = .N,
                           Nham_aa = unique(Nham_aa),WT = unique(WT),STOP = unique(STOP),STOP_readthrough=unique(STOP_readthrough)),aa_seq]
X=merge(all_data_rest[,.(aa_seq,fitness_rest = fitness,Nsyn_rest=Nsyn)],
      all_data[,.(aa_seq,fitness,Nsyn)])
ggplot(X[Nsyn > 1],aes(fitness,fitness_rest)) +
  geom_point()
#not a huge difference in total, but for ~50 variants 
X[abs(fitness-fitness_rest) > 0.05]

## make this the singles dataset
singles_NM_stabilityPCA = all_data_rest


#### NM binding
load("dataset/PDZ/01d-PDZ_NM_bindingPCA/01d-PDZ_NM_bindingPCA_fitness_intermediate.RData")
#normalize by generations
all_variants[,g1 := 5.1]
all_variants[,g2 := 5.3]
all_variants[,g3 := 5.2]
all_variants[,fitness1 := fitness1_uncorr / g1]
all_variants[,fitness2 := fitness2_uncorr / g2]
all_variants[,fitness3 := fitness3_uncorr / g3]
all_variants[,sigma1 := sigma1_uncorr / g1]
all_variants[,sigma2 := sigma2_uncorr / g2]
all_variants[,sigma3 := sigma3_uncorr / g3]
#merge fitness, error and input counts
all_variants[,fitness := rowSums(.SD[,1:3]/(.SD[,(3+1):(2*3)]^2),na.rm=T) / 
               rowSums(1/(.SD[,(3+1):(2*3)]^2),na.rm=T),,
             .SDcols = c(grep("^fitness[123]$",names(all_variants)),
                         grep("^sigma[123]$",names(all_variants)))]
all_variants[,sigma := sqrt(1/rowSums(1/.SD^2,na.rm=T)),,
             .SDcols = grep("^sigma[123]$",names(all_variants))]
all_variants[,mean_count := rowMeans(.SD),,.SDcols = grep("count_e[123]_s0",names(all_variants))]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1,.(aa_seq,count_nt1 = log10(mean_count),fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2,.(aa_seq,count_nt2 = log10(mean_count),fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3,.(aa_seq,count_nt3 = log10(mean_count),fitness_nt3 = fitness)]

X=merge(merge(nt1,nt2,all=T),nt3,all=T)
ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 2 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)) & (count_nt3 > 1 | is.na(count_nt3))))
#if Nham_nt==1, variants below 10^2 read counts are sequencing errors
#if Nham_nt > 1, variants below 10^1 read counts are sequencing errors
#also there's a proportionality between readcounts of the same AA variants that could be taken into account

## merge based on hamming-distance specific read thresholds
singles_NM_bindingPCA = all_variants[Nham_aa == 1 & ((Nham_nt == 1 & mean_count > 100) | (Nham_nt > 1 & mean_count > 10)),
                             .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
                               sigma = sqrt(1/sum(sigma^-2)),
                               mean_count = sum(mean_count),Nsyn = .N,
                               Nham_aa = unique(Nham_aa),WT = unique(WT),STOP = unique(STOP),STOP_readthrough=unique(STOP_readthrough)),aa_seq]

#### epPCR stability
load("dataset/PDZ/01a-PDZ_epPCR_stabilityPCA/01a-PDZ_epPCR_stabilityPCA_fitness_intermediate.RData")
#normalize by generations
all_variants[,g1 := 5]
all_variants[,g2 := 4.8]
all_variants[,g3 := 5]
all_variants[,fitness1 := fitness1_uncorr / g1]
all_variants[,fitness2 := fitness2_uncorr / g2]
all_variants[,fitness3 := fitness3_uncorr / g3]
all_variants[,sigma1 := sigma1_uncorr / g1]
all_variants[,sigma2 := sigma2_uncorr / g2]
all_variants[,sigma3 := sigma3_uncorr / g3]
#merge fitness, error and input counts
all_variants[,fitness := rowSums(.SD[,1:3]/(.SD[,(3+1):(2*3)]^2),na.rm=T) / 
               rowSums(1/(.SD[,(3+1):(2*3)]^2),na.rm=T),,
             .SDcols = c(grep("^fitness[123]$",names(all_variants)),
                         grep("^sigma[123]$",names(all_variants)))]
all_variants[,sigma := sqrt(1/rowSums(1/.SD^2,na.rm=T)),,
             .SDcols = grep("^sigma[123]$",names(all_variants))]
all_variants[,mean_count := rowMeans(.SD),,.SDcols = grep("count_e[123]_s0",names(all_variants))]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1,.(count_nt1 = log10(unlist(.SD[1,1])),fitness_nt1 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2,.(count_nt2 = log10(unlist(.SD[1,1])),fitness_nt2 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3,.(count_nt3 = log10(unlist(.SD[1,1])),fitness_nt3 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]

X=merge(merge(nt1,nt2,all=T),nt3,all=T)
ggpairs(X,columns = grep("count",names(X)))
#there is not Nham_nt>1 variants with sufficient read counts
# compare to NM data instead
X = merge(nt1,singles_NM_stabilityPCA[,.(aa_seq,fitness,mean_count = log10(mean_count))])
ggpairs(X,columns =c("count_nt1","mean_count"))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 1.5 | is.na(count_nt1)) & (mean_count > 1 | is.na(mean_count)) ))
# correlation is good not matter what the input read count of the epPCR library
## >> use all epPCR data from Nham_nt == 1
singles_epPCR_stabilityPCA = all_variants[Nham_aa == 1 & Nham_nt == 1,
                                     .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
                                       sigma = sqrt(1/sum(sigma^-2)),
                                       mean_count = sum(mean_count),Nsyn = .N,
                                       Nham_aa = unique(Nham_aa),WT = unique(WT),STOP = unique(STOP),STOP_readthrough=unique(STOP_readthrough)),aa_seq]

#### epPCR binding
load("dataset/PDZ/01b-PDZ_epPCR_bindingPCA/01b-PDZ_epPCR_bindingPCA_fitness_intermediate.RData")
#normalize by generations
all_variants[,g1 := 5]
all_variants[,g2 := 5]
all_variants[,g3 := 4.8]
all_variants[,fitness1 := fitness1_uncorr / g1]
all_variants[,fitness2 := fitness2_uncorr / g2]
all_variants[,fitness3 := fitness3_uncorr / g3]
all_variants[,sigma1 := sigma1_uncorr / g1]
all_variants[,sigma2 := sigma2_uncorr / g2]
all_variants[,sigma3 := sigma3_uncorr / g3]
#merge fitness, error and input counts
all_variants[,fitness := rowSums(.SD[,1:3]/(.SD[,(3+1):(2*3)]^2),na.rm=T) / 
               rowSums(1/(.SD[,(3+1):(2*3)]^2),na.rm=T),,
             .SDcols = c(grep("^fitness[123]$",names(all_variants)),
                         grep("^sigma[123]$",names(all_variants)))]
all_variants[,sigma := sqrt(1/rowSums(1/.SD^2,na.rm=T)),,
             .SDcols = grep("^sigma[123]$",names(all_variants))]
all_variants[,mean_count := rowMeans(.SD),,.SDcols = grep("count_e[123]_s0",names(all_variants))]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1,.(count_nt1 = log10(unlist(.SD[1,1])),fitness_nt1 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2,.(count_nt2 = log10(unlist(.SD[1,1])),fitness_nt2 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3,.(count_nt3 = log10(unlist(.SD[1,1])),fitness_nt3 = unlist(.SD[1,2])),aa_seq,.SDcols = c("mean_count","fitness")]

X=merge(merge(nt1,nt2,all=T),nt3,all=T)
ggpairs(X,columns = grep("count",names(X)))
#there is not Nham_nt>1 variants with sufficient read counts
# compare to NM data instead
X = merge(nt1,singles_NM_bindingPCA[,.(aa_seq,fitness,mean_count = log10(mean_count))])
ggpairs(X,columns =c("count_nt1","mean_count"))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 1.5 | is.na(count_nt1)) & (mean_count > 2.5 | is.na(mean_count)) ))
# correlation is good not matter what the input read count of the epPCR library
## >> use all epPCR data from Nham_nt == 1
singles_epPCR_bindingPCA = all_variants[Nham_aa == 1 & Nham_nt == 1,
                                          .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
                                            sigma = sqrt(1/sum(sigma^-2)),
                                            mean_count = sum(mean_count),Nsyn = .N,
                                            Nham_aa = unique(Nham_aa),WT = unique(WT),STOP = unique(STOP),STOP_readthrough=unique(STOP_readthrough)),aa_seq]

#####################################
####### merge all single data #######
#####################################
singles = merge(merge(merge(singles_NM_stabilityPCA[,.(aa_seq,STOP,STOP_readthrough,
                                                       mean_count_sNM = mean_count,
                                                       s_fitness_NM = fitness,s_sigma_NM = sigma)],
                            singles_epPCR_stabilityPCA[,.(aa_seq,STOP,STOP_readthrough,
                                                          mean_count_sEP = mean_count,
                                                          s_fitness_EP = fitness,s_sigma_EP = sigma)],all=T),
                      singles_NM_bindingPCA[,.(aa_seq,STOP,STOP_readthrough,
                                               mean_count_bNM = mean_count,
                                               b_fitness_NM = fitness,b_sigma_NM = sigma)],all=T),
                singles_epPCR_bindingPCA[,.(aa_seq,STOP,STOP_readthrough,
                                            mean_count_bEP = mean_count,
                                            b_fitness_EP = fitness,b_sigma_EP = sigma)],all=T)

## redo Pos/WT_AA/Mut stuff
singles[,Pos := which(strsplit(aa_seq,"")[[1]] != wt_AA_seq_split),aa_seq]
singles[,WT_AA := wt_AA_seq_split[Pos],Pos]
singles[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]


### variant statistics
singles[,.N,.(!is.na(s_fitness_NM),!is.na(s_fitness_EP),!is.na(b_fitness_NM),!is.na(b_fitness_EP))]
cutoff=25
singles[,.N,.(mean_count_sNM>cutoff,mean_count_sEP>cutoff,mean_count_bNM>cutoff,mean_count_bEP>cutoff)]
#389 variants in all four datasets
#1057 variants only in NM datasets
# ~150 variants in other combinations

### plot input count distributions between the datasets
ggpairs(singles,aes(alpha=0.1),
        columns = grep("mean_count",names(singles),value=T),
        lower = list(continuous = ggpairs_log10_points),
        diag = list(continuous = ggpairs_log10_diagonal),
        upper = list(continous = ggpairs_log10_cor))
ggsave("results/preprocessing/PDZ_singles_count_distribution_4samples.pdf",width=8,height=8)
# as expected datasets from epPCR or NM are nearly identical
# single nt variants have similar abundance in epPCR and NM but different biases
singles[!is.na(s_fitness_EP),.(.N,mean(mean_count_sEP))]
singles[!is.na(s_fitness_NM),.(.N,mean(mean_count_sNM))]


### fitness versus input counts
X=melt(singles,measure.vars = grep("mean_count",names(singles),value=T))[,.(Pos,Mut,mean_count = value,variable)]
Y=melt(singles,measure.vars = grep("fitness",names(singles),value=T))[,.(Pos,Mut,fitness = value,variable)]
X[,id := strsplit(as.character(variable),"_")[[1]][3],variable]
Y[,id := paste0(strsplit(as.character(variable),"_")[[1]][c(1,3)],collapse=""),variable]
ggplot(merge(X[,.(mean_count,id,Pos,Mut)],Y[,.(fitness,id,Pos,Mut)]),aes(mean_count,fitness)) +
  geom_point() +
  scale_x_log10() +
  geom_vline(xintercept = c(15,50)) +
  facet_wrap(id~.)
ggsave("results/preprocessing/PDZ_mean_count_threshold.pdf",width=8,height=8)
## 15 mean_count seems to be a good cutoff for epPCR
## for NM there has to be a higher cutoff, due to the overabundance of wildtype supplying many seq-error variants in output samples and thus skewing fitness values in lowly abundant variants

#merge fitness and sigma from NM and epPCR
ggplot(singles,aes(s_fitness_EP,s_fitness_NM,color=mean_count_sNM > 50 )) +
  geom_density2d() + 
  geom_point() +
  geom_smooth(method="lm") +
  geom_abline()
ggsave("results/preprocessing/PDZ_stabilityPCA_fitness_comparision1.pdf",width=8,height=8)

## >> remove low fitness NM data with lower than 50 counts
# singles[mean_count_sNM < 50 & s_fitness_NM < -.25,s_sigma_NM := NA]
# singles[mean_count_bNM < 50 & b_fitness_NM < -.25,b_sigma_NM := NA]

ggplot(melt(singles[mean_count_sNM > 50],measure.vars = grep("s_fitness",names(singles),value=T)),aes(value,color=variable)) +
  geom_density()
ggsave("results/preprocessing/PDZ_stabilityPCA_fitness_comparision2.pdf",width=5,height=5)
#distributions don't quite match; but that's also because low fitness variants were selectively excluded in NM data

##### normalize fitness by scaling and shifting
set.seed(1603)
minF = function(p) {
  F_norm = (F_data + matrix(p[(Nsamples+1):(2*Nsamples)],nrow=nrow(F_data),ncol=Nsamples,byrow = T)) * matrix(p[1:Nsamples],nrow=nrow(F_data),ncol=Nsamples,byrow = T)
  F_avg = rowMeans(F_data + matrix(p[(Nsamples+1):(2*Nsamples)],nrow=nrow(F_data),ncol=Nsamples,byrow = T))
  diffF = sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}

### stability data
F_data = singles[mean_count_sNM > 100 & mean_count_sEP > 100,as.matrix(.SD),.SDcols = grep("s_fitness_(NM|EP)$",names(singles))]
Nsamples = 2
x=nlm(minF,rep(c(1,0),each=Nsamples))
print(x)
p = x$estimate
fitness_norm_model = data.table(t(p))
names(fitness_norm_model) = c("scale_sNM","scale_sEP","shift_sNM","shift_sEP")

#wild-type correction such that mean(wild-type) = 0
wt_corr = wildtype[,rowMeans((.SD + unlist(fitness_norm_model[,.SD,,.SDcols = grep("shift_s",names(fitness_norm_model))])) *
                               unlist(fitness_norm_model[,.SD,,.SDcols = grep("scale_s",names(fitness_norm_model))])),
                   ,.SDcols = grep("s_fitness_(NM|EP)$",names(wildtype))]
## normalize fitness values
singles[,s_fitness_NM_norm := (s_fitness_NM + fitness_norm_model$shift_sNM)*fitness_norm_model$scale_sNM - wt_corr]
singles[,s_fitness_EP_norm := (s_fitness_EP + fitness_norm_model$shift_sEP)*fitness_norm_model$scale_sEP - wt_corr] #why is this happening? looks like the scaling should be > 1 to increase spread
singles[,s_sigma_NM_norm := s_sigma_NM*fitness_norm_model$scale_sNM]
singles[,s_sigma_EP_norm := s_sigma_EP*fitness_norm_model$scale_sEP]
ggplot(melt(singles[mean_count_sNM > 100 & mean_count_sEP > 100],measure.vars = grep("s_fitness",names(singles),value=T)),
       aes((value),color=variable)) +
  geom_density() +
  geom_vline(xintercept = 0)
# this is caused by low count variants, so it's fine

# merge fitness and error
singles[,s_fitness := sum(c(s_fitness_NM_norm/s_sigma_NM_norm^2,s_fitness_EP_norm/s_sigma_EP_norm^2),na.rm=T)/sum(c(s_sigma_NM_norm^-2,s_sigma_EP_norm^-2),na.rm=T),.(Pos,Mut)]
singles[,s_sigma := sqrt(1/sum(c(s_sigma_NM_norm^-2,s_sigma_EP_norm^-2),na.rm=T)),.(Pos,Mut)]

#plot comparision
ggpairs(singles,columns=grep("s_fitness(|_NM|_EP)$",names(singles),value=T))
ggsave(filename = "results/preprocessing/PDZ_stability_fitness_merge.pdf",width=6,height=6)


### binding data
F_data = singles[mean_count_bNM > 100 & mean_count_bEP > 100,as.matrix(.SD),.SDcols = grep("b_fitness_(NM|EP)$",names(singles))]
Nsamples = 2
x=nlm(minF,rep(c(1,0),each=Nsamples))
print(x)
p = x$estimate
fitness_norm_model = cbind(fitness_norm_model,data.table(t(p)))
names(fitness_norm_model)[5:8] = c("scale_bNM","scale_bEP","shift_bNM","shift_bEP")

#wild-type correction such that mean(wild-type) = 0
wt_corr = wildtype[,rowMeans((.SD + unlist(fitness_norm_model[,.SD,,.SDcols = grep("shift_b",names(fitness_norm_model))])) *
                               unlist(fitness_norm_model[,.SD,,.SDcols = grep("scale_b",names(fitness_norm_model))])),
                   ,.SDcols = grep("b_fitness_(NM|EP)$",names(wildtype))]
## normalize fitness values
singles[,b_fitness_NM_norm := (b_fitness_NM + fitness_norm_model$shift_bNM)*fitness_norm_model$scale_bNM - wt_corr]
singles[,b_fitness_EP_norm := (b_fitness_EP + fitness_norm_model$shift_bEP)*fitness_norm_model$scale_bEP - wt_corr] #why is this happening? looks like the scaling should be > 1 to increase spread
singles[,b_sigma_NM_norm := b_sigma_NM*fitness_norm_model$scale_bNM]
singles[,b_sigma_EP_norm := b_sigma_EP*fitness_norm_model$scale_bEP]
ggplot(melt(singles[mean_count_bNM > 100 & mean_count_bEP > 100],measure.vars = grep("b_fitness",names(singles),value=T)),aes((value),color=variable)) +
  geom_density() +
  geom_vline(xintercept = 0)
# this is caused by low count variants, so it's fine

#### merge fitness and error
singles[,b_fitness := sum(c(b_fitness_NM_norm/b_sigma_NM_norm^2,b_fitness_EP_norm/b_sigma_EP_norm^2),na.rm=T)/sum(c(b_sigma_NM_norm^-2,b_sigma_EP_norm^-2),na.rm=T),.(Pos,Mut)]
singles[,b_sigma := sqrt(1/sum(c(b_sigma_NM_norm^-2,b_sigma_EP_norm^-2),na.rm=T)),.(Pos,Mut)]
# singles[,b_fitness := sum(c(b_fitness_NM/b_sigma_NM^2,b_fitness_EP/b_sigma_EP^2),na.rm=T)/sum(c(b_sigma_NM^-2,b_sigma_EP^-2),na.rm=T),.(Pos,Mut)]
# singles[,b_sigma := sqrt(1/sum(c(b_sigma_NM^-2,b_sigma_EP^-2),na.rm=T)),.(Pos,Mut)]

#plot comparision; also to McLaughlin data
mclaughlin = fread("dataset/PDZ/McLaughlin2012_DLG4_singlemutants.txt")
mclaughlin[,Pos := as.integer(gsub("[A-Z]","",mutant))-310,mutant]
mclaughlin[,WT_AA := gsub("[0-9]+[A-Z]","",mutant),mutant]
mclaughlin[,Mut := gsub("[A-Z][0-9]+","",mutant),mutant]
singles_RN = merge(singles,mclaughlin[,.(Pos,Mut,b_fitness_mclaughlin=CRIPT)],by=c("Pos","Mut"))
ggpairs(singles_RN,columns=grep("b_fitness(|_NM|_EP|_mclaughlin)$",names(singles_RN),value=T))
ggsave(filename = "results/preprocessing/PDZ_binding_fitness_merge.pdf",width=8,height=8)

p1=ggplot(singles_RN,aes((b_fitness_mclaughlin),(s_fitness),color=(b_fitness)-(s_fitness))) +
  geom_point() +
  scale_color_gradient2(midpoint = 0,mid="grey") +
  labs(x="fitness McLaughlin2012",y="fitness stabilityPCA",color="deltaF ddPCA")
p2=ggplot(singles_RN,aes((b_fitness_mclaughlin),(b_fitness),color=(b_fitness)-(s_fitness))) +
  geom_point() +
  scale_color_gradient2(midpoint = 0,mid="grey") +
  labs(x="fitness McLaughlin2012",y="fitness bindingPCA",color="deltaF ddPCA")
P=grid.arrange(p1,p2,nrow=1)
ggsave(plot=P,filename = "results/PDZ_ddPCA_McLaughlin2012_comparision.pdf",width=9,height=4)

##################################
### check library inequalities ###
##################################

X=melt(singles,id.vars=c("Pos","Mut"),measure.vars = grep("mean_count",names(singles)))
ggplot(X,aes(x=Pos,group=Pos,value)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(variable~.)
ggsave("results/preprocessing/PDZ_library_inequalities_byPos.pdf",width=10,height=6)

Tm_NM = fread("dataset/PDZ/PDZ_Sigma_NNKlibrary.txt")
Tm_NM[,Pos := as.integer(gsub("PDZ_","",Name)),Name]
singles_Tm_NM = merge(singles[,.(mean_count_sNM=mean(mean_count_bNM,na.rm=T)),Pos],Tm_NM[,.(Pos,Tm)])
singles_Tm_NM = merge(singles[!is.na(mean_count_sNM),.(Pos,mean_count_sNM)],Tm_NM[,.(Pos,Tm)])
singles_Tm_NM[,cor(mean_count_sNM,Tm)]
ggplot(singles_Tm_NM,aes(Tm,mean_count_sNM)) +
  geom_point() +
  scale_y_log10() +
  geom_smooth()
ggsave("results/preprocessing/PDZ_library_inequalities_NMvsTm.pdf",width=6,height=6)
# weak dependence on Tm


###############################
### integrate PDB structure ###
###############################

### contact map
pairdistances_from_PDB(input_file = "dataset/PDZ/PDB/pdb1be9.ent",dataset_dir="",
                       given_chainids = c("A","B"),
                       aa_seq = list(paste0(unique(singles[,.(Pos,WT_AA)])$WT_AA,collapse=""),"KQTSV"),
                       idx_pdb_start = c(311,5),
                       idx_DMS_start = c(1,1), debug_this = F)
contactmap_AB = fread("processed_data/PDB_contactmap_pdb1be9_AB.txt")
singles = merge(singles,contactmap_AB[,.(HAmin_ligand = min(HAmin)),Pos1],by.x="Pos",by.y="Pos1",all.x=T)

### load RSA values
RSA_chainA = fread("dataset/PDZ/pdb1be9_A.rsa",skip=9)
names(RSA_chainA) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
RSA_chainAB = fread("dataset/PDZ/pdb1be9_AB.rsa",skip=9)
names(RSA_chainAB) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
# add to singles data.table
singles = merge(singles,RSA_chainA[chain == "A",.(Pos=Pos-310,RSA_unbound = RSA_all_rel)],by="Pos",all.x=T)
singles = merge(singles,RSA_chainAB[chain == "A",.(Pos=Pos-310,RSA_bound = RSA_all_rel)],by="Pos",all.x=T)
# >>> all residues that loose RSA when bound to ligand are on interface with ligand
ggplot(singles,aes(RSA_unbound,RSA_bound,color=HAmin_ligand < 5)) + geom_point(size=3)

## split data into core, surface and binding
ggplot(singles,aes(RSA_unbound,s_fitness)) + geom_point() +geom_smooth()
# >>> set core threshodl as RSA <= 10
threshold_RSA = 10
threshold_ligand_dist = 5 
singles[RSA_unbound <= threshold_RSA,type := "core"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "ligand_binding"]
singles[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist,type := "surface"]
singles[,type := factor(type,levels=c("core","surface","ligand_binding"))]
singles[,.(Npos = length(unique(Pos))),type]



### >>save data files
write.table(singles,file = "processed_data/PDZ_singles_alldata.txt",quote=F,row.names=F)
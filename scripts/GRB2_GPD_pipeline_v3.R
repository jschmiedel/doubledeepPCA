#####################
##### GRB2 data #####
#####################

#GPD / stability data

#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v2/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v2/")

#load required packages
require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)

doplot = F

#############################################################################################
##### preprocess data (calculate fitness scores and errors, set quality thresholds etc) #####
#############################################################################################

###### load DiMSum datasets
#################################
######## stability data #########
#################################
load("dataset/DiMSum_GRB2_GPD_variant_data_merge_Q25.RData")

if (doplot == T) {
  #check subs and indels in library
  X=variant_data_merge[,.(.N,reads = sum(input1_e1_s0_bNA_count)),.(Nmut_nt,Nsub_nt)]
  X[,freq := reads/sum(reads)]
  X[order(-freq)][1:10]
  X[,.(sum_freq = sum(freq)),.(sub_diff=Nmut_nt - Nsub_nt)][order(-sum_freq)][1:10]
  ggplot(X[order(-freq)][1:20],aes(Nmut_nt,freq,fill=factor(Nmut_nt-Nsub_nt))) +
    geom_bar(stat='identity') +
    labs(fill="indels",x="# nt mutations",y="seq read frequency")
  ggsave("dataset/preprocessing/GRB2_GPD_sub_indel_frequency.pdf")
}

#discard indel data and all higher order AA mutated variants
dataset = variant_data_merge[Nins_aa == 0 & Ndel_aa == 0 & Nsub_aa <= 2]

if (doplot == T) {
  ggplot(dataset[between(Nmut_nt,1,4)],aes(input1_e1_s0_bNA_count,..scaled..,color=factor(Nmut_nt))) +  geom_density() +
    scale_x_log10() +
    facet_grid(Nmut_aa ~ Nmut_nt <= 2)
  ggsave("dataset/preprocessing/GRB2_GPD_seqread_distributions.pdf",width=7,height=7)
}


dataset[,.N,STOP]
#these data include the STOP codon at the last AA position!!!! (therefore there's Nmut_aa == 2 & Nmut_nt == 1, which are single NT inserts that shift the last two AA)
# > sort out variants which have this STOP codon mutated
dataset[,last_STOP := strsplit(aa_seq,"")[[1]][57] == "*",aa_seq]
dataset = dataset[last_STOP==T]
#later recode STOP variable

#get rid of variants with Nmut_nt+2 > Nmut_aa
dataset = dataset[Nmut_nt <= Nmut_aa+2]
dataset[,.(variants=.N,total_input=sum(input1_e1_s0_bNA_count)),.(Nmut_aa,Nmut_nt)][order(Nmut_aa,Nmut_nt)]

if (doplot == T) {
  ggplot(dataset,aes((output1_e1_s1_b1_count+0.001)/(input1_e1_s0_bNA_count+0.001),..scaled..,color=factor(Nmut_nt))) +
    geom_density() + 
    scale_x_log10() +
    geom_vline(xintercept = dataset[Nmut_nt==0,(output1_e1_s1_b1_count+0.001)/(input1_e1_s0_bNA_count+0.001)]) +
    facet_grid(Nmut_aa ~ .) +
    labs(x="fitness")
  ggsave("dataset/preprocessing/GRB2_GPD_fitness_higherorder_ntmut.pdf",width=7,height=7)
  #peak at 1 comes from variants observed once in input & output
  ggplot(dataset[Nmut_aa==0 & Nmut_nt == 2],aes(input1_e1_s0_bNA_count,output1_e1_s1_b1_count)) +
    geom_hex() +
    geom_abline(color="red") +
    scale_x_log10() +
    scale_y_log10()
  #but mean fitness is same across all Nmut_nt combinations
  dataset[,mean(output1_e1_s1_b1_count)/mean(input1_e1_s0_bNA_count),.(Nmut_aa,Nmut_nt)][order(Nmut_aa,Nmut_nt)]
}

#### add up counts from synonymous variants
idx = names(dataset)[grep(names(dataset),pattern="_count$")]
for (i in seq_along(idx)) {
  dataset[Nmut_aa == 0,paste0(idx[i],"_agg") := .SD,,.SDcols = idx[i]]
  dataset[Nmut_aa > 0,paste0(idx[i],"_agg") := sum(.SD),aa_seq,.SDcols = idx[i]]
}
dataset[Nmut_aa > 0,Nmut_nt := NA]
dataset_syn = unique(dataset[,.SD,aa_seq,.SDcols = c("Nmut_aa","Nmut_nt","WT","STOP",names(dataset)[grep(names(dataset),pattern="_agg$")])])

if (doplot == T) {
  #count distribution
  ggplot(melt(dataset_syn[Nmut_aa > 0],id.vars = c("aa_seq","Nmut_aa","Nmut_nt","WT","STOP")),aes(value,color=factor(Nmut_aa))) +
    geom_density() +
    scale_x_log10() +
    facet_wrap(variable ~ .)
  ggsave("dataset/preprocessing/GRB2_GPD_seqread_syn_distributions.pdf",width=7,height=7)
}

# calculate fitness and count-based error
doublings = 5.11
for (E in 1:3) {
  c_total_in = variant_data_merge[,sum(.SD[,1]),,
                                  .SDcols = paste0("input",E,"_e",E,"_s0_bNA_count")]
  c_total_out = variant_data_merge[,sum(.SD[,1]),,
                                   .SDcols = paste0("output",E,"_e",E,"_s1_b1_count")]
  f_wt_corr = dataset_syn[WT == T,
                          doublings + log2((.SD[,1]/c_total_out)/(.SD[,2]/c_total_in)),,
                          .SDcols = c(paste0("output",E,"_e",E,"_s1_b1_count_agg"),paste0("input",E,"_e",E,"_s0_bNA_count_agg"))]
  
  ### non-logged fitness scores (similar to PPI from Diss2019)
  dataset_syn[,paste0("fitness",E,"_uncorr") := (doublings + log2(unlist(.SD[,1])/c_total_out / (unlist(.SD[,2])/c_total_in))) / (rep(t(f_wt_corr),nrow(dataset_syn))),
              ,.SDcols = c(paste0("output",E,"_e",E,"_s1_b1_count_agg"),paste0("input",E,"_e",E,"_s0_bNA_count_agg"))]
  
  ### divide by wild-type correction
  dataset_syn[,paste0("sigma",E,"_uncorr") := sqrt(log2(exp(1))) * sqrt(1/unlist(.SD[,1]) + 1/unlist(.SD[,2])) / rep(t(f_wt_corr),nrow(dataset_syn)),
              ,.SDcols = c(paste0("output",E,"_e",E,"_s1_b1_count_agg"),paste0("input",E,"_e",E,"_s0_bNA_count_agg"),paste0("fitness",E,"_uncorr"))]
}

if (doplot == T) {
  DT = rbind(dataset_syn[,.(count_in = input1_e1_s0_bNA_count_agg,sigma = sigma1_uncorr,fitness = fitness1_uncorr,rep = "rep1",Nmut_aa,Nmut_nt)],
             dataset_syn[,.(count_in = input2_e2_s0_bNA_count_agg,sigma = sigma2_uncorr,fitness = fitness2_uncorr,rep = "rep2",Nmut_aa,Nmut_nt)],
             dataset_syn[,.(count_in = input3_e3_s0_bNA_count_agg,sigma = sigma3_uncorr,fitness = fitness3_uncorr,rep = "rep3",Nmut_aa,Nmut_nt)])
  
  ggplot(DT[Nmut_aa != 0 & count_in > 1000],aes(fitness,color=factor(Nmut_aa),linetype = rep)) +
    geom_density()
  
  ggplot(DT[Nmut_aa > 0],aes(count_in,(fitness))) +
    geom_hex() +
    scale_x_log10() +
    facet_grid(rep ~ Nmut_aa)
  ggsave("dataset/preprocessing/GRB2_GPD_read_fitness_2dhex.pdf",width=7,height=7)
  #severe limitations from low input counts on fitness range
  ggplot(DT,aes(fitness,sigma)) +
    geom_hex() +
    scale_y_log10() +
    facet_wrap( ~ rep)
  ggsave("dataset/preprocessing/GRB2_GPD_fitness_sigma_2dhex.pdf",width=7,height=7)
}

#################################################################
## extract position and mutations from protein_variant factor ####
##################################################################
wildtype = dataset_syn[Nmut_aa==0 & Nmut_nt < 2,]
wt_AAseq_split = strsplit(wildtype[WT==T,aa_seq],"")[[1]][1:56]
wt_AAseq = paste0(wt_AAseq_split,collapse = "")

singles = dataset_syn[Nmut_aa==1,]
singles[,Pos := which(strsplit(aa_seq,"")[[1]][1:56] !=wt_AAseq_split),aa_seq]
singles[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
singles[,WT_AA := wt_AAseq_split[Pos],aa_seq]
singles[,mean_count := rowMeans(.SD),,.SDcols = c("input1_e1_s0_bNA_count_agg","input2_e2_s0_bNA_count_agg","input3_e3_s0_bNA_count_agg")]

if (doplot == T) {
  ggplot(singles,aes(mean_count)) + 
    geom_density() + 
    scale_x_log10() +
    geom_vline(xintercept = 500)
  #lower peak is from singles with AA changes 2 nt away
  
  ggplot(singles,aes(input1_e1_s0_bNA_count_agg,fitness1_uncorr)) +
    geom_point() +
    scale_x_log10() +
    geom_vline(xintercept = 500)
  #only use high-count singles to avoid 'expected fitness' misestimates
  # singles = singles[mean_count > 500,.(Pos,WT_AA,Mut,STOP,mean_count,fitness1 = fitness1_uncorr,sigma1 = sigma1_uncorr,
  #                                      fitness2 = fitness2_uncorr,sigma2 = sigma2_uncorr,
  #                                      fitness3 = fitness3_uncorr,sigma3 = sigma3_uncorr,
  #                                      input1_e1_s0_bNA_count_agg,input2_e2_s0_bNA_count_agg,input3_e3_s0_bNA_count_agg,
  #                                      output1_e1_s1_b1_count_agg,output2_e2_s1_b1_count_agg,output3_e3_s1_b1_count_agg)]
  #change this threshold to 20, as these have full fitness range
}
singles = singles[mean_count > 20,.(Pos,WT_AA,Mut,STOP,mean_count,fitness1 = fitness1_uncorr,sigma1 = sigma1_uncorr,
                                    fitness2 = fitness2_uncorr,sigma2 = sigma2_uncorr,
                                    fitness3 = fitness3_uncorr,sigma3 = sigma3_uncorr,
                                    input1_e1_s0_bNA_count_agg,input2_e2_s0_bNA_count_agg,input3_e3_s0_bNA_count_agg,
                                    output1_e1_s1_b1_count_agg,output2_e2_s1_b1_count_agg,output3_e3_s1_b1_count_agg)]


doubles = dataset_syn[Nmut_aa==2]
doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]][1:56] !=wt_AAseq_split)[1],aa_seq]
doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]][1:56] !=wt_AAseq_split)[2],aa_seq]
doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]

doubles = merge(doubles,singles[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
doubles = merge(doubles,singles[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
doubles[,mean_count := rowMeans(.SD),,.SDcols = c("input1_e1_s0_bNA_count_agg","input2_e2_s0_bNA_count_agg","input3_e3_s0_bNA_count_agg")]

doubles[,.N,s1_mean_count > 500 & s2_mean_count > 500]
doubles[,.N,mean_count > 20]
#only few doubles that have individual amino acid changes 2nt away from wild-type
# for epistasis analysis only AAchanges 1nt away from wildtype will be useful

doubles = doubles[mean_count > 20]

############################
##### merge replicates #####
############################

## assume 5% replicate error [[improve]]
rep_err = 0.05

### merge wildtype and singles
wildtype[,fitness := 0]
wildtype[,sigma := sqrt(1/rowSums(1/(.SD^2 + rep_err^2))),,.SDcols = grep("sigma",names(wildtype))]

singles[,fitness := rowSums(.SD[,1:3]/(.SD[,4:6]^2  + rep_err^2 ),na.rm=T)/
          rowSums(1/(.SD[,4:6]^2 + rep_err^2),na.rm=T),,
        .SDcols = c(grep("fitness",names(singles)),grep("sigma",names(singles)))]
singles[,sigma := sqrt(1/rowSums(1/(.SD[,1:3]^2 + rep_err^2),na.rm=T)),,
        .SDcols = grep("sigma",names(singles))]

### merge doubles
f_vec = doubles[,grep("fitness[123]_uncorr",names(doubles)),with=F]
s_vec = doubles[,grep("sigma[123]_uncorr",names(doubles)),with=F]

doubles[,fitness := rowSums(f_vec/(s_vec^2 + rep_err^2),na.rm=T) / rowSums(1/(s_vec^2 + rep_err^2),na.rm=T)]
doubles[,sigma := sqrt(1 / rowSums(1/(s_vec^2 + rep_err^2),na.rm=T))]

if (doplot == T) {
  ggpairs(rbind(doubles[sample(.N,1000),cbind(.SD[,1],log10(.SD[,2]),Nmut_aa=2),,
                        .SDcols = c("fitness","sigma")],
                singles[,cbind(.SD[,1],log10(.SD[,2]),Nmut_aa=1),,
                        .SDcols = c("fitness","sigma")]),
          aes(alpha=0.1,color=factor(Nmut_aa)),columns = c("fitness","sigma"))
  ggsave("dataset/preprocessing/GRB2_GPD_fitness_merged_uncorrVScond.pdf",width=7.5,height=7.5)
}

#define STOPs
singles[,STOP := Mut == "*"]
doubles[,STOP := (Mut1 == "*" | Mut2 == "*")]

#merge single fitness to doubles DT
doubles = doubles[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,STOP,mean_count,fitness,sigma,
                     fitness1_uncorr,sigma1_uncorr,
                     fitness2_uncorr,sigma2_uncorr,
                     fitness3_uncorr,sigma3_uncorr,
                     s1_mean_count,s2_mean_count)]
# doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
# doubles[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
# doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
# doubles[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]

###### write data tables to processed_data folder
write.table(x = wildtype, file = "dataset/GRB2_GPD_wildtype.txt",
            quote = F,row.names = F, col.names = T)
write.table(x = singles, file = "dataset/GRB2_GPD_singles_readT20.txt",
            quote = F,row.names = F, col.names = T)
write.table(x = doubles, file = "dataset/GRB2_GPD_doubles_readT20.txt",
            quote = F,row.names = F, col.names = T)

### import and process GB1 data from Olson2014

setwd("doubledeepPCA/")

## source functions
filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)

# create directory structure if not present
dir.create("processed_data/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/preprocessing/", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
require(foreach)
require(doMC)
theme_set(theme_bw(base_size = 9))





#############################
### from replicate counts ###
#############################

# calculate fitness and count-based error
#load dataset
work_data <- fread("dataset/GB1/GB1.txt")

#how many replicates?
reps = paste0(gsub("input", "", grep("input", names(work_data), value = T)), collapse = "")
reps_num = as.numeric(strsplit(reps, "")[[1]])

#####################
###### fitness ######
#####################
### add pseudocount
for (j in reps_num) {
	work_data[, paste0("input", j) := .SD + 0.5,
		.SDcols = paste0("input", j)]
	work_data[, paste0("output", j) := .SD + 0.5,
		.SDcols = paste0("output", j)]
}



## calculate fitness
for (j in reps_num) {
  wt_corr = as.numeric(work_data[WT == T, log(.SD[, 2] / .SD[, 1]),,
                         .SDcols = c(grep(paste0("^input", j, "$"),
                         	names(work_data)), grep(paste0("^output", j, "$"), names(work_data)))])
  
  work_data[,paste0("fitness", j) := log(.SD[, 2] / .SD[, 1]) - wt_corr,,
            .SDcols = c(grep(paste0("^input", j, "$"),
            	names(work_data)), grep(paste0("^output", j ,"$"), names(work_data)))]
}

# flag variants that don't have reads in all input/output replicates
work_data[,all_reads := rowSums(.SD > 0) == (2 * nchar(reps)),, 
	.SDcols = grep(paste0("put[", reps,"]$"), names(work_data))]

# find input read threshold for full fitness range
input_count_threshold = work_data[all_reads == T,
	exp(-quantile(.SD, probs = 0.01, na.rm = T)),, .SDcols = grep("fitness", names(work_data))]

#define variants above threshold for later use
work_data[, input_above_threshold := rowSums(.SD > input_count_threshold) == nchar(reps),,
	.SDcols = grep(paste0("input[", reps, "]$"), names(work_data))]

# normalize fitness by scaling and shifting
# set.seed(1603)
# minF = function(p) {
#   F_norm = (F_data + matrix(p[(nchar(reps) + 1) : (2 * nchar(reps))], 
#   		nrow = nrow(F_data), ncol = nchar(reps), byrow = T)) * matrix(p[1 : nchar(reps)], 
#   			nrow = nrow(F_data), ncol = nchar(reps), byrow = T)
#   F_avg = rowMeans(F_data + matrix(p[(nchar(reps) + 1) : (2 * nchar(reps))], 
#   		nrow = nrow(F_data), ncol = nchar(reps), byrow = T))
#   diffF = sqrt(rowSums((F_norm - F_avg)^2))
#   return(sum(diffF))
# }
# F_data = work_data[input_above_threshold == T & all_reads ==T,as.matrix(.SD),.SDcols = grep(paste0("fitness[",reps,"]$"),names(work_data))]
# W_data = work_data[input_above_threshold == T & all_reads ==T,as.matrix(1/.SD^2),.SDcols = grep(paste0("cbe[",reps,"]$"),names(work_data))]
# x = nlm(minF, rep(c(1, 0), each = nchar(reps)))
# # print(x)
# p = x$estimate
# print(p)
# p[1 : nchar(reps)] = p[1 : nchar(reps)] / p[1]

p = c(1, 1, 1, 0, 0, 0)
fitness_norm_model = data.table(t(p))
names(fitness_norm_model) = c(paste0("scale_", reps_num), paste0("shift_", reps_num))

#wild-type correction such that mean(wild-type) = 0
wt_corr = work_data[WT == T,rowMeans((.SD + unlist(fitness_norm_model[,.SD,, .SDcols = grep(paste0("shift_[", reps,"]"),names(fitness_norm_model))])) * 
                                      unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("scale_[",reps,"]"),names(fitness_norm_model))])),
                    ,.SDcols = grep(paste0("fitness[",reps,"]"),names(work_data))]
#normalize fitness values
for (j in seq_along(reps_num)) {
  work_data[all_reads ==T,paste0("fitness",reps_num[j]) := (.SD + unlist(fitness_norm_model[,.SD,,.SDcols = paste0("shift_",reps_num[j])])) * unlist(fitness_norm_model[,.SD,,.SDcols = paste0("scale_",reps_num[j])]) - wt_corr,
            ,.SDcols = paste0("fitness",reps_num[j])]
}

#calculate count-based error for each variant and each replicate
for (j in as.numeric(strsplit(reps,"")[[1]])) {
  wt_corr = as.numeric(work_data[WT==T,1/.SD[,2] + 1/.SD[,1],,
                                 .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))])
  
  work_data[,paste0("cbe",j) := sqrt(unlist(fitness_norm_model[,.SD,.SDcols=paste0("scale_",j)])) * sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
            .SDcols = c(grep(paste0("^input",j,"$"),names(work_data)),grep(paste0("^output",j,"$"),names(work_data)))]
}

###########################################################################
### error model with replicate specific additive and multiplicative errors

#note fitness replication factor
Fcorr = unlist(fitness_norm_model[,.SD,.SDcols=grep(paste0("scale_[",reps,"]"),names(fitness_norm_model))])

#estimate additve and multiplicative factors
parameters = error_model_mult_inout_rep(DT = work_data,
                                             reps = reps,
                                             Ncores = 1,
                                             Fcorr = Fcorr)
#  table with error model parameters
error_model = data.table(parameter = rep(c("input","output","reperror"),each = 3),
                             rep = rep(reps_num, 3),
                             mean_value = colMeans(parameters, na.rm = T),
                             sd_value = apply(parameters, 2, sd, na.rm = T),
                             ensemble = sum(!is.na(parameters[,1])),
                             which_reps = reps)


#fitness and estimated error from full error model (MioA - Multiplicative Input/Output terms + additive error terms)
for (j in seq_along(reps_num)) {
	Corr = matrix(unlist(fitness_norm_model[, .SD, .SDcols=paste0("scale_", reps_num[j])]),ncol = 1, nrow = nrow(work_data))
	work_data[,paste0("error_MioA",reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input","output") & rep == reps_num[j], mean_value]),
														nrow = .N, ncol = 2, byrow = T) / .SD) + 
                                                                   matrix(error_model[parameter %in% c("reperror") & rep == reps_num[j], mean_value], 
                                                                   	nrow = .N, ncol = 1, byrow = T)),,
              .SDcols = c(grep(paste0("^input",reps_num[j],"$"),names(work_data)),
                          grep(paste0("^output",reps_num[j],"$"),names(work_data)))]
}
NTreps = 3
#merged fitness values
work_data[, fitness := rowSums(.SD[, 1 : NTreps] / (.SD[, (NTreps + 1) : (2 * NTreps)]^2), na.rm = T) / 
              rowSums(1 / (.SD[, (NTreps + 1) : (2 * NTreps)]^2), na.rm = T),
            ,.SDcols = c(grep(paste0("^fitness[", reps, "]$"), names(work_data)),
                         grep(paste0("^error_MioA[", reps, "]$"), names(work_data)))]

#and merged error
work_data[,error := sqrt(1 / rowSums(1 / .SD^2)),,
            .SDcols = grep(paste0("^error_MioA[", reps, "]$"), names(work_data))]

work_data[Nmut == 0, id1 := Mut]
work_data[Nmut == 1, id1 := paste0(strsplit(Mut, "")[[1]][2 : nchar(Mut)], collapse = ""), Mut]
work_data[Nmut == 1, WT_AA1 := paste0(strsplit(Mut, "")[[1]][1], collapse = ""), Mut]
work_data[Nmut == 1, Mut1 := paste0(strsplit(Mut, "")[[1]][nchar(Mut)], collapse = ""), Mut]
work_data[Nmut == 2, tmp := strsplit(Mut, "-")[[1]][1], Mut]
work_data[Nmut == 2, id1 := paste0(strsplit(tmp, "")[[1]][2 : nchar(tmp)], collapse = ""), tmp]
work_data[Nmut == 2, WT_AA1 := paste0(strsplit(tmp, "")[[1]][1], collapse = ""), tmp]
work_data[Nmut == 2, Mut1 := paste0(strsplit(tmp, "")[[1]][nchar(tmp)], collapse = ""), tmp]
work_data[Nmut == 2, tmp := strsplit(Mut, "-")[[1]][2], Mut]
work_data[Nmut == 2, id2 := paste0(strsplit(tmp, "")[[1]][2 : nchar(tmp)], collapse = ""), tmp]
work_data[Nmut == 2, WT_AA2 := paste0(strsplit(tmp, "")[[1]][1], collapse = ""), tmp]
work_data[Nmut == 2, Mut2 := paste0(strsplit(tmp, "")[[1]][nchar(tmp)], collapse = ""), tmp]
work_data[, tmp := NULL]

#adjust position variable (+1)
work_data[Nmut ==0, Pos1 := 0]
work_data[Nmut > 0, Pos1 := as.integer(gsub("[A-Z]", "", id1)) + 1, id1]
work_data[Nmut == 2, Pos2 := as.integer(gsub("[A-Z]", "", id2)) + 1, id2]
#redo id's
# work_data[Nmut == 0, id1 := "OX"]
work_data[Nmut > 0, id1 := paste0(Pos1, Mut1, collapse = ""), .(Pos1, Mut1)]
work_data[Nmut == 2, id2 := paste0(Pos2, Mut2, collapse = ""), .(Pos2, Mut2)]

all_data <- work_data[,.(Nmut, 
							Pos1, Pos2,
							id1, id2,
              mean_count = input1 + input2 + input3,
							b_fitness = fitness,
							b_sigma = error)]

ggplot(all_data, aes(mean_count, b_fitness)) +
  geom_hex() +
  scale_x_log10() + 
  scale_fill_continuous(trans = "log")
ggsave("results/preprocessing/GB1_input_fitness.pdf")  

ggplot(all_data, aes(b_fitness, color = factor(Nmut), linetype = mean_count > 10^4)) +
  geom_density(adjust = 0.5) +
  geom_vline(xintercept = c(-5.69, -5.09, 0, 0.12))
ggsave("results/preprocessing/GB1_b_fitness_distribution.pdf")


ggplot(work_data, aes(fitness1, color = factor(Nmut), linetype = input1 > 1/3*10^4)) +
  geom_density(adjust = 0.5) +
  geom_vline(xintercept = c(-5.69, -5.09, 0, 0.12))
ggsave("results/preprocessing/GB1_b_fitness_distribution.pdf")

 
ggplot(all_data, aes(b_fitness,b_sigma)) +
  geom_hex() +
  scale_y_log10() +
  scale_fill_continuous(trans="log")
ggsave("results/preprocessing/GB1_b_fitness_error.pdf")


all_data[, id := ifelse(Nmut == 2, paste0(id1, "_", id2), id1), .(id1, id2)]

all_data[Nmut == 1, .N, Pos1][order(Pos1)]
#this is missing positions 1 and 2 which showed weird behaviour between replicates

set.seed(1603)
all_data[, test_set := Nmut == 2 & runif(1) <= 0.01,id]
all_data[, .N, .(test_set, Nmut)]

# save
write.table(all_data,
	file = paste0("processed_data/GB1_dG_dataset.txt"),
	row.names = F,
	quote = F
)



######################################
### derive from DMS2struct project ###
######################################

GB1_wildtype = fread("dataset/GB1/DMS_wildtype.txt")
GB1_singles = fread("dataset/GB1/DMS_singles.txt")
GB1_doubles = fread("dataset/GB1/DMS_doubles.txt")

#regularization for doubles
GB1_doubles_reg = fread("dataset/GB1/DMS_doubles_cond.txt")

GB1_singles[, id1 := paste0(Pos, Mut, collapse = ""), .(Pos, Mut)]

GB1_doubles[, id1 := paste0(Pos1, Mut1, collapse = ""), .(Pos1, Mut1)]
GB1_doubles[, id2 := paste0(Pos2, Mut2, collapse = ""), .(Pos2, Mut2)]


GB1_doubles_reg[, id1 := paste0(Pos1, Mut1, collapse = ""), .(Pos1, Mut1)]
GB1_doubles_reg[, id2 := paste0(Pos2, Mut2, collapse = ""), .(Pos2, Mut2)]


GB1 = rbind(GB1_wildtype[,
                      .(Nmut = 0,
                        Pos1 = 0, Pos2 = NA,
                        id1 = "0X", id2 = NA, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))],
          GB1_singles[is.fitness == T,
                      .(Nmut = 1,
                        Pos1 = Pos, Pos2 = NA,
                        id1, id2 = NA, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))],
      GB1_doubles[is.fitness == T,
                      .(Nmut = 2,
                        Pos1, Pos2,
                        id1, id2, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))])

GB1[, id := ifelse(Nmut == 2, paste0(id1, "_", id2), id1), .(id1, id2)]

set.seed(1603)
GB1[, test_set := Nmut == 2 & runif(1) <= 0.01,id]
GB1[, .N, .(test_set, Nmut)]

write.table(GB1,
  file = paste0("processed_data/GB1_dG_dataset2.txt"),
  row.names = F,
  quote = F
)



GB1_reg = rbind(GB1_wildtype[,
                      .(Nmut = 0,
                        Pos1 = 0, Pos2 = NA,
                        id1 = "0X", id2 = NA, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))],
          GB1_singles[is.fitness == T,
                      .(Nmut = 1,
                        Pos1 = Pos, Pos2 = NA,
                        id1, id2 = NA, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))],
      GB1_doubles_reg[is.fitness == T,
                      .(Nmut = 2,
                        Pos1, Pos2,
                        id1, id2, 
                        mean_count = count_r1_t0,
                        b_fitness = fitness, 
                        b_sigma = sqrt(sigma^2 + 0.01^2))])

GB1_reg[, id := ifelse(Nmut == 2, paste0(id1, "_", id2), id1), .(id1, id2)]

set.seed(1603)
GB1_reg[, test_set := Nmut == 2 & runif(1) <= 0.01,id]
GB1_reg[, .N, .(test_set, Nmut)]

write.table(GB1_reg,
  file = paste0("processed_data/GB1_dG_dataset2reg.txt"),
  row.names = F,
  quote = F
)




#############################
### structural properties ###
#############################

structural_data <- unique(GB1_singles[,.(id = id1, 
  Pos, WT_AA, Mut)])
setkey(structural_data, Pos, Mut)
GB1_seq <- scan("dataset/GB1/GB1_sequence.fasta", what = "character")[2]
GB1_1FCCseq <- "TTYKLVINGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
FC_seq <- "PSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPQVKFNWYVDGVQVHNAKTKPREQQYNSTYRVVSVLTVLHQNWLDGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSL"


## calculate contact maps & secondary structure from PDB file
pairdistances_from_PDB(input_file = "dataset/GB1/1pga.pdb",
                       dataset_dir = "",
                       aa_seq = GB1_seq)

# trans contact map with IgG-FC
pairdistances_from_PDB(input_file = "dataset/GB1/1fcc.pdb",
                       dataset_dir = "",
                       given_chainids = c("C","A"),idx_DMS_start = c(1,1),idx_pdb_start = c(1,238),
                       aa_seq = list(GB1_1FCCseq,FC_seq))
PDB_1fcc_CA <- fread("processed_data/PDB_contactmap_1fcc_CA.txt")
PDB_1fcc_CA[,IgG_FC_id := paste0(Pos2, WT_AA2), .(Pos2, WT_AA2)]

structural_data[, HAmin_ligand := PDB_1fcc_CA[Pos1 %in% Pos, min(HAmin)], Pos]
structural_data[, scHAmin_ligand := PDB_1fcc_CA[Pos1 %in% Pos, min(scHAmin)], Pos]
structural_data[, ligand_proximal := PDB_1fcc_CA[Pos1 %in% Pos, IgG_FC_id[which.min(HAmin)]], Pos]

### RSA
RSA_1pga = fread("dataset/GB1/1pga.rsa",skip=8,nrows = 56)
names(RSA_1pga) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
structural_data = merge(structural_data,
  RSA_1pga[chain == "A", .(
    Pos,
    RSA_unbound = RSA_all_rel)],
  by = "Pos")

RSA_1fcc = fread("dataset/GB1/1fcc.rsa",skip=215,nrows = 56)
names(RSA_1fcc) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
RSA_1fcc
# add to singles data.table

structural_data = merge(structural_data,
  RSA_1fcc[chain == "C", .(
    Pos,
    RSA_bound = RSA_all_rel)],
  by = "Pos")

ggplot(structural_data, aes(RSA_bound, RSA_unbound, color = RSA_unbound < 10)) +
  geom_point(size = 3) +
  geom_abline()
ggsave("results/preprocessing/GB1_RSA.pdf")
#RSA similar between 1pga and 1fcc

threshold_RSA = 10
threshold_ligand_dist = 5 
structural_data[RSA_unbound <= threshold_RSA,type := "core"]
structural_data[RSA_unbound <= threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "core_bind"]
structural_data[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist,type := "surf_bind"]
structural_data[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist,type := "surface"]
structural_data[,type := factor(type,levels=c("core","surface","core_bind","surf_bind"))]
structural_data[,.(Npos = length(unique(Pos))),type]

### integrate Mayo dG data
Nishtal_sdG = fread("dataset/GB1/Nishtal2019_GB1_dG_stability.txt")
Nishtal_sdG[, variant := gsub("^[A-Z]0*", "", Description), Description]

structural_data[, s_ddg_measured := 
  Nishtal_sdG[variant == id & `Assay/Protocol` == "ddG(mAvg)_mean" , -Data],id]
structural_data[, s_ddg_measured_sd := 
  Nishtal_sdG[variant == id & `Assay/Protocol` == "SD of ddG(mAvg)_mean" , Data],id]

structural_data[!is.na(s_ddg_measured) & !is.na(s_ddg_measured_sd) & s_ddg_measured < 4, 
  s_ddg_measured_sdsmooth := loess(s_ddg_measured_sd ~ s_ddg_measured)$fitted]

structural_data[s_ddg_measured == 4, 
  s_ddg_measured_sdsmooth := structural_data[s_ddg_measured < 4, 
    s_ddg_measured_sdsmooth[which.max(s_ddg_measured)]]]

ggplot(structural_data) +
  geom_point(aes(s_ddg_measured, s_ddg_measured_sd)) + 
  geom_line(aes(s_ddg_measured, s_ddg_measured_sdsmooth), color = 'red') +
  geom_smooth(aes(s_ddg_measured, s_ddg_measured_sd))

ggsave("results/preprocessing/GB1_sddg_measured_Nishtal2019.pdf")

write.table(x = structural_data,
  file = "processed_data/GB1_structural_data.txt",
  col.names = T,
  row.names = F,
  quote = F)
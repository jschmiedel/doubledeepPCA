################################################################
############ data preprocessing for GRB2-SH3 domain ############
################################################################


############################################################################
############ import and process DLG4-PDZ3 doubledeepPCA datasets ###########
############################################################################

setwd("doubledeepPCA")

# create directory structure
dir.create("dg_models/PDZ3/", showWarnings = FALSE)
dir.create("dg_models/PDZ3/data", showWarnings = FALSE)
dir.create("dg_models/PDZ3/data/preprocessing", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
require(ggpubr)
theme_set(theme_bw(base_size = 9))

########################################################
### from intermediate DiMSum files before aa merging ###
########################################################

############################################
#### second nicking mutagenesis library ####
############################################
# one round of nicking mutagensis on a pool of 10 different backgrounds

##########################
#### NM2 abundancePCA ####
##########################

all_variants = fread("dataset/DLG4-PDZ3/01e-PDZ_NM2_stabilityPCA/01e-PDZ_NM2_stabilityPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(4.81, 5.31, 4.83)

#error model parameters tmp/errormodel.txt
mult_in = c(3.53, 2.35, 1.01)
mult_out = c(1.04, 2.46, 1.04)
add_rep = c(0.108, 0.082, 1e-4)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]

### identify backgrounds
backgrounds = c("N53R", "G14V", "P1G",  "S10I", "L13R", "S40A", "V5R",  "I49L", "G46A")
backgrounds_codon = c("1ggt","5agg","10att","13cgt","14gtt","40gcg","46gcg","49ctg","53agg")

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
wt_AA_seq <- all_variants[WT == T, aa_seq]
wt_AA_seq_split <- strsplit(wt_AA_seq,"")[[1]]

all_variants[Nham_aa == Nmut_codons, Pos1 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[1], aa_seq]
all_variants[Nham_aa == Nmut_codons, Mut1 := strsplit(aa_seq, "")[[1]][Pos1], aa_seq]
all_variants[Nham_aa == Nmut_codons, Codon1 := paste0(strsplit(nt_seq, "")[[1]][1:3 + (Pos1-1)*3], collapse=""), nt_seq]
all_variants[Nham_aa == Nmut_codons, Nham_nt1 := sum((strsplit(nt_seq, "")[[1]][1:3 + (Pos1-1)*3] != wt_NT_seq_split[1:3 + (Pos1-1)*3])), nt_seq]
all_variants[Nham_aa == Nmut_codons, idc1 := paste0(Pos1, Codon1, collapse = ""), .(Pos1, Codon1)]

all_variants[Nham_aa == 2 & Nmut_codons == 2, Pos2 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[2], aa_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, Mut2 := strsplit(aa_seq, "")[[1]][Pos2], aa_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, Codon2 := paste0(strsplit(nt_seq, "")[[1]][1:3 + (Pos2-1)*3], collapse=""), nt_seq]
all_variants[Nham_aa == Nmut_codons, Nham_nt2 := sum((strsplit(nt_seq, "")[[1]][1:3 + (Pos2-1)*3] != wt_NT_seq_split[1:3 + (Pos2-1)*3])), nt_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, idc2 := paste0(Pos2, Codon2, collapse = ""), .(Pos2, Codon2)]

all_variants[, is.background := !is.na(WT) | idc1 %in% backgrounds_codon & is.na(idc2)]
all_variants[is.background == T]

all_variants[WT == T, background := "WT"]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & is.background, background := idc1]
all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, background := "WT"]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Pos := Pos1]
all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Mut := Mut1]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Nham_nt_2bg := Nham_nt1]

all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, background := idc1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Pos := Pos2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Mut := Mut2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Nham_nt_2bg := Nham_nt2]

all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, background := idc2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Pos := Pos1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Mut := Mut1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Nham_nt_2bg := Nham_nt1]

all_variants[, .N, background]
all_variants[, .N, Nham_nt_2bg]
# >50k variants that don't fit into a background

all_variants[is.na(background)][order(mean_count)]
#mark synonymous variants
all_variants[Nham_aa == 0 & Nmut_codons == 1, is.synonymous := T]
all_variants[, .N, is.synonymous]

all_variants[is.na(background) & is.na(is.synonymous)][order(-mean_count)][1:20]
# mostly sequencing errors (only 40 above mean_count == 10) or synonymous variants with Nmut_codons >= 2 >>> kick out

all_variants <- all_variants[!is.na(background) | !is.na(is.synonymous)]

#check against background mean count
all_variants[background != "WT", background_mc := all_variants[is.background == T & idc1 %in% unlist(.SD), mean_count], background, .SDcols = "background"]
all_variants[background == "WT", background_mc := all_variants[WT == T, mean_count]]

theme_set(theme_bw(base_size = 11))
p <- ggplot(all_variants[background_mc > 1e3], aes(mean_count/background_mc, color = factor(Nham_nt_2bg))) +
  geom_density() +
  scale_x_log10()
ggsave(plot = p, file = "dg_models/PDZ3/data/preprocessing/PDZ_NM2_abundancePCA_Nham_nt_2bg_VSmean_counts.pdf", width = 8, height = 6.5)


# # single AA variants
nt1 <- all_variants[Nham_nt_2bg == 1,
  .(background, Pos, Mut, count_nt1 = log10(mean_count/background_mc), fitness_nt1 = fitness)]
nt2 <- all_variants[Nham_nt_2bg == 2,
  .(background, Pos, Mut, count_nt2 = log10(mean_count/background_mc), fitness_nt2 = fitness)]
nt3 <- all_variants[Nham_nt_2bg == 3,
  .(background, Pos, Mut, count_nt3 = log10(mean_count/background_mc), fitness_nt3 = fitness)]

X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
# X <- merge(nt1, nt2)
theme_set(theme_bw())
p <- ggplot(X, aes(fitness_nt2, fitness_nt1, color = count_nt1 > -3.5 & count_nt2 > -4)) +
  geom_point() +
  facet_wrap(background ~ .) +
  stat_cor(method = "pearson")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_abundancePCA_Nham_nt_count_fitness_dependencies.pdf", width = 6.5, height = 4)

all_variants[,.N, .(Nham_nt_2bg, freq = mean_count > 10 &  (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) | (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4))][order(Nham_nt_2bg, freq)]
all_variants[,.N, .(freq = (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) | (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4.5))][order(freq)]


#how many mutations per background?
all_variants[,.N, background][order(N)]
all_variants[is.background == T, .(background, mean_count)][order(mean_count)]
# >> depends on frequency of background in nicking mutagenesis
ggplot(merge(all_variants[,.(mutations_per_background = .N), background],
              all_variants[is.background == T, .(background, mean_count)]),
        aes(mean_count, mutations_per_background)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "read counts of background", y = "# variants contributed to final dataset")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_abundancePCA_background_mutations_meancount.pdf", width = 5, height = 5)

#how many backgrounds per mutation?
all_variants[,.N, .(Pos, Mut)][,.(backgrounds_per_mut = .N),N][order(N)]
all_variants[,.N, .(Pos, Mut)][,.(backgrounds_per_mut = .N),N][,sum(backgrounds_per_mut)]
# most mutations appear at least onces
ggplot(all_variants[,.N, .(Pos, Mut)],aes(N)) +
  geom_histogram(bins = 10) +
  scale_x_continuous(breaks = seq(1:10)) +
  labs(x = "backgrounds per variant", y = "# variants")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_abundancePCA_mutations_per_background.pdf", width = 5, height = 4)

## input reads against fitness | backgrounds
p <- ggplot(all_variants, aes(mean_count, fitness)) +
  geom_hex() +
  scale_x_log10() +
  scale_fill_continuous(trans = "log10") +
  facet_wrap(background ~ .)
ggsave(plot = p, file = "dg_models/PDZ3/data/preprocessing/PDZ_NM2_abundancePCA_meancount_fitness_background.pdf")


## gather synonymous variants and cut down to fitness, sigma and aa_seq columns
all_variants <- all_variants[mean_count > 10 & (is.background == T |
    (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) | (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4.5)) &
    STOP == FALSE & STOP_readthrough == FALSE,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), .(aa_seq)
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01e-PDZ3_NM2_abundancePCA_thresholded.RData")








########################
#### NM2 bindingPCA ####
########################

all_variants = fread("dataset/DLG4-PDZ3/01f-PDZ_NM2_bindingPCA/01f-PDZ_NM2_bindingPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5.07, 4.93, 4.73)

#error model parameters
mult_in = c(1.21, 1.25, 1.08)
mult_out = c(1, 1.01, 1.05)
add_rep = c(0.0173, 1e-4, 0.056)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]

### identify backgrounds
backgrounds = c("N53R", "G14V", "P1G",  "S10I", "L13R", "S40A", "V5R",  "I49L", "G46A")
backgrounds_codon = c("1ggt","5agg","10att","13cgt","14gtt","40gcg","46gcg","49ctg","53agg")

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
wt_AA_seq <- all_variants[WT == T, aa_seq]
wt_AA_seq_split <- strsplit(wt_AA_seq,"")[[1]]

all_variants[Nham_aa == Nmut_codons, Pos1 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[1], aa_seq]
all_variants[Nham_aa == Nmut_codons, Mut1 := strsplit(aa_seq, "")[[1]][Pos1], aa_seq]
all_variants[Nham_aa == Nmut_codons, Codon1 := paste0(strsplit(nt_seq, "")[[1]][1:3 + (Pos1-1)*3], collapse=""), nt_seq]
all_variants[Nham_aa == Nmut_codons, Nham_nt1 := sum((strsplit(nt_seq, "")[[1]][1:3 + (Pos1-1)*3] != wt_NT_seq_split[1:3 + (Pos1-1)*3])), nt_seq]
all_variants[Nham_aa == Nmut_codons, idc1 := paste0(Pos1, Codon1, collapse = ""), .(Pos1, Codon1)]

all_variants[Nham_aa == 2 & Nmut_codons == 2, Pos2 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[2], aa_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, Mut2 := strsplit(aa_seq, "")[[1]][Pos2], aa_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, Codon2 := paste0(strsplit(nt_seq, "")[[1]][1:3 + (Pos2-1)*3], collapse=""), nt_seq]
all_variants[Nham_aa == Nmut_codons, Nham_nt2 := sum((strsplit(nt_seq, "")[[1]][1:3 + (Pos2-1)*3] != wt_NT_seq_split[1:3 + (Pos2-1)*3])), nt_seq]
all_variants[Nham_aa == 2 & Nmut_codons == 2, idc2 := paste0(Pos2, Codon2, collapse = ""), .(Pos2, Codon2)]

all_variants[, is.background := !is.na(WT) | idc1 %in% backgrounds_codon & is.na(idc2)]
all_variants[is.background == T]

all_variants[WT == T, background := "WT"]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & is.background, background := idc1]
all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, background := "WT"]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Pos := Pos1]
all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Mut := Mut1]

all_variants[Nham_aa == 1 & Nmut_codons == 1 & !is.background, Nham_nt_2bg := Nham_nt1]

all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, background := idc1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Pos := Pos2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Mut := Mut2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc1 %in% backgrounds_codon, Nham_nt_2bg := Nham_nt2]

all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, background := idc2]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Pos := Pos1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Mut := Mut1]
all_variants[Nham_aa == 2 & Nmut_codons == 2 & idc2 %in% backgrounds_codon, Nham_nt_2bg := Nham_nt1]

all_variants[, .N, background]
all_variants[, .N, Nham_nt_2bg]
# >70k variants that don't fit into a background

all_variants[is.na(background)][order(mean_count)]
#mark synonymous variants
all_variants[Nham_aa == 0 & Nmut_codons == 1, is.synonymous := T]
all_variants[, .N, is.synonymous]

all_variants[is.na(background) & is.na(is.synonymous)][order(-mean_count)][1:20]
# mostly sequencing errors (only 40 above mean_count == 10) or synonymous variants with Nmut_codons >= 2 >>> kick out

all_variants <- all_variants[!is.na(background) | !is.na(is.synonymous)]

#check against background mean count
all_variants[background != "WT", background_mc := all_variants[is.background == T & idc1 %in% unlist(.SD), mean_count], background, .SDcols = "background"]
all_variants[background == "WT", background_mc := all_variants[WT == T, mean_count]]

theme_set(theme_bw(base_size = 11))
p <- ggplot(all_variants[background_mc > 1e3], aes(mean_count/background_mc, color = factor(Nham_nt_2bg))) +
  geom_density() +
  scale_x_log10()
ggsave(plot = p, file = "dg_models/PDZ3/data/preprocessing/PDZ_NM2_bindingPCA_Nham_nt_2bg_VSmean_counts.pdf", width = 8, height = 6.5)


# # single AA variants
nt1 <- all_variants[Nham_nt_2bg == 1,
  .(background, Pos, Mut, count_nt1 = log10(mean_count/background_mc), fitness_nt1 = fitness)]
nt2 <- all_variants[Nham_nt_2bg == 2,
  .(background, Pos, Mut, count_nt2 = log10(mean_count/background_mc), fitness_nt2 = fitness)]
nt3 <- all_variants[Nham_nt_2bg == 3,
  .(background, Pos, Mut, count_nt3 = log10(mean_count/background_mc), fitness_nt3 = fitness)]

X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
# X <- merge(nt1, nt2)
theme_set(theme_bw())
p <- ggplot(X, aes(fitness_nt2, fitness_nt1, color = count_nt1 > -3.5 & count_nt2 > -4)) +
  geom_point() +
  facet_wrap(background ~ .) +
  stat_cor(method = "pearson")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_bindingPCA_Nham_nt_count_fitness_dependencies.pdf",
  width = 6.5,
  height = 4)

all_variants[,.N, .(Nham_nt_2bg, freq = mean_count > 10 &  (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) | (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4))][order(Nham_nt_2bg, freq)]
all_variants[,.N, .(freq = (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) | (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4.5))][order(freq)]

#how many mutations per background?
all_variants[,.N, background][order(N)]
all_variants[is.background == T, .(background, mean_count)][order(mean_count)]
# >> depends on frequency of background in nicking mutagenesis
ggplot(merge(all_variants[,.(mutations_per_background = .N), background],
              all_variants[is.background == T, .(background, mean_count)]),
        aes(mean_count, mutations_per_background)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "read counts of background", y = "# variants contributed to final dataset")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_bindingPCA_background_mutations_meancount.pdf", width = 5, height = 5)

#how many backgrounds per mutation?
all_variants[,.N, .(Pos, Mut)][,.(backgrounds_per_mut = .N),N][order(N)]
all_variants[,.N, .(Pos, Mut)][,.(backgrounds_per_mut = .N),N][,sum(backgrounds_per_mut)]
# most mutations appear at least onces
ggplot(all_variants[,.N, .(Pos, Mut)],aes(N)) +
  geom_histogram(bins = 10) +
  scale_x_continuous(breaks = seq(1:10)) +
  labs(x = "backgrounds per variant", y = "# variants")
ggsave("dg_models/PDZ3/data/preprocessing/PDZ_NM2_bindingPCA_mutations_per_background.pdf", width = 5, height = 4)

## input reads against fitness | backgrounds
p <- ggplot(all_variants, aes(mean_count, fitness)) +
  geom_hex() +
  scale_x_log10() +
  scale_fill_continuous(trans = "log10") +
  facet_wrap(background ~ .)
ggsave(plot = p, file = "dg_models/PDZ3/data/preprocessing/PDZ_NM2_bindingPCA_meancount_fitness_background.pdf")


## gather synonymous variants and cut down to fitness, sigma and aa_seq columns
all_variants <- all_variants[mean_count > 10 & (is.background == T |
    (Nham_nt_2bg  == 1 & log10(mean_count/background_mc) > -3.5) |
    (Nham_nt_2bg > 1 & log10(mean_count/background_mc) > -4.5)) &
    STOP == FALSE & STOP_readthrough == FALSE,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), .(aa_seq)
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01f-PDZ3_NM2_bindingPCA_thresholded.RData")




###########################################
#### first nicking mutagenesis library ####
###########################################
# one round of nicking mutagensis on wildtype sequence

##########################
#### NM1 abundancePCA ####
##########################

all_variants = fread("dataset/DLG4-PDZ3/01c-PDZ_NM_stabilityPCA/01c-PDZ_NM_stabilityPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(4.8, 5, 5.1)

#error model parameters
mult_in = c(1.97, 1.77, 1.5)
mult_out = c(1, 1.12, 1.14)
add_rep = c(0.0062, 0.0008, 0.0021)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]


all_variants[, .N,. (Nham_aa, Nham_nt)][order(Nham_aa, Nham_nt)]

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
all_variants[, Nham_codons := length(unique(ceiling((which(strsplit(nt_seq,"")[[1]] != wt_NT_seq_split)) / 3))), nt_seq]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1 & Nham_codons == 1,
    .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2 & Nham_codons == 1,
    .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3 & Nham_codons == 1,
    .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

X = merge(merge(nt1,nt2,all=T),nt3,all=T)
# ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 2 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)) & (count_nt3 > 1 | is.na(count_nt3))))
ggsave("dg_models/PDZ3/data/preprocessing/NM1_abundancePCA_Nhamnt_count_fitness_dependencies.pdf",
  width = 6,
  height = 6)
#if Nham_nt == 1, variants below 10^2 read counts are sequencing errors
#if Nham_nt > 1, variants below 10^1 read counts are sequencing errors
#also there's a proportionality between readcounts of the same AA variants that could be taken into account


# ## merge based on hamming-distance specific read thresholds
all_data_rest = all_variants[Nham_aa == 1 & Nham_codons == 1 & ((Nham_nt == 1 & mean_count > 100) | (Nham_nt > 1 & mean_count > 10)),
             .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
               sigma = sqrt(1/sum(sigma^-2)),
               mean_count = sum(mean_count),
               Nsyn = .N,
               Nham_aa = unique(Nham_aa),
               WT = unique(WT),
               STOP = unique(STOP),
               STOP_readthrough=unique(STOP_readthrough)),
             aa_seq]
all_data = all_variants[Nham_aa == 1 & Nham_codons == 1 & mean_count > 10,
             .(fitness = sum(fitness * sigma^-2)/sum(sigma^-2),
               sigma = sqrt(1/sum(sigma^-2)),
               mean_count = sum(mean_count),
               Nsyn = .N,
               Nham_aa = unique(Nham_aa),
               WT = unique(WT),
               STOP = unique(STOP),
               STOP_readthrough=unique(STOP_readthrough)),
             aa_seq]

X = merge(all_data_rest[, .(aa_seq, fitness_rest = fitness, Nsyn_rest = Nsyn)],
      all_data[, .(aa_seq, fitness, Nsyn)])
ggplot(X[Nsyn > 1],aes(fitness,fitness_rest)) +
  geom_point()
ggsave("dg_models/PDZ3/data/preprocessing/NM1_abundancePCA_Nham_nt_specific_thresholds.pdf",
  width = 6,
  height = 6)
#not a huge difference in total, but for ~50 variants
X[abs(fitness - fitness_rest) > 0.05]

all_variants <- all_variants[
  Nham_aa == 1 & Nham_codons == 1 &
  ((Nham_nt == 1 & mean_count > 100) |
    (Nham_nt > 1 & mean_count > 10)),
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01c-PDZ3_NM1_abundancePCA_thresholded.RData")

########################
#### NM1 bindingPCA ####
########################

all_variants = fread("dataset/DLG4-PDZ3/01d-PDZ_NM_bindingPCA/01d-PDZ_NM_bindingPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5.1, 5.3, 5.2)

#error model parameters
mult_in = c(1, 1, 1)
mult_out = c(1.17, 1, 1)
add_rep = c(0.0014, 0.0055, 0.0046)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]


all_variants[, .N,. (Nham_aa, Nham_nt)][order(Nham_aa, Nham_nt)]

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
all_variants[, Nham_codons := length(unique(ceiling((which(strsplit(nt_seq,"")[[1]] != wt_NT_seq_split)) / 3))), nt_seq]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1 & Nham_codons == 1,
    .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2 & Nham_codons == 1,
    .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3 & Nham_codons == 1,
    .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

X = merge(merge(nt1,nt2,all=T),nt3,all=T)
# ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("fitness",names(X)),
        aes(alpha=0.5,color=(count_nt1 > 2 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)) & (count_nt3 > 1 | is.na(count_nt3))))
ggsave("dg_models/PDZ3/data/preprocessing/NM1_bindingPCA_Nhamnt_count_fitness_dependencies.pdf",
  width = 6,
  height = 6)
#if Nham_nt == 1, variants below 10^2 read counts are sequencing errors
#if Nham_nt > 1, variants below 10^1 read counts are sequencing errors
#also there's a proportionality between readcounts of the same AA variants that could be taken into account

all_variants <- all_variants[
  Nham_aa == 1 & Nham_codons == 1 &
  ((Nham_nt == 1 & mean_count > 100) |
    (Nham_nt > 1 & mean_count > 10)),
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01d-PDZ3_NM1_bindingPCA_thresholded.RData")



################################
#### errorprone PCR library ####
################################

############################
#### epPCR abundancePCA ####
############################

all_variants = fread("dataset/DLG4-PDZ3/01a-PDZ_epPCR_stabilityPCA/01a-PDZ_epPCR_stabilityPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5, 4.8, 5)

#error model parameters
mult_in = c(1.03, 4.16, 1.02)
mult_out = c(1.14, 1, 1.41)
add_rep = c(0.0054, 0.00045, 0.008)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]


all_variants[, .N,. (Nham_aa, Nham_nt)][order(Nham_aa, Nham_nt)]

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
all_variants[, Nham_codons := length(unique(ceiling((which(strsplit(nt_seq,"")[[1]] != wt_NT_seq_split)) / 3))), nt_seq]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1 & Nham_codons == 1,
    .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2 & Nham_codons == 1,
    .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3 & Nham_codons == 1,
    .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

X = merge(merge(nt1,nt2,all=T),nt3,all=T)
# ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("count",names(X)))
ggsave("dg_models/PDZ3/data/preprocessing/epPCR_abundancePCA_Nhamnt_count_fitness_dependencies.pdf",
  width = 6,
  height = 6)
# there is not Nham_nt>1 variants with sufficient read counts

all_variants <- all_variants[
  Nham_aa == 1 & Nham_nt == 1 & Nham_codons == 1,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01a-PDZ3_epPCR_abundancePCA_thresholded.RData")


##########################
#### epPCR bindingPCA ####
##########################

all_variants = fread("dataset/DLG4-PDZ3/01b-PDZ_epPCR_bindingPCA/01b-PDZ_epPCR_bindingPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5, 5, 4.8)

#error model parameters
mult_in = c(1, 1, 1)
mult_out = c(1, 1.12, 1.25)
add_rep = c(1e-4, 0.0149, 0.0143)

for (i in 1:3) {
  correction_factor_reads =  all_variants[,2^gen[i] * sum(.SD[, 1]) / sum(.SD[, 2]),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  wt_corr = all_variants[WT == T, (log(unlist(.SD[, 2]/.SD[, 1])) + log(correction_factor_reads)),
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]

  all_variants[, paste0("fitness", i) := (log(unlist((.SD[, 2] + 0.5)/(.SD[, 1] + 0.5))) +
    log(correction_factor_reads)) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
  all_variants[, paste0("sigma", i) := sqrt(mult_in[i] / (.SD[, 1] + 0.5) +
    mult_out[i] / (.SD[, 2] + 0.5) + add_rep[i]) / wt_corr,
    .SDcols = c(paste0("input",i,"_e",i,"_s0_bNA_count"), paste0("output", i, "_e",i,"_s1_b1_count"))]
}
# merge fitness, error and input counts
all_variants[, fitness := rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD), ,
   .SDcols = grep("input[123]", names(all_variants))]


all_variants[, .N,. (Nham_aa, Nham_nt)][order(Nham_aa, Nham_nt)]

wt_NT_seq <- all_variants[WT == T, nt_seq]
wt_NT_seq_split <- strsplit(wt_NT_seq,"")[[1]]
all_variants[, Nham_codons := length(unique(ceiling((which(strsplit(nt_seq,"")[[1]] != wt_NT_seq_split)) / 3))), nt_seq]

nt1 = all_variants[Nham_aa==1 & Nham_nt ==1 & Nham_codons == 1,
    .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
nt2 = all_variants[Nham_aa==1 & Nham_nt ==2 & Nham_codons == 1,
    .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
nt3 = all_variants[Nham_aa==1 & Nham_nt ==3 & Nham_codons == 1,
    .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

X = merge(merge(nt1,nt2,all=T),nt3,all=T)
# ggpairs(X,columns = grep("count",names(X)))
ggpairs(X,columns = grep("count",names(X)))
ggsave("dg_models/PDZ3/data/preprocessing/epPCR_bindingPCA_Nhamnt_count_fitness_dependencies.pdf",
  width = 6,
  height = 6)
# there is not Nham_nt>1 variants with sufficient read counts

all_variants <- all_variants[
  Nham_aa == 1 & Nham_nt == 1 & Nham_codons == 1,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq
]

## write to file
save(all_variants, file = "dg_models/PDZ3/data/01b-PDZ3_epPCR_bindingPCA_thresholded.RData")










##################################################
### extract structural data from PDB structure ###
##################################################

### compute contact map between PDZ3 domain and ligand
## source functions
filelist <- list.files("functions/")
invisible(sapply(paste0("functions/", filelist), source, .GlobalEnv))

pairdistances_from_PDB(
  input_file = "dataset/DLG4-PDZ3/PDB/pdb1be9.ent", dataset_dir = "",
  given_chainids = c("A", "B"),
  aa_seq = list(paste0(unique(singles[, .(Pos, WT_AA)])$WT_AA, collapse = ""), "KQTSV"),
  idx_pdb_start = c(311, 5),
  idx_DMS_start = c(1, 1), debug_this = F
)
contactmap_AB <- fread("processed_data/PDZ/PDB_contactmap_pdb1be9_AB.txt")
contactmap_AB[, WTAAPos2 := paste0(WT_AA2, Pos2)]
contactmap_AB[, Pos := Pos1]

structural_properties <- SH3_GAB2_distances[,
  .(HAmin_ligand = min(HAmin),
    scHAmin_ligand = min(scHAmin),
    ligand_AA = WTAAPos2[which.min(HAmin)]
  ),
Pos]


### load RSA values
RSA_chainA <- fread("dataset/DLG4-PDZ3/pdb1be9_A.rsa", skip = 9)
names(RSA_chainA) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel", "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel", "RSA_pol_abs", "RSA_pos_rel")
RSA_chainAB <- fread("dataset/DLG4-PDZ3/pdb1be9_AB.rsa", skip = 9)
names(RSA_chainAB) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel", "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel", "RSA_pol_abs", "RSA_pos_rel")
# add to singles data.table
structural_properties <- merge(structural_properties,
  RSA_chainA[chain == "A", .(Pos = Pos - 310, RSA_unbound = RSA_all_rel)], by = "Pos", all.x = T)
structural_properties <- merge(structural_properties, R
  SA_chainAB[chain == "A", .(Pos = Pos - 310, RSA_bound = RSA_all_rel)], by = "Pos", all.x = T)
# >>> all residues that loose RSA when bound to ligand are on interface with ligand
ggplot(structural_properties, aes(RSA_unbound, RSA_bound, color = HAmin_ligand < 5)) + geom_point(size = 3)

## split data into core, surface and binding
ggplot(structural_properties, aes(RSA_unbound, s_fitness)) + geom_point() + geom_smooth()
# >>> set core threshodl as RSA <= 25
threshold_RSA <- 25
threshold_ligand_dist <- 5
structural_properties[RSA_unbound <= threshold_RSA, type := "core"]
structural_properties[RSA_unbound <= threshold_RSA & HAmin_ligand <= threshold_ligand_dist, type := "core_binding"]
structural_properties[RSA_unbound > threshold_RSA, type := "surface"]
structural_properties[RSA_unbound > threshold_RSA & HAmin_ligand <= threshold_ligand_dist, type := "surface_binding"]
structural_properties[, type := factor(type, levels = c("core", "surface", "core_binding", "surface_binding"))]
structural_properties[, .(Npos = length(unique(Pos))), type]

write.table(
  x = structural_properties,
  file = "dg_models/PDZ3/data/structural_properties.txt",
  row.names = F,
  quote = F
)

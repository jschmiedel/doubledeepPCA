################################################################
############ data preprocessing for GRB2-SH3 domain ############
################################################################




############################################################################
############ import and process GRB2-SH3 doubledeepPCA datasets ############
############################################################################

setwd("doubledeepPCA")

## source functions
filelist <- list.files("functions/")
invisible(sapply(paste0("functions/", filelist), source, .GlobalEnv))

# create directory structure
dir.create("dg_models/SH3/", showWarnings = FALSE)
dir.create("dg_models/SH3/data", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
theme_set(theme_bw(base_size = 9))

########################################################
### from intermediate DiMSum files before aa merging ###
########################################################

##################################################
#### nicking mutagenesis library abundancePCA ####
##################################################

all_variants = fread("dataset/GRB2-SH3/01c-GRB2_NM2_stabilityPCA/01c-GRB2_NM2_stabilityPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5, 5.1, 5.4)

#error model parameters tmp/errormodel.txt
mult_in = c(15.89, 12.46, 13.64)
mult_out = c(1, 1, 1)
add_rep = c(0.0074, 0.0028, 0.0014)

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


# # single AA variants
# nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1 & !STOP,
#   .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
# nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1 & !STOP,
#   .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
# nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1 & !STOP,
#   .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

# X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
# theme_set(theme_bw())
# ggpairs(X, columns = grep("count", names(X)))
# ggpairs(X,
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = (count_nt1 > 3 | is.na(count_nt1)) & (count_nt2 > 2 | is.na(count_nt2)) & (count_nt3 > 2 | is.na(count_nt3)))
# )
# ggsave("results/SH3/preprocessing/SH3_Nhamnt_count_fitness_dependencies_stabilityNM.pdf", width = 6, height = 6)
# # if Nham_nt==1, variants below 10^3 read counts are sequencing errors
# # if Nham_nt > 1, variants below 10^2 read counts are sequencing errors
# # also there's a proportionality between readcounts of the same AA variants that could be taken into account; this is probably from sequencing errors as well


# ## merge based on hamming-distance specific read thresholds
# all_data_rest <- all_variants[
#   Nham_aa == 1 & Nmut_codons == 1 & ((Nham_nt == 1 & mean_count > 10^3) | (Nham_nt > 1 & mean_count > 10^2)),
#   .(
#     fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
#     sigma = sqrt(1 / sum(sigma^-2)),
#     mean_count = sum(mean_count), Nsyn = .N,
#     Nham_aa = unique(Nham_aa), WT = unique(WT), STOP = unique(STOP), STOP_readthrough = unique(STOP_readthrough)
#   ), aa_seq
# ]
# all_data <- all_variants[Nham_aa == 1 & Nmut_codons == 1 & mean_count > 10, .(
#   fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
#   sigma = sqrt(1 / sum(sigma^-2)),
#   mean_count = sum(mean_count), Nsyn = .N,
#   Nham_aa = unique(Nham_aa), WT = unique(WT), STOP = unique(STOP), STOP_readthrough = unique(STOP_readthrough)
# ), aa_seq]
# X <- merge(
#   all_data_rest[, .(aa_seq, fitness_rest = fitness, Nsyn_rest = Nsyn)],
#   all_data[, .(aa_seq, fitness, Nsyn)]
# )
# ggplot(X[Nsyn > 1], aes(fitness, fitness_rest)) +
#   geom_point()
# # not a huge difference in total, but for ~30 variants
# X[abs(fitness - fitness_rest) > 0.05]

# ### double aa variants
# aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[1, 1])), fitness_nt2 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[1, 1])), fitness_nt3 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt4 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 4 & !STOP, .(count_nt4 = log10(unlist(.SD[1, 1])), fitness_nt4 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt5 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 5 & !STOP, .(count_nt5 = log10(unlist(.SD[1, 1])), fitness_nt5 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt6 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 6 & !STOP, .(count_nt6 = log10(unlist(.SD[1, 1])), fitness_nt6 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

# X <- merge(merge(merge(merge(aa2_nt2, aa2_nt3, all = T),
#   aa2_nt4,
#   all = T
# ),
# aa2_nt5,
# all = T
# ),
# aa2_nt6,
# all = T
# )

# theme_set(theme_bw())
# ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
# ggpairs(Y, columns = grep("count", names(X)))
# ggpairs(X[sample(x = .N, size = 10^5)],
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = ((count_nt2 > 2 | is.na(count_nt2)) &
#     (count_nt3 > 1.5 | is.na(count_nt3)) &
#     (count_nt4 > 1.5 | is.na(count_nt4)) &
#     (count_nt5 > 1.5 | is.na(count_nt5)) &
#     (count_nt6 > 1.5 | is.na(count_nt6))))
# )
# ggpairs(X[sample(x = .N, size = 10^4)],
#   columns = grep("fitness_nt[23]", names(X)),
#   aes(alpha = 0.5, color = (count_nt2 > 1 | is.na(count_nt2)) &
#     (count_nt3 > 1 | is.na(count_nt3)))
# )
# ggpairs(X[sample(x = .N, size = 10^4)],
#   columns = grep("fitness_nt[34]", names(X)),
#   aes(alpha = 0.5, color = (count_nt4 > 1.5 | is.na(count_nt4)) &
#     (count_nt3 > 1.5 | is.na(count_nt3)))
# )
# ggsave("results/SH3/preprocessing/SH3_doubles_Nhamnt_count_fitness_dependencies_stabilityNM.pdf", width = 6, height = 6)

# X2 <- merge(aa2_nt2, aa2_nt3)
# dt <- data.table(i = rep(seq(0, 3, 0.2), each = 16), j = seq(0, 3, 0.2), c = 0, n = 0)
# for (idx in seq(0, 3, 0.2)) {
#   for (jdx in seq(0, 3, 0.2)) {
#     dt[i == idx & j == jdx, c := X2[count_nt2 > (idx) & count_nt3 > (jdx), cor(fitness_nt2, fitness_nt3)]]
#     dt[i == idx & j == jdx, n := X2[count_nt2 > (idx) & count_nt3 > (jdx), .N]]
#   }
# }
# ggplot(dt[n > 100], aes(i, j, fill = (c))) + geom_raster() + scale_fill_gradient2(midpoint = 0.85, mid = "grey")
# ggplot(dt[n > 100], aes(i, j, fill = log10(n))) + geom_raster()
# # exclude variants below 10^2 reads for aa2_nt2
# # exclude viarants below 10^1.5 reads for all other

all_variants <- all_variants[
  (WT == TRUE |
  (Nham_aa == 1 & ((Nham_nt == 1 & mean_count > 10^3) | (Nham_nt > 1 & mean_count > 10^2))) |
  (Nham_aa == 2 & ((Nham_nt == 2 & mean_count > 10^2) | (Nham_nt > 2 & mean_count > 10^1.5)))) &
  STOP == FALSE & STOP_readthrough == FALSE,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq]

save(all_variants, file = "dg_models/SH3/data/01c-GRB2_NM2_stabilityPCA_aa012_thresholded.RData")


################################################
#### errorpronePCR library from abundancePCA ###
################################################

all_variants <- fread("dataset/GRB2-SH3/01a-GRB2_epPCR_stabilityPCA/01a-GRB2_epPCR_stabilityPCA_variant_data_merge.tsv")
# normalize by generations
gen = c(5.19, 5.07, 5.1)
time = rep(34.1, 3)
mult_in = c(2.23, 4.51, 4.35)
mult_out = c(1, 1, 1)
add_rep = c(0.022, 0.0005, 0.0042)
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
all_variants[, fitness :=
  rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, mean_count := rowMeans(.SD, na.rm = T), ,
  .SDcols = grep("input[123]", names(all_variants))]

# ### single aa variants
# nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1, .(count_nt1 = log10(unlist(.SD[, 1])), fitness_nt1 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

# X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
# ggpairs(X, columns = grep("count", names(X)))
# # there are no nt_3 variants with sufficient read counts
# ggpairs(X,
#   columns = grep("fitness_nt[12]", names(X)),
#   aes(alpha = 0.5, color = (count_nt1 > 1 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)))
# )
# # for nt_2 variants, it seems that only very high thresholds that leave no data contain real variants not swamped by sequencing errors
# # compare nt_1 vairants to NM data
# X <- merge(nt1, singles_NM_stabilityPCA[, .(aa_seq, fitness, mean_count = log10(mean_count))])
# ggpairs(X, columns = c("count_nt1", "mean_count"))
# ggpairs(X,
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = (count_nt1 > 3 | is.na(count_nt1)) & (mean_count > 3.5 | is.na(mean_count)))
# )
# # correlation is good not matter what the input read count of the epPCR library
# ## >> use all epPCR data from Nham_nt == 1

# ### double aa variants
# aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

# X <- merge(aa2_nt2, aa2_nt3)

# ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
# ggpairs(X,
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = ((count_nt2 > 2.5 | is.na(count_nt2))))
# )
# # nt==3 variants are just random
# # compare to NM data
# X <- merge(aa2_nt2, doubles_NM_stabilityPCA[, .(aa_seq, fitness, mean_count = log10(mean_count))])
# ggpairs(X,
#   columns = c("count_nt2", "mean_count"),
#   aes(alpha = 0.5, color = count_nt2 > 1.5)
# )
# ggpairs(X[count_nt2 > 1.5],
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = count_nt2 > 2)
# )

all_variants <- all_variants[
  (WT == TRUE |
  (Nham_aa == 1 & Nham_nt == 1) |
  (Nham_aa == 2 & Nham_nt == 2 & mean_count > 10^1.5)) &
  STOP == FALSE & STOP_readthrough == FALSE,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq]

save(all_variants, file = "dg_models/SH3/data/01a-GRB2_epPCR_stabilityPCA_aa012_thresholded.RData")



##############################################
#### errorpronePCR library from bindingPCA ###
##############################################

all_variants <- fread("dataset/GRB2-SH3/01b-GRB2_epPCR_bindingPCA/01b-GRB2_epPCR_bindingPCA_variant_data_merge.tsv")

# normalize by generations
gen = c(5.21, 5.26, 5.29)
time = rep(20.4, 3)
mult_in = c(1, 1, 1)
mult_out = c(1, 1, 1)
add_rep = c(0.00089, 0.00066, 0.0032)
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
all_variants[, fitness :=
  rowSums(.SD[, 1:3] / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T) /
  rowSums(1 / (.SD[, (3 + 1):(2 * 3)]^2), na.rm = T), ,
.SDcols = c(
  grep("^fitness[123]$", names(all_variants)),
  grep("^sigma[123]$", names(all_variants))
)
]
all_variants[, sigma := sqrt(1 / rowSums(1 / .SD^2, na.rm = T)), ,
  .SDcols = grep("^sigma[123]$", names(all_variants))
]
all_variants[, total_reads_input := input1_e1_s0_bNA_count + input2_e2_s0_bNA_count + input3_e3_s0_bNA_count]
all_variants[, total_reads_output := output1_e1_s1_b1_count + output2_e2_s1_b1_count + output3_e3_s1_b1_count]
all_variants[, mean_count := rowMeans(.SD), , .SDcols = grep("input[123]", names(all_variants))]

# ### single aa variants
# nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1, .(count_nt1 = log10(unlist(.SD[, 1])), fitness_nt1 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

# X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
# ggpairs(X, columns = grep("count", names(X)))
# # there is not Nham_nt == 3 variants with sufficient read counts
# ggpairs(X,
#   columns = grep("fitness_nt[12]", names(X)),
#   aes(alpha = 0.5, color = (count_nt1 > 1 | is.na(count_nt1)) & (count_nt2 > 1.5 | is.na(count_nt2)))
# )
# # for nt_2 variants, this looks better than in the epPCR stability data;
# # >> keep another 150 variants with Nham_nt==2 wit greater 10^1.5 reads

### double aa variants
# aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
# aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

# X <- merge(aa2_nt2, aa2_nt3)

# ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
# ggpairs(X,
#   columns = grep("fitness", names(X)),
#   aes(alpha = 0.5, color = ((count_nt2 > 2.5 | is.na(count_nt2))))
# )
# # nt==3 variants are just random
# # use same standards as epPCR_stabilityPCA
# # use only variants with greater 10^1.5 counts

all_variants <- all_variants[
  (WT == TRUE |
  (Nham_aa == 1 &  (Nham_nt == 1 | (Nham_nt == 2 & mean_count > 10^1.5))) |
  (Nham_aa == 2 & Nham_nt == 2 & mean_count > 10^1.5)) &
  STOP == FALSE & STOP_readthrough == FALSE,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2))
  ), aa_seq]

save(all_variants, file = "dg_models/SH3/data/01b-GRB2_epPCR_bindingPCA_aa012_thresholded.RData")





##################################################
### extract structural data from PDB structure ###
##################################################

### load distances between SH3 and ligand
SH3_GAB2_distances <- fread("dataset/GRB2-SH3/PDB_contactmap_2vwf_AB.txt")
SH3_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
SH3_GAB2_distances[, Pos := Pos1]
# add minimal side-chain heavy atom distance to ligand to singles data.table
# singles <- merge(singles, SH3_GAB2_distances[, .(HAmin_ligand = min(HAmin), scHAmin_ligand = min(scHAmin), GAB2_AA = WTAAPos2[which.min(HAmin)]), Pos1], by.x = "Pos", by.y = "Pos1")
structural_properties <- SH3_GAB2_distances[,
  .(HAmin_ligand = min(HAmin),
    scHAmin_ligand = min(scHAmin),
    GAB2_AA = WTAAPos2[which.min(HAmin)]
  ),
  Pos]


### load RSA values extracted using freeSASA
RSA_chainA <- fread("dataset/GRB2-SH3/2vwf_A.rsa", skip = 8, nrows = 56)
names(RSA_chainA) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel",
  "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel",
  "RSA_pol_abs", "RSA_pos_rel")
RSA_chainAB <- fread("dataset/GRB2-SH3/2vwf_AB.rsa", skip = 8, nrows = 56 + 14)
names(RSA_chainAB) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel",
  "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel",
  "RSA_pol_abs", "RSA_pos_rel")
# add to singles data.table
structural_properties <- merge(structural_properties,
  RSA_chainA[chain == "A", .(Pos, RSA_unbound = RSA_all_rel)], by = "Pos")
structural_properties <- merge(structural_properties,
  RSA_chainAB[chain == "A", .(Pos, RSA_bound = RSA_all_rel)], by = "Pos")
# >>> all residues that loose RSA when bound to ligand are on interface with ligand
# ggplot(structural_properties,aes(RSA_unbound,RSA_bound,color=HAmin_ligand < 5)) + geom_point(size=3)

## split data into core, surface and binding
# ggplot(structural_properties,aes(RSA_unbound,stability)) + geom_point() +geom_smooth()
# >>> set core threshold as RSA <= 25
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
  file = "dg_models/SH3/data/structural_properties.txt",
  row.names = F,
  quote = F
)

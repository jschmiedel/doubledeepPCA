### import and process GRB2 data
# singles from epPCR both assays and NM stability

## source functions
filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)

# create directory structure
dir.create("processed_data/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/preprocessing/", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
theme_set(theme_classic(base_size = 9))


### load the different data files
## wildtype
wildtype_epPCR_stabilityPCA <- 
  fread("dataset/GRB2/01a-GRB2_epPCR_stabilityPCA/fitness_wildtype.txt")
wildtype_epPCR_bindingPCA <- 
  fread("dataset/GRB2/01b-GRB2_epPCR_bindingPCA/fitness_wildtype.txt")
wildtype_NM_stabilityPCA <- 
  fread("dataset/GRB2/01c-GRB2_NM2_stabilityPCA/fitness_wildtype.txt")
wt_AA_seq <- wildtype_NM_stabilityPCA$aa_seq
wt_AA_seq_split <- strsplit(wt_AA_seq, "")[[1]]
# merge
wildtype <- merge(
  merge(
    wildtype_epPCR_stabilityPCA[, .(
      Pos = NA, WT, STOP, STOP_readthrough, Nham_nt, Nham_aa, Nmut_codons,
      mean_count_sEP = mean_count, 
      s_fitness_EP = 0,
      s_sigma_EP = sqrt(1 / (sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))
    )],
    wildtype_epPCR_bindingPCA[, .(
      Pos = NA, WT, STOP, STOP_readthrough, Nham_nt, Nham_aa, Nmut_codons,
      mean_count_bEP = mean_count, 
      total_reads_input = count_e1_s0 + count_e2_s0 + count_e3_s0, #for Otwinowski method
      total_reads_output = count_e1_s1 + count_e2_s1 + count_e3_s1, #for Otwinowski method
      b_fitness_EP = 0,
      b_sigma_EP = sqrt(1 / (sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))
    )]
  ),
  wildtype_NM_stabilityPCA[, .(
    Pos = NA, WT, STOP, STOP_readthrough, Nham_nt, Nham_aa, Nmut_codons,
    mean_count_sNM = mean_count, 
    s_fitness_NM = 0,
    s_sigma_NM = sqrt(1 / (sigma1_uncorr^-2 + sigma2_uncorr^-2 + sigma3_uncorr^-3))
  )]
)

##############################################
### intermediate data files before merging ###
##############################################
######################
#### NM stability ####
######################
load("dataset/GRB2/01c-GRB2_NM2_stabilityPCA/01c-GRB2_NM2_stabilityPCA_fitness_intermediate.RData")
# normalize by generations
all_variants[, g1 := 5.01]
all_variants[, g2 := 5.1]
all_variants[, g3 := 5.04]
all_variants[, fitness1 := log2(exp(1)) * fitness1_uncorr / g1]
all_variants[, fitness2 := log2(exp(1)) * fitness2_uncorr / g2]
all_variants[, fitness3 := log2(exp(1)) * fitness3_uncorr / g3]
all_variants[, sigma1 := log2(exp(1)) * sigma1_uncorr / g1]
all_variants[, sigma2 := log2(exp(1)) * sigma2_uncorr / g2]
all_variants[, sigma3 := log2(exp(1)) * sigma3_uncorr / g3]
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
all_variants[, mean_count := rowMeans(.SD), , .SDcols = grep("count_e[123]_s0", names(all_variants))]

# single AA variants
nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1 & !STOP,
  .(aa_seq, count_nt1 = log10(mean_count), fitness_nt1 = fitness)]
nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1 & !STOP,
  .(aa_seq, count_nt2 = log10(mean_count), fitness_nt2 = fitness)]
nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1 & !STOP,
  .(aa_seq, count_nt3 = log10(mean_count), fitness_nt3 = fitness)]

X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
theme_set(theme_bw())
ggpairs(X, columns = grep("count", names(X)))
ggpairs(X,
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = (count_nt1 > 3 | is.na(count_nt1)) & (count_nt2 > 2 | is.na(count_nt2)) & (count_nt3 > 2 | is.na(count_nt3)))
)
ggsave("results/preprocessing/GRB2_Nhamnt_count_fitness_dependencies_stabilityNM.pdf", width = 6, height = 6)
# if Nham_nt==1, variants below 10^3 read counts are sequencing errors
# if Nham_nt > 1, variants below 10^2 read counts are sequencing errors
# also there's a proportionality between readcounts of the same AA variants that could be taken into account; this is probably from sequencing errors as well


## merge based on hamming-distance specific read thresholds
all_data_rest <- all_variants[
  Nham_aa == 1 & Nmut_codons == 1 & ((Nham_nt == 1 & mean_count > 10^3) | (Nham_nt > 1 & mean_count > 10^2)),
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count), Nsyn = .N,
    Nham_aa = unique(Nham_aa), WT = unique(WT), STOP = unique(STOP), STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]
all_data <- all_variants[Nham_aa == 1 & Nmut_codons == 1 & mean_count > 10, .(
  fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
  sigma = sqrt(1 / sum(sigma^-2)),
  mean_count = sum(mean_count), Nsyn = .N,
  Nham_aa = unique(Nham_aa), WT = unique(WT), STOP = unique(STOP), STOP_readthrough = unique(STOP_readthrough)
), aa_seq]
X <- merge(
  all_data_rest[, .(aa_seq, fitness_rest = fitness, Nsyn_rest = Nsyn)],
  all_data[, .(aa_seq, fitness, Nsyn)]
)
ggplot(X[Nsyn > 1], aes(fitness, fitness_rest)) +
  geom_point()
# not a huge difference in total, but for ~30 variants
X[abs(fitness - fitness_rest) > 0.05]

## make this the singles dataset
singles_NM_stabilityPCA <- all_data_rest


### double aa variants
aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[1, 1])), fitness_nt2 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[1, 1])), fitness_nt3 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt4 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 4 & !STOP, .(count_nt4 = log10(unlist(.SD[1, 1])), fitness_nt4 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt5 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 5 & !STOP, .(count_nt5 = log10(unlist(.SD[1, 1])), fitness_nt5 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt6 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 6 & !STOP, .(count_nt6 = log10(unlist(.SD[1, 1])), fitness_nt6 = unlist(.SD[1, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

X <- merge(merge(merge(merge(aa2_nt2, aa2_nt3, all = T),
  aa2_nt4,
  all = T
),
aa2_nt5,
all = T
),
aa2_nt6,
all = T
)

theme_set(theme_bw())
ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
ggpairs(Y, columns = grep("count", names(X)))
ggpairs(X[sample(x = .N, size = 10^5)],
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = ((count_nt2 > 2 | is.na(count_nt2)) &
    (count_nt3 > 1.5 | is.na(count_nt3)) &
    (count_nt4 > 1.5 | is.na(count_nt4)) &
    (count_nt5 > 1.5 | is.na(count_nt5)) &
    (count_nt6 > 1.5 | is.na(count_nt6))))
)
ggpairs(X[sample(x = .N, size = 10^4)],
  columns = grep("fitness_nt[23]", names(X)),
  aes(alpha = 0.5, color = (count_nt2 > 1 | is.na(count_nt2)) &
    (count_nt3 > 1 | is.na(count_nt3)))
)
ggpairs(X[sample(x = .N, size = 10^4)],
  columns = grep("fitness_nt[34]", names(X)),
  aes(alpha = 0.5, color = (count_nt4 > 1.5 | is.na(count_nt4)) &
    (count_nt3 > 1.5 | is.na(count_nt3)))
)
ggsave("results/preprocessing/GRB2_doubles_Nhamnt_count_fitness_dependencies_stabilityNM.pdf", width = 6, height = 6)

X2 <- merge(aa2_nt2, aa2_nt3)
dt <- data.table(i = rep(seq(0, 3, 0.2), each = 16), j = seq(0, 3, 0.2), c = 0, n = 0)
for (idx in seq(0, 3, 0.2)) {
  for (jdx in seq(0, 3, 0.2)) {
    dt[i == idx & j == jdx, c := X2[count_nt2 > (idx) & count_nt3 > (jdx), cor(fitness_nt2, fitness_nt3)]]
    dt[i == idx & j == jdx, n := X2[count_nt2 > (idx) & count_nt3 > (jdx), .N]]
  }
}
ggplot(dt[n > 100], aes(i, j, fill = (c))) + geom_raster() + scale_fill_gradient2(midpoint = 0.85, mid = "grey")
ggplot(dt[n > 100], aes(i, j, fill = log10(n))) + geom_raster()
# exclude variants below 10^2 reads for aa2_nt2
# exclude viarants below 10^1.5 reads for all other

doubles_NM_stabilityPCA <- all_variants[
  Nham_aa == 2 & Nmut_codons == 2 &
      ((Nham_nt == 2 & mean_count > 10^2) | 
        (Nham_nt > 2 & mean_count > 10^1.5)),
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count),
    Nsyn = .N,
    Nham_aa = unique(Nham_aa),
    WT = unique(WT),
    STOP = unique(STOP),
    STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]

########################
#### epPCR stability ###
########################
load("dataset/GRB2/01a-GRB2_epPCR_stabilityPCA/01a-GRB2_epPCR_stabilityPCA_fitness_intermediate.RData")
# normalize by generations
all_variants[, g1 := 5.19]
all_variants[, g2 := 5.07]
all_variants[, g3 := 5.1]
all_variants[, fitness1 := log2(exp(1)) * fitness1_uncorr / g1]
all_variants[, fitness2 := log2(exp(1)) * fitness2_uncorr / g2]
all_variants[, fitness3 := log2(exp(1)) * fitness3_uncorr / g3]
all_variants[, sigma1 := log2(exp(1)) * sigma1_uncorr / g1]
all_variants[, sigma2 := log2(exp(1)) * sigma2_uncorr / g2]
all_variants[, sigma3 := log2(exp(1)) * sigma3_uncorr / g3]
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
all_variants[, mean_count := rowMeans(.SD), , 
  .SDcols = grep("count_e[123]_s0", names(all_variants))]


### single aa variants
nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1, .(count_nt1 = log10(unlist(.SD[, 1])), fitness_nt1 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
ggpairs(X, columns = grep("count", names(X)))
# there are no nt_3 variants with sufficient read counts
ggpairs(X,
  columns = grep("fitness_nt[12]", names(X)),
  aes(alpha = 0.5, color = (count_nt1 > 1 | is.na(count_nt1)) & (count_nt2 > 1 | is.na(count_nt2)))
)
# for nt_2 variants, it seems that only very high thresholds that leave no data contain real variants not swamped by sequencing errors
# compare nt_1 vairants to NM data
X <- merge(nt1, singles_NM_stabilityPCA[, .(aa_seq, fitness, mean_count = log10(mean_count))])
ggpairs(X, columns = c("count_nt1", "mean_count"))
ggpairs(X,
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = (count_nt1 > 3 | is.na(count_nt1)) & (mean_count > 3.5 | is.na(mean_count)))
)
# correlation is good not matter what the input read count of the epPCR library
## >> use all epPCR data from Nham_nt == 1
singles_epPCR_stabilityPCA <- all_variants[
  Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count), 
    Nsyn = .N,
    Nham_aa = unique(Nham_aa),
    WT = unique(WT),
    STOP = unique(STOP),
    STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]

### double aa variants
aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

X <- merge(aa2_nt2, aa2_nt3)

ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
ggpairs(X,
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = ((count_nt2 > 2.5 | is.na(count_nt2))))
)
# nt==3 variants are just random
# compare to NM data
X <- merge(aa2_nt2, doubles_NM_stabilityPCA[, .(aa_seq, fitness, mean_count = log10(mean_count))])
ggpairs(X,
  columns = c("count_nt2", "mean_count"),
  aes(alpha = 0.5, color = count_nt2 > 1.5)
)
ggpairs(X[count_nt2 > 1.5],
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = count_nt2 > 2)
)

# use only variants with greater 10^1.5 counts
doubles_epPCR_stabilityPCA <- all_variants[
  Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & mean_count > 10^1.5,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count), 
    Nsyn = .N,
    Nham_aa = unique(Nham_aa), 
    WT = unique(WT), 
    STOP = unique(STOP), 
    STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]






#######################
#### epPCR binding ####
#######################
load("dataset/GRB2/01b-GRB2_epPCR_bindingPCA/01b-GRB2_epPCR_bindingPCA_fitness_intermediate.RData")
# normalize by generations
all_variants[, g1 := 5.21]
all_variants[, g2 := 5.26]
all_variants[, g3 := 5.29]
all_variants[, fitness1 := log2(exp(1)) * fitness1_uncorr / g1]
all_variants[, fitness2 := log2(exp(1)) * fitness2_uncorr / g2]
all_variants[, fitness3 := log2(exp(1)) * fitness3_uncorr / g3]
all_variants[, sigma1 := log2(exp(1)) * sigma1_uncorr / g1]
all_variants[, sigma2 := log2(exp(1)) * sigma2_uncorr / g2]
all_variants[, sigma3 := log2(exp(1)) * sigma3_uncorr / g3]
all_variants[, total_reads_input := count_e1_s0 + count_e2_s0 + count_e3_s0]
all_variants[, total_reads_output := count_e1_s1 + count_e2_s1 + count_e3_s1]
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
all_variants[, mean_count := rowMeans(.SD), , .SDcols = grep("count_e[123]_s0", names(all_variants))]


### single aa variants
nt1 <- all_variants[Nham_aa == 1 & Nham_nt == 1 & Nmut_codons == 1, .(count_nt1 = log10(unlist(.SD[, 1])), fitness_nt1 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
nt2 <- all_variants[Nham_aa == 1 & Nham_nt == 2 & Nmut_codons == 1, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
nt3 <- all_variants[Nham_aa == 1 & Nham_nt == 3 & Nmut_codons == 1, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

X <- merge(merge(nt1, nt2, all = T), nt3, all = T)
ggpairs(X, columns = grep("count", names(X)))
# there is not Nham_nt == 3 variants with sufficient read counts
ggpairs(X,
  columns = grep("fitness_nt[12]", names(X)),
  aes(alpha = 0.5, color = (count_nt1 > 1 | is.na(count_nt1)) & (count_nt2 > 1.5 | is.na(count_nt2)))
)
# for nt_2 variants, this looks better than in the epPCR stability data;
# >> keep another 150 variants with Nham_nt==2 wit greater 10^1.5 reads
singles_epPCR_bindingPCA <- all_variants[
  Nham_aa == 1 & Nmut_codons == 1 & (Nham_nt == 1 | (Nham_nt == 2 & mean_count > 10^1.5)),
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count),
    total_reads_input = sum(total_reads_input, na.rm = T),
    total_reads_output = sum(total_reads_output, na.rm = T),
    Nsyn = .N,
    Nham_aa = unique(Nham_aa),
    WT = unique(WT),
    STOP = unique(STOP),
    STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]

### double aa variants
aa2_nt2 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & !STOP, .(count_nt2 = log10(unlist(.SD[, 1])), fitness_nt2 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]
aa2_nt3 <- all_variants[Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 3 & !STOP, .(count_nt3 = log10(unlist(.SD[, 1])), fitness_nt3 = unlist(.SD[, 2])), aa_seq, .SDcols = c("mean_count", "fitness")]

X <- merge(aa2_nt2, aa2_nt3)

ggpairs(X, columns = grep("count", names(X)), lower = list(continuous = "density"))
ggpairs(X,
  columns = grep("fitness", names(X)),
  aes(alpha = 0.5, color = ((count_nt2 > 2.5 | is.na(count_nt2))))
)
# nt==3 variants are just random
# use same standards as epPCR_stabilityPCA
# use only variants with greater 10^1.5 counts
doubles_epPCR_bindingPCA <- all_variants[
  Nham_aa == 2 & Nmut_codons == 2 & Nham_nt == 2 & mean_count > 10^1.5,
  .(
    fitness = sum(fitness * sigma^-2) / sum(sigma^-2),
    sigma = sqrt(1 / sum(sigma^-2)),
    mean_count = sum(mean_count),
    total_reads_input = sum(total_reads_input, na.rm = T),
    total_reads_output = sum(total_reads_output, na.rm = T),
    Nsyn = .N,
    Nham_aa = unique(Nham_aa),
    WT = unique(WT),
    STOP = unique(STOP),
    STOP_readthrough = unique(STOP_readthrough)
  ), aa_seq
]


#####################################
####### merge all single data #######
#####################################
singles <- merge(
  merge(
    singles_NM_stabilityPCA[, .(
      aa_seq, 
      STOP, 
      STOP_readthrough,
      mean_count_sNM = mean_count,
      s_fitness_NM = fitness, 
      s_sigma_NM = sigma
    )],
    singles_epPCR_stabilityPCA[, .(
      aa_seq,
      STOP,
      STOP_readthrough,
      mean_count_sEP = mean_count,
      s_fitness_EP = fitness, 
      s_sigma_EP = sigma
    )],
    all = T
  ),
  singles_epPCR_bindingPCA[, .(
    aa_seq, 
    STOP, 
    STOP_readthrough,
    mean_count_b = mean_count,
    total_reads_input,
    total_reads_output,
    b_fitness = fitness, 
    b_sigma = sigma
  )],
  all = T
)

## redo Pos/WT_AA/Mut stuff
singles[, Pos := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split), aa_seq]
singles[, WT_AA := wt_AA_seq_split[Pos], Pos]
singles[, Mut := strsplit(aa_seq, "")[[1]][Pos], aa_seq]

### doubles
doubles <- merge(
  merge(
    doubles_NM_stabilityPCA[, .(
      aa_seq, 
      STOP, 
      STOP_readthrough,
      mean_count_sNM = mean_count,
      s_fitness_NM = fitness, 
      s_sigma_NM = sigma
    )],
    doubles_epPCR_stabilityPCA[, .(
      aa_seq, 
      STOP, 
      STOP_readthrough,
      mean_count_sEP = mean_count,
      s_fitness_EP = fitness, 
      s_sigma_EP = sigma
    )],
    all = T
  ),
  doubles_epPCR_bindingPCA[, .(
    aa_seq, 
    STOP, 
    STOP_readthrough,
    mean_count_b = mean_count,
    total_reads_input,
    total_reads_output,
    b_fitness = fitness, 
    b_sigma = sigma
  )],
  all = T
)

## redo Pos/WT_AA/Mut stuff
doubles[, Pos1 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[1], aa_seq]
doubles[, Pos2 := which(strsplit(aa_seq, "")[[1]] != wt_AA_seq_split)[2], aa_seq]
doubles[, WT_AA1 := wt_AA_seq_split[Pos1], Pos1]
doubles[, WT_AA2 := wt_AA_seq_split[Pos2], Pos2]
doubles[, Mut1 := strsplit(aa_seq, "")[[1]][Pos1], aa_seq]
doubles[, Mut2 := strsplit(aa_seq, "")[[1]][Pos2], aa_seq]


### variant statistics
singles[, .N, .(!is.na(s_fitness_NM), !is.na(s_fitness_EP), !is.na(b_fitness))]
doubles[, .N, .(!is.na(s_fitness_NM), !is.na(s_fitness_EP), !is.na(b_fitness))]
cutoff <- 25
singles[, .N, .(mean_count_sNM > cutoff, mean_count_sEP > cutoff, mean_count_b > cutoff)]
# 326 variants in all four datasets
# 679 variants only in NM dataset
# retaining nt2 variants from epPCR binding data adds another 100 usable variants

### plot input count distributions between the datasets
ggpairs(singles, aes(alpha = 0.1),
  columns = grep("mean_count", names(singles), value = T),
  lower = list(continuous = ggpairs_log10_points),
  diag = list(continuous = ggpairs_log10_diagonal),
  upper = list(continous = ggpairs_log10_cor)
)
ggsave("results/preprocessing/GRB2_singles_count_distribution_3samples.pdf", width = 8, height = 8)
# as expected datasets from epPCR or NM are nearly identical
# single nt variants have similar abundance in epPCR and NM but different biases
singles[!is.na(s_fitness_EP), .(.N, mean(mean_count_sEP))]
singles[!is.na(s_fitness_NM), .(.N, mean(mean_count_sNM))]


### fitness versus input counts
X <- melt(singles, measure.vars = grep("mean_count", names(singles), value = T))[, .(Pos, Mut, mean_count = value, variable)]
Y <- melt(singles, measure.vars = grep("fitness", names(singles), value = T))[, .(Pos, Mut, fitness = value, variable)]
X[, id := strsplit(as.character(variable), "_")[[1]][3], variable]
Y[, id := paste0(strsplit(as.character(variable), "_")[[1]][c(1, 3)], collapse = ""), variable]
ggplot(merge(X[, .(mean_count, id, Pos, Mut)], Y[, .(fitness, id, Pos, Mut)]), aes(mean_count, fitness)) +
  geom_point() +
  scale_x_log10() +
  geom_vline(xintercept = c(15, 50)) +
  facet_wrap(id ~ .)
ggsave("results/preprocessing/GRB2_mean_count_threshold.pdf", width = 8, height = 5)
## 15 mean_count seems to be a good cutoff for epPCR
## for NM there has to be a higher cutoff, due to the overabundance of wildtype supplying many seq-error variants in output samples and thus skewing fitness values in lowly abundant variants

# merge fitness and sigma from NM and epPCR
ggplot(singles, aes(s_fitness_EP, s_fitness_NM)) +
  geom_density2d() +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline()
ggsave("results/preprocessing/GRB2_stabilityPCA_fitness_comparision1.pdf", width = 8, height = 8)

ggplot(melt(singles, measure.vars = grep("s_fitness", names(singles), value = T)), aes(value, color = variable)) +
  geom_density()
ggsave("results/preprocessing/GRB2_stabilityPCA_fitness_comparision2.pdf", width = 5, height = 5)

##### normalize fitness by scaling and shifting
set.seed(1603)
minF <- function(p) {
  F_norm <- (F_data + matrix(p[(Nsamples + 1):(2 * Nsamples)],
                            nrow = nrow(F_data),
                            ncol = Nsamples,
                            byrow = T)) * 
              matrix(p[1:Nsamples], 
                    nrow = nrow(F_data),
                    ncol = Nsamples,
                    byrow = T)
  F_avg <- rowMeans(F_data + matrix(p[(Nsamples + 1):(2 * Nsamples)],
                                  nrow = nrow(F_data),
                                  ncol = Nsamples,
                                  byrow = T))
  diffF <- sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}

### stability data
F_data <- singles[mean_count_sNM > 100 & mean_count_sEP > 100, 
            as.matrix(.SD), .SDcols = grep("s_fitness_(NM|EP)$", names(singles))]
Nsamples <- 2
x <- nlm(minF, rep(c(1, 0), each = Nsamples))
print(x)
p <- x$estimate
fitness_norm_model <- data.table(t(p))
names(fitness_norm_model) <- c("scale_sNM", "scale_sEP", "shift_sNM", "shift_sEP")

# wild-type correction such that mean(wild-type) = 0
wt_corr <- wildtype[, rowMeans((.SD + unlist(fitness_norm_model[, .SD, , .SDcols = grep("shift_s", names(fitness_norm_model))])) *
  unlist(fitness_norm_model[, .SD, , .SDcols = grep("scale_s", names(fitness_norm_model))])), ,
.SDcols = grep("s_fitness_(NM|EP)$", names(wildtype))
]
## normalize fitness values
singles[, s_fitness_NM_norm := 
  (s_fitness_NM + fitness_norm_model$shift_sNM) * 
  fitness_norm_model$scale_sNM - wt_corr]
singles[, s_fitness_EP_norm := 
  (s_fitness_EP + fitness_norm_model$shift_sEP) * 
  fitness_norm_model$scale_sEP - wt_corr]
doubles[, s_fitness_NM_norm := 
  (s_fitness_NM + fitness_norm_model$shift_sNM) * 
  fitness_norm_model$scale_sNM - wt_corr]
doubles[, s_fitness_EP_norm := 
  (s_fitness_EP + fitness_norm_model$shift_sEP) * 
  fitness_norm_model$scale_sEP - wt_corr]
## normalize error values
singles[, s_sigma_NM_norm := s_sigma_NM * fitness_norm_model$scale_sNM]
singles[, s_sigma_EP_norm := s_sigma_EP * fitness_norm_model$scale_sEP]
doubles[, s_sigma_NM_norm := s_sigma_NM * fitness_norm_model$scale_sNM]
doubles[, s_sigma_EP_norm := s_sigma_EP * fitness_norm_model$scale_sEP]
ggplot(
  melt(singles[mean_count_sNM > 100 & mean_count_sEP > 100], measure.vars = grep("s_fitness", names(singles), value = T)),
  aes((value), color = variable)
) +
  geom_density() +
  geom_vline(xintercept = 0)
# this is caused by low count variants, so it's fine

# merge fitness and error
singles[, s_fitness := sum(c(s_fitness_NM_norm / s_sigma_NM_norm^2, 
                             s_fitness_EP_norm / s_sigma_EP_norm^2), na.rm = T) /
                      sum(c(s_sigma_NM_norm^-2, s_sigma_EP_norm^-2), na.rm = T),
        .(Pos, Mut)]
singles[, s_sigma := sqrt(1 / sum(c(s_sigma_NM_norm^-2, s_sigma_EP_norm^-2), na.rm = T)),
        .(Pos, Mut)]
doubles[, s_fitness := sum(c(s_fitness_NM_norm / s_sigma_NM_norm^2, 
                             s_fitness_EP_norm / s_sigma_EP_norm^2), na.rm = T) /
                      sum(c(s_sigma_NM_norm^-2, s_sigma_EP_norm^-2), na.rm = T),
        .(Pos1, Mut1, Pos2, Mut2)]
doubles[, s_sigma := sqrt(1 / sum(c(s_sigma_NM_norm^-2, s_sigma_EP_norm^-2), na.rm = T)),
        .(Pos1, Mut1, Pos2, Mut2)]

# plot comparision
ggpairs(singles, columns = grep("[bs]_fitness(|_NM|_EP)$", names(singles), value = T))
ggsave(filename = "results/preprocessing/GRB2_fitness_merge.pdf", width = 6, height = 6)


### write table of processed single mutants
write.table(wildtype, file = "processed_data/GRB2_wildtype_alldata.txt",
  quote = F, row.names = F)
write.table(singles, file = "processed_data/GRB2_singles_alldata.txt",
  quote = F, row.names = F)
### write table of processed doubles mutants
write.table(doubles, file = "processed_data/GRB2_doubles_alldata.txt",
  quote = F, row.names = F)



### reintegrate this into analysis script
# ## add test for affecting stability or binding phenotype
# singles[, diff_bs_fitness := b_fitness - s_fitness]
# singles[, diff_bs_sigma := sqrt(b_sigma^2 + s_sigma^2)]

# # add p-values using t-test for the 3 replicates
# ttest <- function(av, se, df = 2, mu = 1) {
#   tstat <- (av - mu) / se
#   # Using T-dist
#   pval <- 2 * pt(abs(tstat), df, lower = FALSE)
#   return(pval)
# }
# singles[, pval_bind := p.adjust(ttest(b_fitness, b_sigma, mu = 0), method = "fdr")]
# singles[, pval_stab := p.adjust(ttest(s_fitness, s_sigma, mu = 0), method = "fdr")]
# singles[, pval_diff := p.adjust(ttest(diff_bs_fitness, diff_bs_sigma, mu = 0), method = "fdr")]

# pval_sig_thrs <- 0.1


# # mutations that affect binding (growth phenotype)
# singles[, affects_binding := pval_bind < pval_sig_thrs & b_fitness < 0]
# # mutations that affect stability (both biochem and growth)
# singles[, affects_stability := pval_stab < pval_sig_thrs & s_fitness < 0]
# # mutations effects binding (biochemical phenotype)
# singles[, affects_binding_biochem := pval_bind < pval_sig_thrs & pval_diff < pval_sig_thrs & diff_bs_fitness < 0]

# singles[, .N, .(affects_binding, affects_stability, affects_binding_biochem)]

##################################
### check library inequalities ###
##################################

X <- melt(singles, id.vars = c("Pos", "Mut"), measure.vars = grep("mean_count", names(singles)))
ggplot(X, aes(x = Pos, group = Pos, value)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(variable ~ .)
ggsave("results/preprocessing/GRB2_library_inequalities_byPos.pdf", width = 10, height = 6)

Tm_NM <- fread("dataset/GRB2/GRB2_Sigma_NNKlibrary.txt")
Tm_NM[, Pos := as.integer(gsub("GRB2_", "", Name)), Name]
singles_Tm_NM <- merge(singles[!is.na(mean_count_sNM), .(Pos, mean_count_sNM)], Tm_NM[, .(Pos, Tm)])
singles_Tm_NM[, cor(mean_count_sNM, Tm)]
ggplot(singles_Tm_NM, aes(Tm, mean_count_sNM)) +
  geom_point() +
  scale_y_log10() +
  geom_smooth() +
  labs(title = paste0("R=", singles_Tm_NM[, round(cor(mean_count_sNM, Tm), 2)]))
ggsave("results/preprocessing/GRB2_library_inequalities_NMvsTm.pdf", width = 6, height = 6)
# weak dependence on Tm


###############################
### integrate PDB structure ###
###############################

### load distances between GRB2 and ligand
GRB2_GAB2_distances <- fread("dataset/GRB2/PDB_contactmap_2vwf_AB.txt")
GRB2_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
GRB2_GAB2_distances[, Pos := Pos1]
# add minimal side-chain heavy atom distance to ligand to singles data.table
# singles <- merge(singles, GRB2_GAB2_distances[, .(HAmin_ligand = min(HAmin), scHAmin_ligand = min(scHAmin), GAB2_AA = WTAAPos2[which.min(HAmin)]), Pos1], by.x = "Pos", by.y = "Pos1")
variant_structuralproperties <- GRB2_GAB2_distances[, 
  .(HAmin_ligand = min(HAmin), 
    scHAmin_ligand = min(scHAmin), 
    GAB2_AA = WTAAPos2[which.min(HAmin)]
  ), 
  Pos]


### load RSA values
RSA_chainA <- fread("dataset/GRB2/2vwf_A.rsa", skip = 8, nrows = 56)
names(RSA_chainA) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel", 
  "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel", 
  "RSA_pol_abs", "RSA_pos_rel")
RSA_chainAB <- fread("dataset/GRB2/2vwf_AB.rsa", skip = 8, nrows = 56 + 14)
names(RSA_chainAB) <- c("restype", "aa", "chain", "Pos", "RSA_all_abs", "RSA_all_rel", 
  "RSA_ts_abs", "RSA_ts_rel", "RSA_mc_abs", "RSA_mc_rel", "RSA_np_abs", "RSA_np_rel", 
  "RSA_pol_abs", "RSA_pos_rel")
# add to singles data.table
variant_structuralproperties <- merge(variant_structuralproperties, 
  RSA_chainA[chain == "A", .(Pos, RSA_unbound = RSA_all_rel)], by = "Pos")
variant_structuralproperties <- merge(variant_structuralproperties, 
  RSA_chainAB[chain == "A", .(Pos, RSA_bound = RSA_all_rel)], by = "Pos")
# >>> all residues that loose RSA when bound to ligand are on interface with ligand
# ggplot(variant_structuralproperties,aes(RSA_unbound,RSA_bound,color=HAmin_ligand < 5)) + geom_point(size=3)

## split data into core, surface and binding
# ggplot(variant_structuralproperties,aes(RSA_unbound,stability)) + geom_point() +geom_smooth()
# >>> set core threshodl as RSA <= 10
threshold_RSA <- 10
threshold_ligand_dist <- 5
variant_structuralproperties[RSA_unbound <= threshold_RSA, type := "core"]
variant_structuralproperties[RSA_unbound > threshold_RSA & HAmin_ligand < threshold_ligand_dist, type := "ligand_binding"]
variant_structuralproperties[RSA_unbound > threshold_RSA & HAmin_ligand > threshold_ligand_dist, type := "surface"]
variant_structuralproperties[, type := factor(type, levels = c("core", "surface", "ligand_binding"))]
variant_structuralproperties[, .(Npos = length(unique(Pos))), type]

write.table(
  x = variant_structuralproperties,
  file = "processed_data/GRB2_variant_structuralproperties.txt",
  row.names = F,
  quote = F
)

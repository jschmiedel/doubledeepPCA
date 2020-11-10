# compare 3 state and 4 state thermodynamic models for SH3 domain

## source functions
filelist <- list.files("functions/")
invisible(sapply(paste0("functions/", filelist), source, .GlobalEnv))

#list of packages to install/load
packages = c("tempura","data.table","Matrix","ggplot2","GGally","ggpubr","stringr", "gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))
col_purple = "#9161A8"
col_orange = "#F7941E"


dataset_folder = "dg_models/SH3"
model_name3 = "three_state"
model_name4 = "four_state"

dir.create(file.path(dataset_folder, "analysis"))

# load varlist
load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
# load parlist
load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
# load estimated models
load(file = file.path(dataset_folder, model_name3, "model_results.RData"))
res3 <- model_results
load(file = file.path(dataset_folder, model_name4, "model_results.RData"))
res4 <- model_results

if (file.exists(file.path(dataset_folder, "data/structural_properties.txt"))) {
    structural_properties <- fread(file.path(dataset_folder, "data/structural_properties.txt"))
}

res3[["variant_data"]][, f_ddg3 := as.numeric(varlist[["varxmut"]] %*% res3[["avg_model"]][grep("f_ddg", parameter), boot_mean])]
res3[["variant_data"]][, b_ddg3 := as.numeric(varlist[["varxmut"]] %*% res3[["avg_model"]][grep("b_ddg", parameter), boot_mean])]
res3[["variant_data"]][, b_ddg_sd := as.numeric(varlist[["varxmut"]] %*% res3[["avg_model"]][grep("b_ddg", parameter), boot_sd])]
res3[["variant_data"]][, b_ddg_fdr := stats::p.adjust(2*stats::pnorm(-abs(b_ddg3 / b_ddg_sd)),
      method = "fdr")]
res4[["variant_data"]][, fA_ddg4 := as.numeric(varlist[["varxmut"]] %*% res4[["avg_model"]][grep("fA_ddg", parameter), boot_mean])]
res4[["variant_data"]][, fB_ddg4 := as.numeric(varlist[["varxmut"]] %*% res4[["avg_model"]][grep("fB_ddg", parameter), boot_mean])]
res4[["variant_data"]][, b_ddg4 := as.numeric(varlist[["varxmut"]] %*% res4[["avg_model"]][grep("b_ddg", parameter), boot_mean])]

res3[["variant_data"]][!grepl("_", aa_subs) & grepl("[0-9]", aa_subs), Pos := as.integer(paste0(strsplit(aa_subs,"")[[1]][1:(nchar(aa_subs)-1)], collapse = "")), aa_subs]

vd <- res3[["variant_data"]]
vd <- merge(vd,
                structural_properties,
                by = "Pos")
vd <- merge(vd, res4[["variant_data"]][,.(aa_seq, fA_ddg4, fB_ddg4, b_ddg4)], by = "aa_seq")



x <- vd[, .(fA_ddg = mean(fA_ddg4),
  fB_ddg = mean(fB_ddg4),
  b_ddg4 = mean(b_ddg4),
  b_ddg3 = mean(b_ddg3)),
  .(Pos, WT_AA, type, RSA_unbound, scHAmin_ligand)]

p_fAfB <- ggplot(x, aes(fA_ddg, fB_ddg, fill = scHAmin_ligand, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(linetype = 2) +
  geom_label_repel(aes(label = Pos)) +
  scale_fill_gradient2(midpoint = 8,
    low = col_purple,
    high = col_orange,
    mid = "grey") +
  scale_shape_manual(values = c(21, 23, 24, 25)) +
  labs(x = "folding ddG state A", y = "folding ddG state B")
p_fAb <- ggplot(x, aes(fA_ddg, b_ddg4, fill = scHAmin_ligand, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_label_repel(aes(label = Pos)) +
  scale_fill_gradient2(midpoint = (8),
    low = col_purple,
    high = col_orange,
    mid = "grey") +
  scale_shape_manual(values = c(21, 23, 24, 25)) +
  labs(x = "folding ddG state A", y = "binding ddG")
p_fBb <- ggplot(x, aes(fB_ddg, b_ddg4, fill = scHAmin_ligand, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_label_repel(aes(label = Pos)) +
  scale_fill_gradient2(midpoint = (8),
    low = col_purple,
    high = col_orange,
    mid = "grey") +
  scale_shape_manual(values = c(21, 23, 24, 25)) +
  labs(x = "folding ddG state B", y = "binding ddG")
p_b3b4 <- ggplot(x, aes(b_ddg3, b_ddg4, fill = scHAmin_ligand, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(linetype = 2) +
  geom_label_repel(aes(label = Pos), size = 2) +
  scale_fill_gradient2(midpoint = 8,
    low = col_purple,
    high = col_orange,
    mid = "grey") +
  scale_shape_manual(values = c(21, 23, 24, 25)) +
  labs(x = "binding ddG 3state model", y = "binding ddG 4state model")
ggsave(grid.arrange(p_fAfB, p_fAb, p_fBb, p_b3b4, nrow = 2, ncol = 2),
  file = file.path(dataset_folder, "analysis", "dG_scatter_perposition_4state.pdf"),
  width = 13, height = 10)








p <- ggplot(vd[,.(f_ddg3=mean(f_ddg3), b_ddg3=mean(b_ddg3), b_ddg4=mean(b_ddg4)),.(Pos, type)], aes(b_ddg3, b_ddg4, color = f_ddg3)) +
    geom_point() +
    # geom_smooth(method = "lm") +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    # scale_color_brewer(palette = "Set1") +
    labs(x = "three state model: ddG binding",
      y = "four state model: ddG binding")
ggplot2::ggsave(p,
  file = file.path(dataset_folder, "analysis", "comparision_3vs4_state_ddG_binding.pdf"),
  width = 5,
  height = 4)

p <- ggplot(vd[b_ddg3 != 0], aes(fB_ddg4-fA_ddg4, b_ddg4, color = type)) +
    geom_point() +
    # geom_smooth(se = F) +
    geom_hline(yintercept = 0, linetype = 2)
ggplot2::ggsave(p,
  file = file.path(dataset_folder, "analysis", "comparision_3vs4_state_ddG_binding_scHAmin.pdf"),
  width = 5,
  height = 4)


### investigate differences between cis-distance of mutations and prediction error
res3$variant_data <- get_pairdistance(res3$variant_data)

res3$variant_data[, mpw_bin := findInterval(min_pairdist, seq(4,30,4))]

ggplot(res3$variant_data[!is.na(min_pairdist)],
    aes(x = min_pairdist, y = (f1_fitness - f1_pred))) +
  geom_boxplot(aes(group = mpw_bin)) +
  geom_hline(yintercept = 0, linetype = 2)

ggplot(res3$variant_data[!is.na(min_pairdist)],
    aes(x = min_pairdist, y = (f2_fitness - f2_pred))) +
  geom_boxplot(aes(group = mpw_bin)) +
  geom_hline(yintercept = 0, linetype = 2)


ggplot(res3$variant_data[!is.na(min_pairdist)],
  aes(x = min_pairdist, y = (b1_fitness - b1_pred))) +
  geom_boxplot(aes(group = mpw_bin)) +
  geom_hline(yintercept = 0, linetype = 2)

# no shift detected in prediction performance as function of pairwise distances


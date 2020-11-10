# analyse and compare SH3 3state dg model

## source functions
filelist <- list.files("functions/")
invisible(sapply(paste0("functions/", filelist), source, .GlobalEnv))

#list of packages to install/load
packages = c("tempura","data.table","Matrix","ggplot2","GGally","ggpubr","stringr", "gridExtra", "ggrepel")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))
col_purple = "#9161A8"
col_blue =  "#0066CC"
col_orange = "#F7941E"
col_red = "#EF4136"


dataset_folder = "dg_models/SH3"
model_name = "three_state"

dir.create(file.path(dataset_folder, "analysis"))
dir.create(file.path(dataset_folder, "pymol"))

# load varlist
load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
# load parlist
load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
# load estimated models
load(file = file.path(dataset_folder, model_name, "model_results.RData"))
# load structural properties
structural_properties <- fread(file.path(dataset_folder, "data/structural_properties.txt"))

######## analyse ddG values
vd <- model_results[["variant_data"]]
vd[, f_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("f_ddg", parameter), boot_mean])]
vd[, f_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("f_ddg", parameter), boot_sd])]
vd[, b_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), boot_mean])]
vd[, b_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), boot_sd])]

vd <- vd[(!is.na(f1_fitness) | !is.na(f2_fitness)) & !is.na(b1_fitness)]

vd[!grepl("_", aa_subs) & grepl("[0-9]", aa_subs), Pos := as.integer(paste0(strsplit(aa_subs,"")[[1]][1:(nchar(aa_subs)-1)], collapse = "")), aa_subs]
vd <- merge(vd,
  structural_properties,
  by = "Pos")

vd[, Mut := gsub("[0-9]", "", aa_subs), aa_subs]
aa_details <- data.table(order = strsplit("RHKDESTNQCUGPAILMFWYV","")[[1]],
    type = c(rep("pos", 3),
             rep("neg", 2),
             rep("polar", 4),
             rep("special", 4),
             rep("hydro", 8)))

vd[, WT_AA := factor(WT_AA, levels = aa_details$order)]
vd[, Mut := factor(Mut, levels = aa_details$order)]
vd[, WT_AA_type := aa_details[order %in% WT_AA, type], WT_AA]
vd[, Mut_type := aa_details[order %in% Mut, type], Mut]
vd[, WT_Mut_type_x := paste0(WT_AA_type, "_", Mut_type), .(WT_AA_type, Mut_type)]

vd_perposition <- vd[, .(f_ddg = mean(f_ddg),
  b_ddg = mean(b_ddg)),
  .(Pos, WT_AA, type, RSA_unbound, scHAmin_ligand)]

### calculate state probabilities for mutants
rt = 1.99e-3 * 310.15
av <- model_results$avg_model
vd[,
  prob_f2 := exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg)/rt) /
      (1 + exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg)/rt))]
vd[,
  prob_f := exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg)/rt) /
      (1 + exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg)/rt) +
    exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg +
      av[parameter == "b_dgwt", boot_mean] + b_ddg)/rt))]
vd[,
  prob_fb := exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg +
      av[parameter == "b_dgwt", boot_mean] + b_ddg)/rt) /
    (1 + exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg)/rt) +
    exp(-(av[parameter == "f_dgwt", boot_mean] + f_ddg +
      av[parameter == "b_dgwt", boot_mean] + b_ddg)/rt))]

p <- ggpairs(vd,
  columns = c("f1_fitness", "b1_fitness",
    "f_ddg", "b_ddg",
    "prob_f2", "prob_f", "prob_fb"),
  ggplot2::aes(color = type, alpha = 0.3))

ggsave(p,
  file = file.path(dataset_folder, "analysis", "state_probabilities.pdf"))


### folding versus binding ddGs per position

p <- ggplot(vd_perposition, aes(f_ddg, b_ddg, fill = scHAmin_ligand, shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_label_repel(aes(label = Pos)) +
  scale_fill_gradient2(midpoint = (8),
    low = col_purple,
    high = col_orange,
    mid = "grey") +
  scale_shape_manual(values = c(21, 23, 24, 25)) +
  labs(x = "folding ddG", y = "binding ddG")
ggsave(p,
  file = file.path(dataset_folder, "analysis", "dG_scatter_perposition.pdf"))



######## linear model predictions
##### folding
attributes <- c("RSA_unbound",  "scHAmin_ligand", "WT_AA_type", "Mut_type", "WT_Mut_type_x")
idx <- t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1))
idx <- idx[,order(colSums(idx))]
idx <- idx[, idx[5,] == 0 | colSums(idx[3:4,]) == 0] #only use either WT_AA type and/or Mut type or their interaction
rownames(idx) <- attributes
lm_f_results <- copy(vd)
lm_f <- list()
adjR2 <- c()
for (i in 2:ncol(idx)) {
  formula0 <- paste0("f_ddg ~ ", paste0(attributes[which(idx[ ,i]==1)], collapse = " + "))
  lm_f[[i]] <- lm(formula0, data = lm_f_results)
  lm_f_results[, paste0("model_", paste0(idx[ ,i], collapse = "")) := lm_f[[i]]$fitted.values]
  adjR2[i] <- summary(lm_f[[i]])$adj.r.squared
}
rbind(idx,
  adjR2,
  predictors = colSums(idx),
  idx = 1:ncol(idx))[,order(colSums(idx),adjR2)]
# predictive value RSA ~ WT_Mut_interaction > WT_AA > Mut >> lig_dist = 0

attr_idx <- c(5, 1, 2)
model_idx <- c(6, 13, 19)
X <- data.table(x = "model",
  R2_adj = c(adjR2[model_idx[1]], diff(adjR2[model_idx])),
  attribute = factor(attributes[attr_idx], levels = attributes[attr_idx]))

p <- ggplot(X, aes(x = x, y = R2_adj, fill = attribute)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "adjusted R2 of folding ddGs")
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_fddg_R2contr.pdf"),
  width = 3,
  height = 4)


# best tradeoff between predictive power and complexity is from model with just RSA
# and WT_Mut_Interaction, R2_adj ~ 0.35

p <- ggplot(melt(lm_f_results, id.vars = c("f_ddg", "type", "Pos")),
    aes(value, f_ddg)) +
  geom_point(aes(color = type)) +
  geom_abline() +
  facet_wrap(variable ~ .) +
  stat_cor(aes(label = ..rr.label..)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'predicted', y = 'observed')
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_fddg.pdf"))

# investigate model using RSA and WT_Mut_interaction
summary(lm_f[[13]])
p <- ggplot(lm_f_results,
    aes(model_10001, f_ddg)) +
  geom_point(aes(color = type), shape = 1) +
  geom_abline() +
  geom_smooth() +
  geom_point(inherit.aes = F,
    data = lm_f_results[,.(x = mean(model_10001), y = mean(f_ddg)),. (Pos, type)],
    aes(x, y, color = type), size = 3, shape = 17) +
  scale_color_brewer(palette = "Set1") +
  stat_cor(aes(label = ..rr.label..))
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_fddg_RSA_WTMutX.pdf"))



##### binding
attributes <- c("RSA_unbound", "scHAmin_ligand", "WT_AA_type", "Mut_type", "WT_Mut_type_x", "f_ddg")
idx <- t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1))
idx <- idx[,order(colSums(idx))]
rownames(idx) <- attributes
idx <- idx[, idx[5,] == 0 | colSums(idx[3:4,]) == 0] #only use either WT_AA type and/or Mut type or their interaction

lm_f_results <- copy(vd[b_ddg != 0])
lm_f <- list()
adjR2 <- c()
for (i in 2:ncol(idx)) {
  formula0 <- paste0("b_ddg ~ ", paste0(attributes[which(idx[ ,i]==1)], collapse = " + "))
  lm_f[[i]] <- lm(formula0, data = lm_f_results)
  lm_f_results[, paste0("model_", paste0(idx[ ,i], collapse = "")) := lm_f[[i]]$fitted.values]
  adjR2[i] <- summary(lm_f[[i]])$adj.r.squared
}
rbind(idx,
  adjR2,
  predictors = colSums(idx),
  idx = 1:ncol(idx))[,order(colSums(idx),adjR2)]
# predictive value individually: lig_dist >> WT_AA_type ~ WT_Mut_type_x > RSA >> Mut_type ~ f_ddg = 0

attr_idx <- c(2, 5, 1, 6)
model_idx <- c(3, 15, 25, 39)
X <- data.table(x = "model",
  R2_adj = c(adjR2[model_idx[1]], diff(adjR2[model_idx])),
  attribute = factor(attributes[attr_idx], levels = attributes[attr_idx]))

p <- ggplot(X, aes(x = x, y = R2_adj, fill = attribute)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "adjusted R2 of binding ddGs")
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_bddg_R2contr.pdf"),
  width = 3,
  height = 4)


# best tradeoff between predictive power and complexity is from model with just ligand distance
# and WT_Mut_Interaction, R2_adj ~ 0.29

p <- ggplot(melt(lm_f_results, measure.vars = grep("model", names(lm_f_results), value = T)),
    aes(value, b_ddg)) +
  geom_point(aes(color = type)) +
  geom_abline() +
  facet_wrap(variable ~ .) +
  stat_cor(aes(label = ..rr.label..)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = 'predicted', y = 'observed')
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_bddg.pdf"))

#model with RSA and WT_Mut_interaction
summary(lm_f[[15]])
p <- ggplot(lm_f_results,
    aes(model_010010, b_ddg)) +
  geom_point(aes(color = type), shape = 1) +
  geom_abline() +
  geom_smooth() +
  geom_point(inherit.aes = F,
    data = lm_f_results[,.(x = mean(model_010010), y = mean(b_ddg)),. (Pos, type)],
    aes(x, y, color = type), size = 3) +
  geom_label_repel(inherit.aes = F,
    data = lm_f_results[,.(x = mean(model_010010), y = mean(b_ddg)),. (Pos, type)],
    aes(x, y, label = Pos, fill = type)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  stat_cor(aes(label = ..rr.label..))
ggsave(p,
  file = file.path(dataset_folder, "analysis", "lm_predict_bddg_ligdist_WTMutX.pdf"))



################ pymol figures

#### visualize b_ddg
script_file = file.path(dataset_folder, "pymol", "001_plot_ddGs.txt")
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
pymol_script[length(pymol_script)+1] = paste0("png 001_SH3_GAB2_reference_0.png, dpi=600")
# write pymol script in a .txt
write(x = pymol_script,file = file.path(dataset_folder, "pymol", "001_plot_ref.txt"))


# color GRB2 by folding ddG values
pymol_script[length(pymol_script)+1] = "set label_position, (1.7,1.7,1.7)"
pymol_script[length(pymol_script)+1] = "label chain A and n. ca, resi"
for (i in vd_perposition$Pos) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=", vd_perposition[Pos == i, f_ddg])
}
pymol_script[length(pymol_script)+1] = "select CA, chain A and name CA"
pymol_script[length(pymol_script)+1] = "show spheres, CA"
pymol_script[length(pymol_script)+1] = "set sphere_scale, 0.75, CA"
pymol_script[length(pymol_script)+1] = "set_bond stick_radius, 0.14, chain A"
pymol_script[length(pymol_script)+1] = "show sticks, chain A"
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white_blue, chain A, minimum=-2, maximum=2'
pymol_script[length(pymol_script)+1] = "set ray_opaque_background, 0"
pymol_script[length(pymol_script)+1] = "set ray_shadow, 0"
pymol_script[length(pymol_script)+1] = "set ray_trace_fog, 0"
pymol_script[length(pymol_script)+1] = "set antialias, 1"
pymol_script[length(pymol_script)+1] = "bg_color white"
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 002_SH3_GAB2_fddg_spheres.png, dpi=600")
write(x = pymol_script,file = file.path(dataset_folder, "pymol", "002_plot_fddg.txt"))

# color GRB2 by binding ddGs
for (i in vd_perposition$Pos) {
  pymol_script[length(pymol_script)+1]  = paste0("alter 2vwf and chain A and resid ",i,", b=", vd_perposition[Pos == i, b_ddg])
}
pymol_script[length(pymol_script)+1] = 'spectrum b, red_white_blue, chain A, minimum=-2, maximum=2'
pymol_script[length(pymol_script)+1] = "ray 2400,2400"
pymol_script[length(pymol_script)+1] = paste0("png 003_GRB2_GAB2_bddg_spheres.png, dpi=600")
write(x = pymol_script,file = file.path(dataset_folder, "pymol", "003_plot_bddg.txt"))

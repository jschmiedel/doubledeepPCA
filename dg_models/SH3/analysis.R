# analyse and compare SH3 dg models

require(tempura)
require(data.table)
require(Matrix)
require(ggplot2)
ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))

dataset_folder = "doubledeepPCA/dg_models/SH3"
# dataset_folder = "/home/joern/Dropbox/Science/CRG/doubledeepPCA/dg_models/SH3/"
model_name3 = "three_state"
model_name4 = "four_state"

setwd("~/")
dir.create(file.path(dataset_folder, "analysis"))

# load varlist
load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
# load parlist
# load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
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
                structural_properties[, list(Pos, type)],
                by = "Pos")
vd <- merge(vd, res4[["variant_data"]][,.(aa_seq, fA_ddg4, fB_ddg4, b_ddg4)], by = "aa_seq")

p <- ggplot(vd, aes(b_ddg3, b_ddg4, color = type)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "three state model: ddG binding",
      y = "four state model: ddG binding")

ggplot2::ggsave(p,
  file = file.path(dataset_folder, "analysis", "comparision_3vs4_state_ddG_binding.pdf"),
  width = 5,
  height = 4)

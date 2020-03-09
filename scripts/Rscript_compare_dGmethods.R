# compare dG methods

require(data.table)
require(ggplot2)
require(GGally)
theme_set(theme_bw(base_size = 9))

# GRB2 data
load("processed_data/dG_method1_GRB2_150models_10iterations.Rdata")
extract_bestmodel <- function(models) {
  parameters <- t(sapply(X = 1:length(models), FUN = function(X) {
    models[[X]]$par
  }))
  objective <- sapply(X = 1:length(models), FUN = function(X) {
    models[[X]]$value
  })
  L <- (length(models[[1]]$par) - 6) / 2
  energies <- data.table(
    mut = names(models[[which.min(objective)]]$par[1:L]),
    s_dG = models[[which.min(objective)]]$par[1:L],
    b_dG = models[[which.min(objective)]]$par[L + 1:L]
  )
  global_par_idx <- grep("[sb]_", colnames(parameters))
  global_par <- data.table(t(parameters[which.min(objective), global_par_idx]))
  return(list(energies, global_par))
}
method1 <- extract_bestmodel(models)
load("processed_data/dG_method2_GRB2_150models_10iterations.Rdata")
method2 <- extract_bestmodel(models)
load("processed_data/dG_method3_GRB2_150models_10iterations.Rdata")
method3 <- extract_bestmodel(models)

# using all variants
load("processed_data/dG_method1_GRB2_allvars_150models_10iterations.Rdata")
method1_allvars <- extract_bestmodel(models)

# otwin data
Otwin_energies <- fread("processed_data/GRB2_dG_Otwin_3state_energies.csv")
Otwin_global_par <- fread("processed_data/GRB2_dG_Otwin_3state_globalpar.csv")
Otwin_energies[2:.N, mut := paste0(pos, aa), .(pos, aa)]

Otwin_allvars_energies <- fread("processed_data/GRB2_dG_allvars_Otwin_3state_energies.csv")
Otwin_allvars_global_par <- fread("processed_data/GRB2_dG_allvars_Otwin_3state_globalpar.csv")
Otwin_allvars_energies[2:.N, mut := paste0(pos, aa), .(pos, aa)]


# compare energies
all_data <- merge(
  merge(
    merge(
      merge(
        merge(method1[[1]][, .(mut,
          s_dG1 = s_dG,
          b_dG1 = b_dG
        )],
        method2[[1]][, .(mut,
          s_dG2 = s_dG,
          b_dG2 = b_dG
        )],
        all = T
        ),
        method3[[1]][, .(mut,
          s_dG3 = s_dG,
          b_dG3 = b_dG
        )],
        all = T
      ),
      Otwin_energies[, .(mut, pos,
        s_dGO = efi,
        b_dGO = ebi
      )],
      all = T
    ),
    Otwin_allvars_energies[, .(mut,
      s_dGO_allvars = efi,
      b_dGO_allvars = ebi
    )],
    all = T
  ),
  method1_allvars[[1]][, .(mut,
    s_dG1_allvars = s_dG,
    b_dG1_allvars = b_dG
  )],
  all = T
)


struct <- fread(paste0("processed_data/GRB2_variant_structuralproperties.txt"))
all_data[, type := struct[Pos == pos, type], pos]

p = ggpairs(all_data,
  columns = grep("s_dG", names(all_data), value = T),
  # aes(color = type, alpha = 0.1),
  lower = list(continuous = "density")
)
ggsave(plot = p,
  filename = "results/dG/compare_methods_stability_energies.pdf")
p = ggpairs(all_data,
  columns = grep("b_dG", names(all_data), value = T),
  lower = list(continuous = "density")
)
ggsave(plot = p,
  filename = "results/dG/compare_methods_binding_energies.pdf")


# compare global parameters
all_pars = rbind(method1[[2]][, cbind(.SD, method = "1")],
  method1_allvars[[2]][, cbind(.SD, method = "1av")],
  method2[[2]][, cbind(.SD, method = "2")],
  method3[[2]][, cbind(.SD, method = "3")],
  Otwin_global_par[, cbind(.SD, method = "O")],
  Otwin_allvars_global_par[, cbind(.SD, method = "Oav")], fill = TRUE)

p = ggpairs(all_pars,
  column = grep("^[sb]", names(all_pars), value = T),
  aes(color = method)
)
ggsave(plot = p,
  filename = "results/dG/compare_methods_globalpar.pdf")
 

###### GB1 data Otwinowski
GB1_energies_org <- fread("Otwinoski2018_ProteinGthermo-master/Otwinoski2018_energies.csv")
GB1_energies_full <- fread("Otwinoski2018_ProteinGthermo-master/Ot_GB1_full_3state_energies.csv")
GB1_energies_ds1 <- fread("Otwinoski2018_ProteinGthermo-master/Ot_GB1_ds1_3state_energies.csv")
GB1_energies_ds1boot <- fread("Otwinoski2018_ProteinGthermo-master/Ot_GB1_ds1_boot_3state_energies.csv")

# load measurements from Nishal
GB1_Nishtal <- fread("dataset/dG/Nishtal2019_GB1_dG_stability.txt")
GB1_Nishtal[, pos := as.integer(paste0(strsplit(Description, "")[[1]][2:3], collapse = "")), Description]
GB1_Nishtal[, aa := strsplit(Description, "")[[1]][4], Description]
GB1_Nishtal[`Assay/Protocol` == "ddG(mAvg)_mean"][order(pos)]

X <- reshape(GB1_Nishtal[, .(pos, aa, Data, `Assay/Protocol`)], idvar = c("pos", "aa"), timevar = "Assay/Protocol", direction = "wide")
names(X) <- gsub("Data.", "", names(X))
names(X) <- gsub(" ", "", names(X))

Y <- X[`ddG(mAvg)_mean` > -4]
Y[, cor(.SD, use = "na.or.complete"), .SDcols = names(X)[!grepl("^(pos|aa|SD)", names(X))]]
Y[, cor(.SD, use = "pairwise.complete.obs"), .SDcols = names(X)[grepl("^ddG", names(X))]]

GB1_energies <- merge(
  merge(
    merge(
      merge(GB1_energies_full[, .(
        pos = pos + 1,
        aa,
        efi_full = efi,
        ebi_full = ebi
      )],
      GB1_energies_org[, .(pos,
        aa,
        efi_org = efi,
        ebi_org = ebi
      )],
      all = T
      ),
      GB1_energies_ds1[, .(
        pos = pos + 1,
        aa,
        efi_ds1 = efi,
        ebi_ds1 = ebi
      )],
      all = T
    ),
    GB1_energies_ds1boot[, .(
      pos = pos + 1,
      aa,
      efi_ds1boot = efi,
      ebi_ds1boot = ebi
    )],
    all = T
  ),
  Y[, .(pos,
    aa,
    dG_folding = `ddG(mAvg)_mean`
  )],
  all = T
)

ggpairs(GB1_energies,
  columns = grep("^[de]", names(GB1_energies), value = T),
  lower = list(continuous = "density")
)
GB1_energies[, cor(.SD, use = "pairwise.complete.obs"), .SDcols = grep("^[de]", names(GB1_energies))]

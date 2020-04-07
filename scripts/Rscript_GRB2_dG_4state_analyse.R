require(data.table)
require(ggplot2)
require(GGally)
require(gridExtra)
require(ggrepel)
theme_set(theme_bw(base_size = 9))

setwd("doubledeepPCA/")

filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)



dataset_name = "GRB2_dG_dataset"
method_vec <- c(2, 2)
predict_binding0_vec <- c(FALSE, TRUE)

for (i in 1) {
  method <- method_vec[i]
  predict_binding0 <- predict_binding0_vec[i]


  initialmodels <- fread(paste0("processed_data/dG/", 
                                  dataset_name, "_4state_method", method, 
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_initialmodels.txt"))

  boot_par_models <- fread(paste0("processed_data/dG/", 
                                  dataset_name, "_4state_method", method, 
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_boot_startpar.txt"))

  best_model <- boot_par_models[which.min(objective)]
  global_par <- best_model[, .SD, .SDcols = grep("^[sb]",names(best_model))]
  
  s1_ddg <- data.table(variant = sapply(X = grep("s1_ddg", names(best_model)),
                                       FUN = function(X) {strsplit(names(best_model)[X], "_")[[1]][1]}),
                      s1_ddg = best_model[, unlist(.SD), .SDcols = grep("s1_ddg", names(best_model))])
  s2_ddg <- data.table(variant = sapply(X = grep("s2_ddg", names(best_model)),
                                       FUN = function(X) {strsplit(names(best_model)[X], "_")[[1]][1]}),
                      s2_ddg = best_model[, unlist(.SD), .SDcols = grep("s2_ddg", names(best_model))])
  b_ddg <- data.table(variant = sapply(X = grep("b_ddg", names(best_model)),
                                       FUN = function(X) {strsplit(names(best_model)[X], "_")[[1]][1]}),
                      b_ddg = best_model[, unlist(.SD), .SDcols = grep("b_ddg", names(best_model))])

  table_dg <- merge(s1_ddg, s2_ddg, all = T)
  table_dg <- merge(table_dg, b_ddg, all = T)

  all_data <- fread(paste0("processed_data/", dataset_name, ".txt"))
  all_data[, s1_ddg1 := table_dg[variant == id1,s1_ddg],id1]
  all_data[, s2_ddg1 := table_dg[variant == id1,s2_ddg],id1]
  all_data[, b_ddg1 := table_dg[variant == id1,b_ddg],id1]
  all_data[, s1_ddg2 := table_dg[variant == id2,s1_ddg],id2]
  all_data[, s2_ddg2 := table_dg[variant == id2,s2_ddg],id2]
  all_data[, b_ddg2 := table_dg[variant == id2,b_ddg],id2]

  # calculate stabilityPCA fitness from dG values
  all_data[!is.na(s_fitness), s_fitness_pred := function_dG_4state_folding_dg2fitness(
    s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
    s2_ddg = sum(c(s2_ddg1, s2_ddg2), na.rm = T),
    s1_dgwt = global_par$s1_dgwt,
    s2_dgwt = global_par$s2_dgwt,
    s_f0 = ifelse("s_f0" %in% names(global_par), 
                  global_par$s_f0, 0), 
    s_fwt = ifelse("s_fwt" %in% names(global_par), 
                   global_par$s_fwt, 1)
  ), 
  .(s1_ddg1, s1_ddg2, s2_ddg1, s2_ddg2)
  ]

  # calculate bindingPCA fitness with estimated b_ddG values
  all_data[!is.na(b_fitness), 
           b_fitness_pred := function_dG_4state_binding_dg2fitness(
             b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = sum(c(s2_ddg1, s2_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(b_ddg1, b_ddg2, s1_ddg1, s1_ddg2, s2_ddg1, s2_ddg2)
           ]

  # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[!is.na(b_fitness), 
           b_fitness_pred_bdg0 := function_dG_4state_binding_dg2fitness(
             b_ddg = 0,
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = sum(c(s2_ddg1, s2_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(s1_ddg1, s1_ddg2, s2_ddg1, s2_ddg2)
           ]

  # calculate bindingPCA fitness if s_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[!is.na(b_fitness), 
           b_fitness_pred_sdg0 := function_dG_4state_binding_dg2fitness(
             b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
             s1_ddg = 0,
             s2_ddg = 0,
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(b_ddg1, b_ddg2)
           ]
  all_data[!is.na(b_fitness), 
           b_fitness_pred_s1eqs2 := function_dG_4state_binding_dg2fitness(
             b_ddg = 0,
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(s1_ddg1, s1_ddg2)
           ]
  all_data[!is.na(b_fitness), 
           s_fitness_pred_s1eqs2 := function_dG_4state_folding_dg2fitness(
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           ),
           .(s1_ddg1, s1_ddg2)
           ]

  all_data[!is.na(b_fitness), 
           b_fitness_pred_onlys2 := function_dG_4state_binding_dg2fitness(
             b_ddg = 0,
             s1_ddg = 0,
             s2_ddg = sum(c(s2_ddg1, s2_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(s2_ddg1, s2_ddg2)
           ]

  all_data[!is.na(b_fitness), 
           s_fitness_pred_onlys2 := function_dG_4state_folding_dg2fitness(
             s1_ddg = 0,
             s2_ddg = sum(c(s2_ddg1, s2_ddg2), na.rm = T),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           ),
           .(s2_ddg1, s2_ddg2)
           ]
  all_data[!is.na(b_fitness), 
           b_fitness_pred_onlys1 := function_dG_4state_binding_dg2fitness(
             b_ddg = 0,
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = 0,
             b_dgwt = global_par$b_dgwt,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(s1_ddg1, s1_ddg2)
           ]

  all_data[!is.na(b_fitness), 
           s_fitness_pred_onlys1 := function_dG_4state_folding_dg2fitness(
             s1_ddg = sum(c(s1_ddg1, s1_ddg2), na.rm = T),
             s2_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           ),
           .(s1_ddg1, s1_ddg2)
           ]
  ## relationship between parameters and fitness
  # load structural property file to color residues by type
  variant_structuralproperties = fread(paste0("processed_data/GRB2_variant_structuralproperties.txt"))
  # all_data[,res_type1 := variant_structuralproperties[Pos == Pos1,type],Pos1]
  all_data[,RSA_unbound := variant_structuralproperties[Pos == Pos1,RSA_unbound],Pos1]
  all_data[,HAmin_ligand := variant_structuralproperties[Pos == Pos1,HAmin_ligand],Pos1]
  all_data[RSA_unbound <= 20 & HAmin_ligand > 5,res_type1 := "core"]
  all_data[RSA_unbound <= 20 & HAmin_ligand <= 5,res_type1 := "core_bind"]
  all_data[RSA_unbound > 20 & HAmin_ligand > 5,res_type1 := "surf"]
  all_data[RSA_unbound > 20 & HAmin_ligand <= 5,res_type1 := "surf_bind"]
  all_data[Nmut == 1, .N, res_type1]



  ##############################
  ### plot model performance ###
  ##############################
  all_models <- rbind(initialmodels[, 
        cbind(.SD, type = "initial"), 
        .SDcols = grep("^[osb]",names(initialmodels))],
      boot_par_models[, 
        cbind(.SD, type = "bootpar"), 
        .SDcols = grep("^[osb]",names(boot_par_models))])
  all_models[, objective := log10(objective)]
  P <- ggpairs(all_models[objective < 5.43],
    columns = grep("^[osb]",names(all_models)))
  ggsave(P,
      file = paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""), 
        "_modelperformance.pdf"),
      width = 10,
      height = 10)


  ###########################
  ### state probabilities ###
  ###########################
  #stabilityPCA
  Pf1 = (exp(-global_par$s1_dgwt/rt)) / 
      (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt))
  Pf2 = (exp(-global_par$s2_dgwt/rt)) / 
      (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt))
  # P1_P2 = (exp(-global_par$s2_dgwt/rt)+exp(-global_par$s1_dgwt/rt)) / 
      # (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt))

  # Pu = 1 - Pf1 - Pf2
  # print(paste0("StabilityPCA: folded state = ", round(Pf*100,1), "%"))
  global_par_dt = data.table(experiment = "stabilityPCA",
                              state = c("f1", "f2"),
                              P = c(Pf1, Pf2))
  # print(paste0("StabilityPCA: folded state 1 = ", round(P1*100,1), "%"))
  # print(paste0("StabilityPCA: folded state 2 = ", round(P2*100,1), "%"))
  # print(paste0("StabilityPCA: both folded states = ", round(P1_P2*100,1), "%"))
  #bindingPCA
  Pf1 = exp(-global_par$s1_dgwt/rt) / 
      (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt) + exp(-(global_par$s2_dgwt + global_par$b_dgwt)/rt))
  Pf2 = exp(-global_par$s2_dgwt/rt) / 
      (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt) + exp(-(global_par$s2_dgwt + global_par$b_dgwt)/rt))
  Pf2b = exp(-(global_par$s2_dgwt + global_par$b_dgwt)/rt) / 
      (1 + exp(-global_par$s1_dgwt/rt) + exp(-global_par$s2_dgwt/rt) + exp(-(global_par$s2_dgwt + global_par$b_dgwt)/rt))
  # Pu = 1 - Pf1 - Pf2 - Pf2b
  global_par_dt = rbind(global_par_dt,
    data.table(experiment = "bindingPCA",
                              state = c("f1", "f2", "f2b"),
                              P = c(Pf1, Pf2, Pf2b)))

  global_par_dt[, experiment := factor(experiment, levels = c("stabilityPCA", "bindingPCA"))]
  global_par_dt[, state := factor(state, levels = c("f1", "f2b", "f2"))]
  P <- ggplot(global_par_dt, aes(x = experiment, group = experiment, y = P, fill = state)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1")
  ggsave(P, 
    file = paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""),
         "_stateprobs.pdf"),
    width = 4, height = 4)
  # print(paste0("bindingPCA: folded state 1 = ", round(P1b*100,1), "%"))
  # print(paste0("bindingPCA: folded state 2 = ", round(P2b*100,1), "%"))
  # print(paste0("bindingPCA: state 2 bound to GAB2 = ", round(Pb*100,1), "%"))
  # print(paste0("bindingPCA: all folded states = ", round((P1b + P2b + Pb)*100,1), "%"))

  ######################
  ### R2 on test-set ###
  ######################

  R2_stab = cor(all_data[test_set==TRUE,.(s_fitness,s_fitness_pred)])[2,1]^2
  R2_stab_w = cov.wt(all_data[test_set==TRUE,.(s_fitness,s_fitness_pred)],wt = all_data[test_set==TRUE,1/s_sigma^2],cor = T)$cor[2,1]^2
  x = matrix(all_data[test_set == TRUE,s_fitness + rnorm(100*.N,mean = 0, sd = s_sigma)],
             nrow = all_data[test_set == TRUE,.N], ncol = 100)
  tmp = cor(x)
  R2_stab_nontech = mean(tmp[lower.tri(tmp)])^2
  R2_stab/R2_stab_nontech
  tmp = cov.wt(x,wt = all_data[test_set==TRUE,1/s_sigma^2],cor = T)
  R2_stab_nontech_w = mean(tmp$cor[lower.tri(tmp$cor)])^2
  R2_stab_w/R2_stab_nontech_w

  ps <- ggplot(all_data[test_set == TRUE],aes(s_fitness_pred,s_fitness)) +
    geom_point() +
    # geom_density2d(aes(color = (id1 %in% idlist | id2 %in% idlist)), show.legend = F) +
    geom_abline(color = "red") +
    geom_smooth() +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    geom_label(inherit.aes = F,
               data = data.table(x = 0,
                                 y = c(1, 0.9),
                                 label = c(paste0("R2 = ",round(R2_stab,2), " (", round(R2_stab/R2_stab_nontech*100,1), "% of total)"),
                                           paste0("R2_weighted = ",round(R2_stab_w,2), " (", round(R2_stab_w/R2_stab_nontech_w*100,1), "% of total)"))),
               aes(x,y,label=label),
               hjust = 0, size = 2.5) +
    labs(y = "fitness stabilityPCA",
         x = "predicted fitness")

  R2_bind = cor(all_data[test_set==TRUE,.(b_fitness,b_fitness_pred)])[2,1]^2
  R2_bind_w = cov.wt(all_data[test_set==TRUE,.(b_fitness,b_fitness_pred)],wt = all_data[test_set==TRUE,1/b_sigma^2],cor = T)$cor[2,1]^2
  x = matrix(all_data[test_set == TRUE,b_fitness + rnorm(100*.N,mean = 0, sd = b_sigma)],
             nrow = all_data[test_set == TRUE,.N], ncol = 100)
  tmp = cor(x)
  R2_bind_nontech = mean(tmp[lower.tri(tmp)])^2
  R2_bind/R2_bind_nontech
  tmp = cov.wt(x,wt = all_data[test_set==TRUE,1/b_sigma^2],cor = T)
  R2_bind_nontech_w = mean(tmp$cor[lower.tri(tmp$cor)])^2
  R2_bind_w/R2_bind_nontech_w
  pb <- ggplot(all_data[test_set == TRUE],aes(b_fitness_pred,b_fitness)) +
    geom_point() +
    # geom_density2d(aes(color = (id1 %in% idlist | id2 %in% idlist)), show.legend = F) +
    geom_abline(color = "red") +
    geom_smooth() +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    geom_label(inherit.aes = F,
               data = data.table(x = 0,
                                 y = c(1, 0.9),
                                 label = c(paste0("R2 = ",round(R2_bind,2), " (", round(R2_bind/R2_bind_nontech*100,1), "% of total)"),
                                           paste0("R2_weighted = ",round(R2_bind_w,2), " (", round(R2_bind_w/R2_bind_nontech_w*100,1), "% of total)"))),
               aes(x,y,label=label),
               hjust = 0, size = 2.5) +
    labs(y = "fitness bindingPCA",
         x = "predicted fitness")
  P <- grid.arrange(ps, pb, nrow=1)
  ggsave(P, 
    file = paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""),
         "_R2.pdf"),
    width = 8, height = 4)


  ###############################################################
  ### dG value relation with each other and to fitness values ###
  ###############################################################


  X = copy(all_data[Nmut == 1])

  X[,s_ddg1_diff := s1_ddg1 - s2_ddg1]
  #higher in s2_ddg1 than s1_ddg1
  X[s1_ddg1 - s2_ddg1 < -1.5 & s1_ddg1 > 1 & Nid > 5,.(Pos1, id1, Nid, 
    s_fitness, b_fitness, 
    s_fitness_pred, b_fitness_pred, 
    s1_ddg1, s2_ddg1, b_ddg1,
    RSA_unbound, HAmin_ligand, res_type1)][order(Pos1)]
  # positions 11, 2*19, 25, 27, 32, 36, 42, 3*45, 47 with varying b_ddg1
  #positions 3, 36, 45, 52 of which 3, 36 and 52 cluster on structure

  unique(X[Nid >= 3,.(
    RSA_unbound, HAmin_ligand, res_type1,
    Ns2_high = sum(s1_ddg1 - s2_ddg1 < -0.2),
    Ns1_high = sum(s1_ddg1 - s2_ddg1 > -0.2)),Pos1])[order(Ns2_high)]


  Y = X[Nid >= 3,.(RSA = unique(RSA_unbound),
    HAmin_ligand = unique(HAmin_ligand), 
    res_type1 = unique(res_type1),
    Ns2_high = sum(s1_ddg1 - s2_ddg1 < -0.2),
    Ns1_high = sum(s1_ddg1 - s2_ddg1 > -0.2),
    mean_diff = mean(s1_ddg1 - s2_ddg1 + 0.2),
    mean_diff_abs = mean(abs(s1_ddg1 - s2_ddg1 + 0.2))),Pos1][order(mean_diff_abs)]

  ggplot(Y, aes(mean_diff_abs, mean_diff, color = HAmin_ligand < 11)) +
    geom_point(size = 3) +
    # geom_label_repel(aes(label = Pos1)) +
    # scale_color_gradient2(midpoint = 12, mid = "grey")
    labs(x = "avg. abs(s1_ddg1 - s2_ddg1)", y = "s2_ddg higher <--- [[avg. (s1_ddg - s2_ddg)]] ---> s1_ddg higher")
  ggsave(paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""),
        "_diff_s1_s2.pdf"))

  #higher in s1_ddg1 than s2_ddg1
  X[s1_ddg1 - s2_ddg1 > 1.5 & s2_ddg1 > 1 & Nid > 5,.(Pos1, id1, Nid, 
    s_fitness, b_fitness, 
    s_fitness_pred, b_fitness_pred, 
    s1_ddg1, s2_ddg1, b_ddg1,
    RSA_unbound, HAmin_ligand, res_type1)][order(b_ddg1)]
  # positions 19, 22, 29, 2*31, 2*40, 47 with low b_ddg
  # positions 9, 25, 27 & 52 with high b_ddg
  
  X[between(s1_ddg1 - s2_ddg1,0, 2.5) & Nid > 15,.(Pos1, id1, Nid, 
    s_fitness, b_fitness, 
    s_fitness_pred, b_fitness_pred, 
    s1_ddg1, s2_ddg1, b_ddg1,
    RSA_unbound, HAmin_ligand, res_type1)][order(Pos1)]

  P <- ggpairs(X[Nmut == 1 & !is.na(b_fitness) & (Nid + Nbid) > 50 & (Nid + Nsid) > 50 & s_fitness > 0.4 & s1_ddg1 - s2_ddg1 < 1.5],
          columns = c(grep("ddg1", names(X)),
                      grep("fitness_pred$", names(X)),
                      grep("fitness$", names(X))),
          aes(color = res_type1, alpha = 0.25),
          # aes(color = s_fitness > 0.4 , alpha = 0.5),
          # aes(color = s1_ddg1 - s2_ddg1 < -1.5 & s1_ddg1 > 1, alpha = 0.25),
          # aes(color = between(s1_ddg1 - s2_ddg1,-1.5, -0.5), alpha = 0.25),
          # aes(color = between(s1_ddg1 - s2_ddg1,0, 1.5), alpha = 0.25),
          # aes(color = b_ddg1 > 2, alpha = 0.25),
          # aes(color = s1_ddg1 - s2_ddg1 > 1.5, alpha = 0.25),
          upper = list(continuous = wrap("cor", size = 3)))
  ggsave(P, 
    file = paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""),
         "_dG_vs_fitness.pdf"),
    width = 12, height = 12)
  
# ### follow some mutants
#   idlist = X[Nmut == 1 & 
#     !is.na(b_fitness) & (
#     Nid + Nbid) > 50 & 
#     (Nid + Nsid) > 50 & 
#     s_fitness > 0.4 & 
#     s1_ddg1 - s2_ddg1 > 1.5, 
#     id1]

#   all_data[test_set == TRUE & (id1 %in% idlist | id2 %in% idlist)]
# # R2_stab_w = 
# cov.wt(all_data[test_set == F & !is.na(s_fitness) & Nmut == 2, .(s_fitness, s_fitness_pred)],
#   wt = all_data[test_set == F & !is.na(s_fitness) & Nmut == 2, 1 / s_sigma^2],
#   cor = T)$cor[2,1]^2
# cov.wt(all_data[test_set == F & !is.na(s_fitness) & Nmut == 2 & (id1 %in% idlist | id2 %in% idlist), .(s_fitness, s_fitness_pred)],
#   wt = all_data[test_set == F & !is.na(s_fitness) & Nmut == 2 & (id1 %in% idlist | id2 %in% idlist), 1 / s_sigma^2],
#   cor = T)$cor[2,1]^2
# # R2_bind_w = 
# cov.wt(all_data[test_set == F & !is.na(b_fitness) & Nmut == 2, .(b_fitness, b_fitness_pred)],
#   wt = all_data[test_set == F & !is.na(b_fitness) & Nmut == 2, 1 / b_sigma^2], 
#   cor = T)$cor[2,1]^2
# cov.wt(all_data[test_set == F & !is.na(b_fitness) & Nmut == 2 & (id1 %in% idlist | id2 %in% idlist), .(b_fitness, b_fitness_pred)],
#   wt = all_data[test_set == F & !is.na(b_fitness) & Nmut == 2 & (id1 %in% idlist | id2 %in% idlist), 1 / b_sigma^2], 
#   cor = T)$cor[2,1]^2


  p1 <- ggplot(X[Nmut == 1 & !is.na(b_fitness)],
          aes(s1_ddg1, s2_ddg1)) +
        # geom_point(aes(color = res_type1)) +
        geom_density2d() +
        geom_abline(intercept = c(-1,1)) +
        geom_smooth(method = "lm") +
        geom_abline(color = "red", intercept = 0.25)
  p2 <- ggplot(X[Nmut == 1 & !is.na(b_fitness)],
          aes(s1_ddg1, b_ddg1, color = res_type1)) +
        geom_point() +
        geom_abline(color = "red")
  p3 <- ggplot(X[Nmut == 1 & !is.na(b_fitness)],
          aes(s2_ddg1, b_ddg1, color = res_type1)) +
        geom_point() +
        geom_abline(color = "red")
  P <- grid.arrange(p1, p2, p3, nrow = 2)
  ggsave(P, 
    file = paste0("results/dG/", dataset_name, "_4state_method", method, 
        ifelse(predict_binding0 == TRUE, "_b0", ""),
        "_dG_vs_dist.pdf"),
    width = 8, height = 8)

  ###############################################################
  ### dG value relation with each other and to fitness values ###
  ###############################################################
  
  ds1_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(s1_ddg1 + global_par$s1_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = c(0.85, 0.75)) +
    labs(x = "dG folding1", color = "")

  ds2_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(s2_ddg1 + global_par$s2_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$s2_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = c(0.85, 0.75)) +
    labs(x = "dG folding2", color = "")

  ds1_fs <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
      aes(s1_ddg1 + global_par$s1_dgwt, s_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], 
      aes(s1_ddg1 + global_par$s1_dgwt, s_fitness, color = res_type1)) +
    geom_line(data = data.table(x = global_par$s1_dgwt + seq(-.75,5,0.05),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = seq(-.75,5,0.05),
             s2_ddg = seq(-.75,5,0.05) * sqrt(2),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y),color = "black") +
    geom_line(data = data.table(x = global_par$s1_dgwt + seq(-.75,5,0.05),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = seq(-.75,5,0.05),
             s2_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y),color = "black", linetype = 3) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$s_f0,global_par$s_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) +
    labs(x = "dG folding1", y = "fitness stabilityPCA", color = "") + theme(legend.position = "none")

  ds2_fs <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
      aes(s2_ddg1 + global_par$s2_dgwt, s_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], 
      aes(s2_ddg1 + global_par$s2_dgwt, s_fitness, color = res_type1)) +
    geom_line(data = data.table(x = global_par$s2_dgwt + seq(-1,5,0.05),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = seq(-1,5,0.05) / sqrt(2),
             s2_ddg = seq(-1,5,0.05),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y),color = "black") +
    geom_line(data = data.table(x = global_par$s2_dgwt + seq(-1,5,0.05),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = 0,
             s2_ddg = seq(-1,5,0.05),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y),color = "black", linetype = 3) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$s_f0,global_par$s_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s2_dgwt,linetype = 2) +
    labs(x = "dG folding2", y = "fitness stabilityPCA", color = "") + theme(legend.position = "none")


    s1VSs2 <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s1_ddg1 + global_par$s1_dgwt, 
                                                            s2_ddg1 + global_par$s2_dgwt)) +
    geom_density2d(color = "black") +
    # geom_point(aes(color = res_type1)) +
    # scale_color_brewer(palette = "Set1") +
    geom_point(aes(color = s_fitness)) +
    geom_abline(slope =  sqrt(2), intercept = -0.2) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_y_continuous(breaks = seq(-2,6,1)) +
    geom_abline(color = "red", intercept = global_par[, s2_dgwt - s1_dgwt]) +
    geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$s2_dgwt,linetype = 2) +
    labs(x = "dG folding1", y = "dG folding2", color = "") + theme(legend.position = "none")


  P <- grid.arrange(ds1_dist, ds1_fs, ds2_fs,
                    s1VSs2, ds2_dist,
                   nrow = 2)
  ggsave(P, file = paste0("results/dG/", dataset_name, "_4state_method", method, 
    ifelse(predict_binding0 == TRUE, "_b0", ""), 
    "_dG_sfitness.pdf"), 
         width = 11, height = 7)

  db_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    coord_flip() +
    labs(x = "dG binding") + theme(legend.position = "none")

  s1VSb <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s1_ddg1 + global_par$s1_dgwt, 
                                                            b_ddg1 + global_par$b_dgwt)) +
    geom_density2d(color = 'black') +
    geom_point(aes(color = res_type1)) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_y_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding1", y = "dG binding", color = "") + theme(legend.position = "none")

  s2VSb <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s2_ddg1 + global_par$s2_dgwt, 
                                                            b_ddg1 + global_par$b_dgwt)) +
    geom_density2d(color = 'black') +
    geom_point(aes(color = res_type1)) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_y_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s2_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding2", y = "dG binding", color = "") + theme(legend.position = "none")

  ds1_fb <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
        aes(s1_ddg1 + global_par$s1_dgwt, b_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], aes(s1_ddg1 + global_par$s1_dgwt, 
                                                             b_fitness, color = res_type1)) +
    geom_line(data = data.table(x = global_par$s1_dgwt + seq(-.75,5,0.05),
        y = function_dG_4state_binding_dg2fitness(
             s1_ddg = seq(-.75,5,0.05),
             s2_ddg = seq(-.75,5,0.05) * sqrt(2),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_line(data = data.table(x = global_par$s1_dgwt + seq(-.75,5,0.05),
        y = function_dG_4state_binding_dg2fitness(
             s1_ddg = seq(-.75,5,0.05),
             s2_ddg = 0,
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black", linetype = 3) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.15)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) +
    labs(x = "dG folding1", y = "fitness bindingPCA", color = "") + theme(legend.position = "none")

  ds2_fb <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
        aes(s2_ddg1 + global_par$s2_dgwt, b_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], aes(s2_ddg1 + global_par$s2_dgwt, 
                                                             b_fitness, color = res_type1)) +
    geom_line(data = data.table(x = global_par$s2_dgwt + seq(-1,5,0.05),
        y = function_dG_4state_binding_dg2fitness(
             s1_ddg = seq(-1,5,0.05) / sqrt(2),
             s2_ddg = seq(-1,5,0.05),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_line(data = data.table(x = global_par$s2_dgwt + seq(-1,5,0.05),
        y = function_dG_4state_binding_dg2fitness(
             s1_ddg = 0,
             s2_ddg = seq(-1,5,0.05),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black", linetype = 3) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.15)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s2_dgwt,linetype = 2) +
    labs(x = "dG folding2", y = "fitness bindingPCA", color = "") + theme(legend.position = "none")

  db_fb <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
        aes(b_ddg1  + global_par$b_dgwt, b_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], aes(b_ddg1 + global_par$b_dgwt, 
                                                             b_fitness, color = res_type1)) +
    geom_line(data = all_data[Nmut==1 & !is.na(b_fitness)],
              aes(b_ddg1 + global_par$b_dgwt, 
                  b_fitness_pred_sdg0),color = "black") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0,1.15)) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG binding", y = "fitness bindingPCA", color = "") + theme(legend.position = "none")

  P <- grid.arrange(db_dist, s1VSb, s2VSb,
                  ds1_fb,ds2_fb,db_fb,
                   nrow = 2)
  ggsave(P, 
    file = paste0("results/dG/", dataset_name, "_4state_method", method, 
    ifelse(predict_binding0 == TRUE, "_b0", ""),
    "_dG_bfitness.pdf"), 
         width = 11, height = 7)

  

  

  

  
  
  # P <- grid.arrange(ds_dist, ds_fs, ds_fb, sVSb, db_dist, db_fb,
  #                  nrow = 2)
  # ggsave(P, file = paste0("results/dG/GRB2_dG_method", method, 
  #   ifelse(predict_binding0 == TRUE, "_b0", ""), "_dG_fitness.pdf"), 
  #        width = 11, height = 7)


  ##############################
  ### s_fitness vs b_fitness ###
  ##############################
  
  change_bdgwt = 0
  p <- ggplot() +
    geom_point(data = all_data[Nmut == 1 & !is.na(b_fitness)], aes(
      b_fitness, 
      s_fitness,
      color = res_type1,
      # size = ifelse(abs(b_ddg1) < 0.5, 0, abs(b_ddg1)),
      shape = factor(ifelse(abs(b_ddg1) < 0.5, 0, sign(b_ddg1)))
    )) +
    geom_line(data = data.table(
        x = function_dG_4state_binding_dg2fitness(
             s1_ddg = seq(-1,5,0.05) / sqrt(2),
             s2_ddg = seq(-1,5,0.05),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt + change_bdgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = seq(-1,5,0.05) / sqrt(2),
             s2_ddg = seq(-1,5,0.05),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y), color = "black") +
    geom_line(data = data.table(
        x = function_dG_4state_binding_dg2fitness(
             s1_ddg = 0,
             s2_ddg = seq(-1,5,0.05),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt + change_bdgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = 0,
             s2_ddg = seq(-1,5,0.05),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y), color = "black", linetype = 2) +
    geom_line(data = data.table(
        x = function_dG_4state_binding_dg2fitness(
             s1_ddg = seq(-1,5,0.05),
             s2_ddg = seq(-1,5,0.05),
             b_ddg = 0,
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             b_dgwt = global_par$b_dgwt + change_bdgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
        y = function_dG_4state_folding_dg2fitness(
             s1_ddg = seq(-1,5,0.05),
             s2_ddg = seq(-1,5,0.05),
             s1_dgwt = global_par$s1_dgwt,
             s2_dgwt = global_par$s2_dgwt,
             s_f0 = global_par$s_f0,
             s_fwt = global_par$s_fwt
           )),
      aes(x,y), color = "black", linetype = 2) +
    labs(
      x = "fitness bindingPCA", 
      y = "fitness stabilityPCA",
      color = "residue type",
      shape = "significant dG binding"
    ) +
    scale_x_continuous(limits = c(0.15, 1.15)) +
    # scale_size_continuous(range = c(1, 3)) +
    scale_color_brewer(palette = "Set1")
  ggsave(plot = p, 
         file = paste0("results/dG/GRB2_dG_4state_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_fitness_globalrelationship.pdf"), 
         width = 6.5, 
         height = 4.5)



  ##### uncertainty of fitness values
  ### bootstrapped fitness values
  boot_fitness_models <- fread(paste0("processed_data/dG/", 
                                      dataset_name, "_4state_method", method, 
                                      ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                      "_boot_fitness.txt"))
  par_sd = rep(0, ncol(boot_fitness_models))
  for (i in 1:ncol(boot_fitness_models)) {
    par_sd[i] = boot_fitness_models[, sd(unlist(.SD)), .SDcols = names(boot_fitness_models)[i]]
  }
  boot_par_sd = data.table(parameter = names(boot_fitness_models),sd = par_sd)

  boot_par_sd[grep("ddg", parameter), variant := strsplit(parameter,"_")[[1]][1], parameter]

  global_par_long = melt(global_par)
  names(global_par_long) = c("parameter", "estimate")
  global_par_long = merge(global_par_long, boot_par_sd)

  table_dg <- merge(table_dg, boot_par_sd[grep("s1_ddg", parameter),.(variant, s1_ddg_sd = sd)], all = T)
  table_dg <- merge(table_dg, boot_par_sd[grep("s2_ddg", parameter),.(variant, s2_ddg_sd = sd)], all = T)
  table_dg <- merge(table_dg, boot_par_sd[grep("b_ddg", parameter),.(variant, b_ddg_sd = sd)], all = T)
  table_dg
  # how does sd relate to # doubles seen?
  all_data[!is.na(s_fitness),N_stability := .N + sum(all_data[!is.na(s_fitness),id2] == id1, na.rm = T), id1]
  all_data[!is.na(b_fitness),N_binding := .N + sum(all_data[!is.na(b_fitness),id2] == id1, na.rm = T), id1]
  table_dg = merge(table_dg,all_data[Nmut == 1,.(variant = id, N_stability, N_binding)], all.x = T)

  p1 = ggplot(table_dg[N_stability > 10 & N_binding > 5], aes(s1_ddg, s2_ddg)) +
    geom_density2d() +
    geom_point(aes(color = s1_ddg_sd)) +
    # geom_errorbarh(aes(xmin = s_ddg - 3*s_ddg_sd, xmax = s_ddg + 3*s_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding1", y = "dG folding2", color = "sd(s1_ddg)")

  p2 = ggplot(table_dg[N_stability > 10& N_binding > 5], aes(s1_ddg, s2_ddg)) +
    geom_density2d() +
    geom_point(aes(color = s2_ddg_sd)) +
    # geom_errorbarh(aes(xmin = s_ddg - 3*s_ddg_sd, xmax = s_ddg + 3*s_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding1", y = "dG folding2", color = "sd(s2_ddg)")

  p3 = ggplot(table_dg[N_stability > 10& N_binding > 5], aes(s1_ddg, b_ddg)) +
    geom_density2d() +
    geom_point(aes(color = b_ddg_sd)) +
    # geom_errorbar(aes(ymin = b_ddg - 3*b_ddg_sd, ymax = b_ddg + 3*b_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding1", y = "dG binding", color = "sd(b_ddg)")
  p4 = ggplot(table_dg[N_stability > 10& N_binding > 5], aes(s2_ddg, b_ddg)) +
    geom_density2d() +
    geom_point(aes(color = b_ddg_sd)) +
    # geom_errorbar(aes(ymin = b_ddg - 3*b_ddg_sd, ymax = b_ddg + 3*b_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding2", y = "dG binding", color = "sd(b_ddg)")

  p = grid.arrange(p1, p2, p3, p4, nrow = 2)
  ggsave(plot = p, 
         file = paste0("results/dG/GRB2_dG_4state_method", method, 
                ifelse(predict_binding0 == TRUE, "_b0", ""), 
                "_dGpar_bootsd.pdf"), 
         width = 8, 
         height = 6.5)


  ##########################
  ### profile likelihood ###
  ##########################

  pll_bounds <- fread(paste0("processed_data/dG/", 
                             dataset_name, "_4state_method", method,
                             ifelse(predict_binding0 == TRUE, "_b0", ""), 
                             "_parameter_pll_df", ncol(best_model) - 2, ".txt"))
  pll_bounds[grep("ddg", par), variant := strsplit(par,"_")[[1]][1], par]
  table_pll_dg <- merge(
    merge(pll_bounds[grep("s1_ddg", par),
                               .(variant, 
                                 s1_ddg = value, 
                                 s1_ddg_upper = upper, 
                                 s1_ddg_lower = lower)],
          pll_bounds[grep("s2_ddg", par),
                               .(variant, 
                                 s2_ddg = value, 
                                 s2_ddg_upper = upper, 
                                 s2_ddg_lower = lower)], all = T),
          pll_bounds[grep("b_ddg", par),
                               .(variant, 
                                 b_ddg = value, 
                                 b_ddg_upper = upper, 
                                 b_ddg_lower = lower)], all = T)
  # p1 = ggplot(table_pll_dg, aes(s_ddg, b_ddg, 
  #     color = ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0))) +
  #   geom_point() +
  #   # scale_color_gradient(low = "black", high = "red") +
  #   # geom_errorbarh(aes(xmin = s_ddg_lower, xmax = s_ddg_upper)) +
  #   labs(x = "dG folding", y = "dG binding", color = "dG folding log(90% conf)")

  # p2 = ggplot(table_pll_dg, aes(s_ddg, b_ddg, 
  #     color = ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0))) +
  #   geom_point() +
  #   # scale_color_gradient(low = "black", high = "red") +
  #   # geom_errorbar(aes(ymin = b_ddg_lower, ymax = b_ddg_upper)) +
  #   labs(x = "dG folding", y = "dG binding", color = "dG binding log(90% conf)")

  color_restype <- setNames(c("red2", "forestgreen", "blue2", "mediumpurple"), c("core", "surf", "core_bind", "surf_bind"))

  sig_sddg1 <- table_pll_dg[ifelse(s1_ddg > (s2_ddg / sqrt(2)), 
        s1_ddg_lower > (s2_ddg / sqrt(2)), s1_ddg_upper < (s2_ddg / sqrt(2))),variant]
  sig_sddg2 <- table_pll_dg[ifelse(s2_ddg > (s1_ddg * sqrt(2)), 
        s2_ddg_lower > (s1_ddg * sqrt(2)), s2_ddg_upper < (s1_ddg * sqrt(2))),variant]
  sig_bddg <- table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]

  scaling <- 2
  sig_sddg1 <- table_pll_dg[ifelse(s1_ddg > (s2_ddg / sqrt(2)), 
        (s1_ddg - (s1_ddg - s1_ddg_lower)/scaling) > (s2_ddg / sqrt(2)), 
        (s1_ddg + (s1_ddg_upper - s1_ddg)/scaling) < (s2_ddg / sqrt(2))),
    variant]
  sig_sddg2 <- table_pll_dg[ifelse(s2_ddg > (s1_ddg * sqrt(2)), 
        (s2_ddg - (s2_ddg - s2_ddg_lower)/scaling) > (s1_ddg * sqrt(2)), 
        (s2_ddg + (s2_ddg_upper - s2_ddg)/scaling) < (s1_ddg * sqrt(2))),
    variant]
  sig_bddg <- table_pll_dg[ifelse(b_ddg > 0, 
        (b_ddg - (b_ddg - b_ddg_lower)/scaling) > 0, 
        (b_ddg + (b_ddg_upper - b_ddg)/scaling) < 0),
    variant]

  p1 <- ggplot() +
    geom_density2d(data = table_pll_dg[!is.na(b_ddg)], aes(s1_ddg + global_par$s1_dgwt, s2_ddg + global_par$s2_dgwt), color = "black") +
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & !(id1 %in% sig_sddg1)],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & id1 %in% sig_sddg1],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt, color = res_type1)) +
    geom_abline(slope = sqrt(2), intercept = -0.2, linetype = 2) +
    scale_color_manual(values  = color_restype) +
    theme(legend.position = "none") +
    labs(x = "dG folding1", y = "dG folding2", title = "sig. dG folding1 != wt")

  p2 <- ggplot() +
    geom_density2d(data = table_pll_dg[!is.na(b_ddg)], aes(s1_ddg + global_par$s1_dgwt, s2_ddg + global_par$s2_dgwt), color = "black") +
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & !(id1 %in% sig_sddg2)],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & id1 %in% sig_sddg2],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt, color = res_type1)) +
    geom_abline(slope = sqrt(2), intercept = -0.2, linetype = 2) +
    scale_color_manual(values  = color_restype) +
    theme(legend.position = "none") +
    labs(x = "dG folding1", y = "dG folding2", title = "sig. dG folding2 != wt")

  p3 <- ggplot() +
    geom_density2d(data = table_pll_dg[!is.na(b_ddg)], aes(s1_ddg + global_par$s1_dgwt, s2_ddg + global_par$s2_dgwt), color = "black") +
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & (!(id1 %in% sig_sddg1) | !(id1 %in% sig_sddg2))],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & id1 %in% sig_sddg1 & id1 %in% sig_sddg2],
      aes(s1_ddg1 + global_par$s1_dgwt, s2_ddg1 + global_par$s2_dgwt, color = res_type1)) +
    geom_abline(slope = sqrt(2), intercept = -0.2, linetype = 2) +
    scale_color_manual(values  = color_restype) +
    theme(legend.position = "none") +
    labs(x = "dG folding1", y = "dG folding2", title = "sig. dG folding 1 & 2 != wt")

  p7 <- ggplot() +
    geom_density2d(data = table_pll_dg[!is.na(b_ddg)], 
      aes(s1_ddg + global_par$s1_dgwt,
      -log(exp(-(s1_ddg + global_par$s1_dgwt)) + exp(-(s2_ddg + global_par$s2_dgwt)))), 
      color = "black") +
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & (!(id1 %in% sig_sddg1) | !(id1 %in% sig_sddg2))],
      aes(s1_ddg1 + global_par$s1_dgwt, 
        -log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt)))),
      alpha = 0.2) + 
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & id1 %in% sig_sddg1 & id1 %in% sig_sddg2],
      aes(s1_ddg1 + global_par$s1_dgwt, 
        -log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt))),
         color = res_type1)) +
    geom_abline(slope = sqrt(2), intercept = -0.2, linetype = 2) +
    scale_color_manual(values  = color_restype) +
    theme(legend.position = "none") +
    labs(x = "dG folding 1", y = "dG folding 1&2", title = "sig. dG folding 1 & 2 != wt")

  p8 <- ggplot() +
    geom_density2d(data = table_pll_dg[!is.na(b_ddg)], 
      aes(s2_ddg + global_par$s2_dgwt,
        -log(exp(-(s1_ddg + global_par$s1_dgwt)) + exp(-(s2_ddg + global_par$s2_dgwt)))), 
    color = "black") +
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & (!(id1 %in% sig_sddg1) | !(id1 %in% sig_sddg2))],
      aes(s2_ddg1 + global_par$s2_dgwt, 
        -log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt)))), 
      alpha = 0.2) + 
    geom_point(data = all_data[!is.na(b_ddg1) & Nmut == 1 & id1 %in% sig_sddg1 & id1 %in% sig_sddg2],
      aes(s2_ddg1 + global_par$s2_dgwt,
        -log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt))), 
        color = res_type1)) +
    geom_abline(slope = sqrt(2), intercept = -0.2, linetype = 2) +
    scale_color_manual(values  = color_restype) +
    theme(legend.position = "none") +
    labs(x = "dG folding 2", y = "dG folding 1&2", title = "sig. dG folding 1 & 2 != wt")
  # p1 <- ggplot() +
  #   geom_density2d(data = table_pll_dg, aes(s1_ddg + global_par$s1_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
  #   geom_point(data = all_data[Nmut == 1 & !(id1 %in% table_pll_dg[ifelse(s1_ddg > 0, s1_ddg_lower > 0, s1_ddg_upper < 0),variant])],
  #     aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
  #   geom_point(data = all_data[Nmut == 1 & id1 %in% table_pll_dg[ifelse(s1_ddg > 0, s1_ddg_lower > 0, s1_ddg_upper < 0),variant]],
  #     aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
  #   scale_color_brewer(palette = "Set1") +
  #   geom_vline(xintercept = global_par$s1_dgwt,linetype = 2) + theme(legend.position = "none") +
  #   labs(x = "dG folding1", y = "dG binding", color = "sig. dG folding1 != wt")

  # p2 <- ggplot() +
  #   geom_density2d(data = table_pll_dg, aes(s2_ddg + global_par$s2_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
  #   geom_point(data = all_data[Nmut == 1 & !(id1 %in% table_pll_dg[ifelse(s2_ddg > 0, s2_ddg_lower > 0, s2_ddg_upper < 0),variant])],
  #     aes(s2_ddg1 + global_par$s2_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
  #   geom_point(data = all_data[Nmut == 1 & id1 %in% table_pll_dg[ifelse(s2_ddg > 0, s2_ddg_lower > 0, s2_ddg_upper < 0),variant]],
  #     aes(s2_ddg1 + global_par$s2_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
  #   scale_color_brewer(palette = "Set1") +
  #   geom_vline(xintercept = global_par$s2_dgwt,linetype = 2) + theme(legend.position = "none") +
  #   labs(x = "dG folding2", y = "dG binding", color = "sig. dG folding2 != wt")

  p4 <- ggplot() +
    geom_density2d(data = table_pll_dg, aes(s1_ddg + global_par$s1_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & !(id1 %in% sig_bddg)],
      aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & id1 %in% sig_bddg],
      aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_manual(values  = color_restype) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) + theme(legend.position = "none") +
    labs(x = "dG folding1", y = "dG binding", color = "sig. dG binding != wt")

  p5 <- ggplot() +
    geom_density2d(data = table_pll_dg, aes(s2_ddg + global_par$s2_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & !(id1 %in% sig_bddg)],
      aes(s2_ddg1 + global_par$s2_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & id1 %in% sig_bddg],
      aes(s2_ddg1 + global_par$s2_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_manual(values  = color_restype) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) + theme(legend.position = "none") +
    labs(x = "dG folding2", y = "dG binding", color = "sig. dG binding != wt")

  p6 <- ggplot() +
    geom_density2d(data = table_pll_dg, 
      aes(-log(exp(-(s1_ddg + global_par$s1_dgwt)) + exp(-(s2_ddg + global_par$s2_dgwt))), 
      b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & (!(id1 %in% sig_bddg) | !(id1 %in% sig_sddg1) | !(id1 %in% sig_sddg2))],
      aes(-log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt))), 
        b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & id1 %in% sig_bddg & id1 %in% sig_sddg1 & id1 %in% sig_sddg2],
      aes(-log(exp(-(s1_ddg1 + global_par$s1_dgwt)) + exp(-(s2_ddg1 + global_par$s2_dgwt))), 
        b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_manual(values  = color_restype) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) + theme(legend.position = "none") +
    labs(x = "dG folding", y = "dG binding", color = "sig. dG binding != wt")

  # p4 <- ggplot() +
  #   geom_density2d(data = table_pll_dg, aes(s1_ddg + global_par$s1_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
  #   geom_point(data = all_data[Nmut == 1 & 
  #       (!(id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]) |
  #       !(id1 %in% table_pll_dg[ifelse(s1_ddg > 0, s1_ddg_lower > 0, s1_ddg_upper < 0),variant]))],
  #     aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
  #   geom_point(data = all_data[Nmut == 1 & 
  #       id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]  &
  #       id1 %in% table_pll_dg[ifelse(s1_ddg > 0, s1_ddg_lower > 0, s1_ddg_upper < 0),variant]],
  #     aes(s1_ddg1 + global_par$s1_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
  #   scale_color_brewer(palette = "Set1") +
  #   geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
  #   geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
  #   # labs(x = "dG folding", y = "dG binding", color = "sig. dG binding != wt & dG folding != wt")
  #   labs(x = "dG folding1", y = "dG binding", color = "")

  p = grid.arrange(grobs = list(p1, p2, p3,
      p4, p5, p6,
      p7, p8), 
      nrow = 3)
      # widths = c(1.33,0.33,1),
      # layout_matrix = rbind(c(1,2,2),
                            # c(3,3,NA)))
  ggsave(plot = p, 
         file = paste0("results/dG/", dataset_name,"_4state_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_dGpar_pll", ncol(best_model)-2, ".pdf"), 
         width = 12.5, 
         height = 10.5)


  # all_data[Nmut == 1 & 
  #       id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]  &
  #       id1 %in% table_pll_dg[ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0),variant] & res_type1 == "core"]
  #global parameters
  print(pll_bounds[grep("^[sb]",par)])


}
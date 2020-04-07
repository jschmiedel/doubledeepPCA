require(data.table)
require(ggplot2)
require(GGally)
require(gridExtra)
theme_set(theme_bw(base_size = 9))

setwd("doubledeepPCA/")

filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)



dataset_name = "GRB2_dG_dataset"
method_vec <- c(2)
predict_binding0_vec <- FALSE

i = 1
for (i in 1:5) {

  method <- method_vec[i]
  predict_binding0 <- predict_binding0_vec[i]


  initialmodels <- boot_par_models <- fread(paste0("processed_data/dG/", 
                                  dataset_name, "_sdgwt_method", method, 
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_initialmodels.txt"))

  boot_par_models <- fread(paste0("processed_data/dg/", 
                                  dataset_name, "_sdgwt_method", method, 
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_boot_startpar.txt"))

  ## plot model performance
  all_models <- rbind(initialmodels[, 
        cbind(.SD, type = "initial"), 
        .SDcols = grep("^[osb]",names(initialmodels))],
      boot_par_models[, 
        cbind(.SD, type = "bootpar"), 
        .SDcols = grep("^[osb]",names(boot_par_models))])
  all_models[, objective := log10(objective)]
  P <- ggpairs(all_models,#[objective < 5.6],
    columns = grep("^[osb]",names(all_models)))
  ggsave(P,
    file = paste0("results/dG/GRB2_dG_sdgwt_method", method, 
      ifelse(predict_binding0 == TRUE, "_b0", ""), "_modelperformance.pdf"),
    width = 10,
    height = 10)

  best_model <- boot_par_models[which.min(objective)]
  best_model <- initialmodels[which.min(objective)]
  global_par <- best_model[, .SD, .SDcols = grep("^[sb]",names(best_model))]

  #stabilityPCA
  rt = 1.99e-3 * 310.15
  Pf = (exp(-global_par$ss_dgwt/rt)) / 
      (1 + exp(-global_par$ss_dgwt/rt))
  Pu = 1-Pf
  # print(paste0("StabilityPCA: folded state = ", round(Pf*100,1), "%"))
  global_par_dt = data.table(experiment = "stabilityPCA",
                              state = c("f"),
                              P = c(Pf))
  #bindingPCA
  Pf = exp(-global_par$sb_dgwt/rt) / 
      (1 + exp(-global_par$sb_dgwt/rt) + exp(-(global_par$sb_dgwt + global_par$b_dgwt)/rt))
  Pfb = exp(-(global_par$sb_dgwt + global_par$b_dgwt)/rt) / 
      (1 + exp(-global_par$sb_dgwt/rt) + exp(-(global_par$sb_dgwt + global_par$b_dgwt)/rt))
  Pu = 1 - Pf - Pfb

  global_par_dt = rbind(global_par_dt,
    data.table(experiment = "bindingPCA",
                              state = c("f", "fb"),
                              P = c(Pf, Pfb)))
  global_par_dt[, experiment := factor(experiment, levels = c("stabilityPCA", "bindingPCA"))]
  global_par_dt[, state := factor(state, levels = c("f", "fb"))]
  P <- ggplot(global_par_dt, aes(x = experiment, group = experiment, y = P, fill = state)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1")
  ggsave(P, 
    file = paste0("results/dG/GRB2_dG_sdgwt_method", method, 
    ifelse(predict_binding0 == TRUE, "_b0", ""), "_stateprobs.pdf"),
    width = 4, height = 4)


  s_ddg <- data.table(variant = sapply(X = grep("s_ddg", names(best_model)),
                                       FUN = function(X) {strsplit(names(best_model)[X], "_")[[1]][1]}),
                      s_ddg = best_model[, unlist(.SD), .SDcols = grep("s_ddg", names(best_model))])
  b_ddg <- data.table(variant = sapply(X = grep("b_ddg", names(best_model)),
                                       FUN = function(X) {strsplit(names(best_model)[X], "_")[[1]][1]}),
                      b_ddg = best_model[, unlist(.SD), .SDcols = grep("b_ddg", names(best_model))])

  table_dg <- merge(s_ddg, b_ddg, all = T)

  all_data <- fread(paste0("processed_data/", dataset_name, ".txt"))
  all_data[, s_ddg1 := table_dg[variant == id1,s_ddg],id1]
  all_data[, b_ddg1 := table_dg[variant == id1,b_ddg],id1]
  all_data[, s_ddg2 := table_dg[variant == id2,s_ddg],id2]
  all_data[, b_ddg2 := table_dg[variant == id2,b_ddg],id2]

  # calculate stabilityPCA fitness from dG values
  all_data[, s_fitness_pred := function_folding_dG2F(
    s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
    s_dgwt = global_par$ss_dgwt,
    s_f0 = ifelse("s_f0" %in% names(global_par), 
                  global_par$s_f0, 0), 
    s_fwt = ifelse("s_fwt" %in% names(global_par), 
                   global_par$s_fwt, 1)
  ), 
  .(s_ddg1, s_ddg2)
  ]

  # calculate bindingPCA fitness with estimated b_ddG values
  all_data[, 
           b_fitness_pred := function_binding_dG2F(
             b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
             s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s_dgwt = global_par$sb_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ),
           .(b_ddg1, s_ddg1, b_ddg2, s_ddg2)
           ]

  # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[, 
           b_fitness_pred_bdg0 := function_binding_dG2F(
             b_ddg = 0,
             s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
             b_dgwt = global_par$b_dgwt,
             s_dgwt = global_par$sb_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ), 
           .(s_ddg1, s_ddg2)
           ]

  # calculate bindingPCA fitness if s_ddG = 0; i.e. global relationship governed by b_ddG
  all_data[, 
           b_fitness_pred_sdg0 := function_binding_dG2F(
             b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
             s_ddg = 0,
             b_dgwt = global_par$b_dgwt,
             s_dgwt = global_par$sb_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           ), 
           .(b_ddg1, b_ddg2)
           ]

  ## relationship between parameters and fitness
  # load structural property file to color residues by type
  variant_structuralproperties = fread(paste0("processed_data/GRB2_variant_structuralproperties.txt"))
  all_data[,RSA_unbound := variant_structuralproperties[Pos == Pos1,RSA_unbound],Pos1]
  all_data[,HAmin_ligand := variant_structuralproperties[Pos == Pos1,HAmin_ligand],Pos1]
  all_data[RSA_unbound <= 20 & HAmin_ligand > 5,res_type1 := "core"]
  all_data[RSA_unbound <= 20 & HAmin_ligand <= 5,res_type1 := "core_bind"]
  all_data[RSA_unbound > 20 & HAmin_ligand > 5,res_type1 := "surf"]
  all_data[RSA_unbound > 20 & HAmin_ligand <= 5,res_type1 := "surf_bind"]
  all_data[Nmut == 1, .N, res_type1]

  ######################
  ### R2 on test-set ###
  ######################

  
  # if (method != 3) {
  #   all_data[test_set == TRUE, Nsid_min := min(all_data[Nmut == 1 & id1 %in% .SD[,1], Nid + Nsid],
  #                                           all_data[Nmut == 1 & id1 %in% .SD[,2], Nid + Nsid]), 
  #                           .(id1, id2),
  #                           .SDcols = c("id1","id2")]
  # } else {
  #   all_data[test_set == TRUE, Nsid_min := min(all_data[Nmut == 1 & id1 %in% .SD[,1], Nid + Nbid],
  #                                           all_data[Nmut == 1 & id1 %in% .SD[,2], Nid + Nbid]), 
  #                           .(id1, id2),
  #                           .SDcols = c("id1","id2")]
  # }
  # all_data[test_set == TRUE, Nbid_min := min(all_data[Nmut == 1 & id1 %in% .SD[,1], Nid + Nbid],
  #                                           all_data[Nmut == 1 & id1 %in% .SD[,2], Nid + Nbid]), 
  #                           .(id1, id2),
  #                           .SDcols = c("id1","id2")]

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

  ps <- ggplot(all_data[test_set == TRUE],aes(s_fitness_pred, s_fitness)) +
    geom_point() +
    geom_abline(color = "red") +
    geom_smooth() +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    geom_label(inherit.aes = F,
               data = data.table(x = 0.1,
                                 y = c(1, 0.9),
                                 label = c(paste0("R2 = ",round(R2_stab,2), " (", round(R2_stab/R2_stab_nontech*100,1), "% of total)"),
                                           paste0("R2_weighted = ",round(R2_stab_w,2), " (", round(R2_stab_w/R2_stab_nontech_w*100,1), "% of total)"))),
               aes(x,y,label = label),
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
  pb <- ggplot(all_data[test_set == TRUE],aes(b_fitness_pred, b_fitness)) +
    geom_point() +
    geom_abline(color = "red") +
    geom_smooth() +
    scale_x_continuous(breaks = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    geom_label(inherit.aes = F,
               data = data.table(x = 0.1,
                                 y = c(1, 0.9),
                                 label = c(paste0("R2 = ",round(R2_bind,2), " (", round(R2_bind/R2_bind_nontech*100,1), "% of total)"),
                                           paste0("R2_weighted = ",round(R2_bind_w,2), " (", round(R2_bind_w/R2_bind_nontech_w*100,1), "% of total)"))),
               aes(x,y,label = label),
               hjust = 0, size = 2.5) +
    labs(y = "fitness bindingPCA",
         x = "predicted fitness")
  P = grid.arrange(ps,pb,nrow=1)
  ggsave(P, 
    file = paste0("results/dG/GRB2_dG_sdgwt_method", method, 
    ifelse(predict_binding0 == TRUE, "_b0", ""), "_R2.pdf"),
    width = 9.5, height = 4)


  ###############################################################
  ### dG value relation with each other and to fitness values ###
  ###############################################################
  
  ds_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(s_ddg1 + global_par$ss_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$ss_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = c(0.85, 0.75)) +
    labs(x = "dG folding", color = "")

  ds_fs <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
      aes(s_ddg1 + global_par$ss_dgwt, s_fitness), color = 'black') +
    # geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], 
      # aes(s_ddg1 + global_par$s_dgwt, s_fitness, color = res_type1)) +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], 
      aes(s_ddg1 + global_par$ss_dgwt, s_fitness, color = res_type1)) +
    geom_line(data = all_data[Nmut==1 & !is.na(b_fitness)],
      aes(s_ddg1 + global_par$ss_dgwt,s_fitness_pred),color = "black") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$s_f0,global_par$s_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$ss_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "fitness stabilityPCA", color = "") + theme(legend.position = "none")


  ds_fb <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
        aes(s_ddg1 + global_par$sb_dgwt, b_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], aes(s_ddg1 + global_par$sb_dgwt, 
                                                             b_fitness, color = res_type1)) +
    geom_line(data = all_data[Nmut==1 & !is.na(b_fitness)],
              aes(s_ddg1 + global_par$sb_dgwt, 
                  b_fitness_pred_bdg0),color = "black") +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.15)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$ss_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "fitness bindingPCA", color = "") + theme(legend.position = "none")

  sVSb <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s_ddg1 + global_par$ss_dgwt, 
                                                            b_ddg1 + global_par$b_dgwt,
                                                            color = res_type1)) +
    geom_density2d() +
    # geom_point(alpha=0.5) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_y_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$ss_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "") + theme(legend.position = "none")

  db_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    coord_flip() +
    labs(x = "dG binding") + theme(legend.position = "none")

  db_fb <- ggplot() +
    geom_density2d(data = all_data[Nmut==1 & !is.na(b_fitness)],
        aes(b_ddg1  + global_par$b_dgwt, b_fitness), color = 'black') +
    geom_point(data = all_data[Nmut==1 & !is.na(b_fitness)], aes(b_ddg1 + global_par$b_dgwt, 
                                                             b_fitness, color = res_type1)) +
    geom_line(data = all_data[Nmut==1 & !is.na(b_fitness)],
              aes(b_ddg1 + global_par$b_dgwt, 
                  b_fitness_pred_sdg0),color = "black") +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1.15)) +
    scale_x_continuous(breaks = seq(-2,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG binding", y = "fitness bindingPCA", color = "") + theme(legend.position = "none")
  
  P <- grid.arrange(ds_dist, ds_fs, ds_fb, sVSb, db_dist, db_fb,
                   nrow = 2)
  ggsave(P, file = paste0("results/dG/GRB2_dG_sdgwt_method", method, 
    ifelse(predict_binding0 == TRUE, "_b0", ""), "_dG_fitness.pdf"), 
         width = 11, height = 7)


  ##############################
  ### s_fitness vs b_fitness ###
  ##############################
  p <- ggplot(all_data[is.na(id2) & !is.na(b_fitness)]) +
    geom_point(aes(
      b_fitness, 
      s_fitness,
      color = res_type1
    )) +
    geom_line(aes(
      b_fitness_pred_bdg0, 
      s_fitness_pred), 
      color = "black",
      size = 1.5
    ) +
    geom_segment(aes(
      x = b_fitness_pred_bdg0, 
      y = s_fitness_pred, 
      xend = b_fitness, 
      yend = s_fitness), 
      alpha = 0.3
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "fitness bindingPCA", 
      y = "fitness stabilityPCA",
      color = "residue type"
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0.15, 1.15)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25))
  ggsave(plot = p, 
         file = paste0("results/dG/GRB2_dG_sdgwt_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_fitness_globalrelationship.pdf"), 
         width = 6, 
         height = 4.5)



  ##### uncertainty of fitness values
  ### bootstrapped fitness values
  boot_fitness_models <- fread(paste0("processed_data/", 
                                      dataset_name, "_method", method, 
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

  table_dg <- merge(table_dg, boot_par_sd[grep("s_ddg", parameter),.(variant, s_ddg_sd = sd)], all = T)
  table_dg <- merge(table_dg, boot_par_sd[grep("b_ddg", parameter),.(variant, b_ddg_sd = sd)], all = T)
  table_dg
  # how does sd relate to # doubles seen?
  all_data[!is.na(s_fitness),N_stability := .N + sum(all_data[!is.na(s_fitness),id2] == id1, na.rm = T), id1]
  all_data[!is.na(b_fitness),N_binding := .N + sum(all_data[!is.na(b_fitness),id2] == id1, na.rm = T), id1]
  table_dg = merge(table_dg,all_data[Nmut == 1,.(variant = id, N_stability, N_binding)], all.x = T)

  p1 = ggplot(table_dg, aes(s_ddg, b_ddg, color = s_ddg_sd)) +
    geom_point() +
    # geom_errorbarh(aes(xmin = s_ddg - 3*s_ddg_sd, xmax = s_ddg + 3*s_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding", y = "dG binding", color = "sd(s_ddg est.)")

  p2 = ggplot(table_dg, aes(s_ddg, b_ddg, color = b_ddg_sd)) +
    geom_point() +
    # geom_errorbar(aes(ymin = b_ddg - 3*b_ddg_sd, ymax = b_ddg + 3*b_ddg_sd)) +
    scale_color_gradient(high = "black", low = "red") +
    labs(x = "dG folding", y = "dG binding", color = "sd(b_ddg est.)")

  p = grid.arrange(p1, p2, nrow = 1)
  ggsave(plot = p, 
         file = paste0("results/dG/GRB2_dG_method", method, 
                ifelse(predict_binding0 == TRUE, "_b0", ""), 
                "_dGpar_bootsd.pdf"), 
         width = 8, 
         height = 3.5)


  ##########################
  ### profile likelihood ###

  # pll_bounds <- fread(paste0("processed_data/", 
  #                            dataset_name, "_method", method,
  #                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
  #                            "_parameter_pll_df1.txt"))
  # pll_bounds[grep("ddg", par), variant := strsplit(par,"_")[[1]][1], par]

  if (method %in% c(2,3)) {


  pll_bounds <- fread(paste0("processed_data/", 
                             dataset_name, "_method", method,
                             ifelse(predict_binding0 == TRUE, "_b0", ""), 
                             "_parameter_pll_df", ncol(best_model)-2, ".txt"))
  pll_bounds[grep("ddg", par), variant := strsplit(par,"_")[[1]][1], par]
  table_pll_dg <- merge(pll_bounds[grep("s_ddg", par),
                               .(variant, 
                                 s_ddg = value, 
                                 s_ddg_upper = upper, 
                                 s_ddg_lower = lower)],
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

  p1 <- ggplot() +
    geom_density2d(data = table_pll_dg, aes(s_ddg + global_par$s_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & !(id1 %in% table_pll_dg[ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0),variant])],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & id1 %in% table_pll_dg[ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0),variant]],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "sig. dG folding != wt")

  p2 <- ggplot() +
    geom_density2d(data = table_pll_dg, aes(s_ddg + global_par$s_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & !(id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant])],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "sig. dG binding != wt")

  p3 <- ggplot() +
    geom_density2d(data = table_pll_dg, aes(s_ddg + global_par$s_dgwt, b_ddg + global_par$b_dgwt), color = "black") +
    geom_point(data = all_data[Nmut == 1 & 
        (!(id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]) |
        !(id1 %in% table_pll_dg[ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0),variant]))],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt), alpha = 0.2) + 
    geom_point(data = all_data[Nmut == 1 & 
        id1 %in% table_pll_dg[ifelse(b_ddg > 0, b_ddg_lower > 0, b_ddg_upper < 0),variant]  &
        id1 %in% table_pll_dg[ifelse(s_ddg > 0, s_ddg_lower > 0, s_ddg_upper < 0),variant]],
      aes(s_ddg1 + global_par$s_dgwt, b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "sig. dG binding != wt & dG folding != wt")

  p = grid.arrange(grobs = list(p1, p2, p3), 
      nrow = 2,
      widths = c(1.25,0.25,1),
      layout_matrix = rbind(c(1,2,2),
                            c(3,3,NA)))
  ggsave(plot = p, 
         file = paste0("results/dG/GRB2_dG_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_dGpar_pll", ncol(best_model)-2, ".pdf"), 
         width = 10, 
         height = 7)




 ## pll recordings for s_dgwt
  pll_models <- fread(paste0("processed_data/", 
                             dataset_name, "_method", method,
                             ifelse(predict_binding0 == TRUE, "_b0", ""), 
                             "_parameter_pll_df", ncol(best_model)-2, "_s_dgwt.txt"))
  global_par_pll <- pll_models[, .SD, .SDcols = grep("^[sb]", names(pll_models))]
  P <- ggpairs(global_par_pll)
  ggsave(plot = P, 
         file = paste0("results/dG/GRB2_dG_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_dGpar_pll", ncol(best_model)-2, "_s_dgwt.pdf"), 
         width = 8, 
         height = 8)

## pll recordings for b_dgwt
  pll_models <- fread(paste0("processed_data/", 
                             dataset_name, "_method", method,
                             ifelse(predict_binding0 == TRUE, "_b0", ""), 
                             "_parameter_pll_df", ncol(best_model)-2, "_b_dgwt.txt"))
  global_par_pll <- pll_models[, .SD, .SDcols = grep("^[sb]", names(pll_models))]
  P <- ggpairs(global_par_pll)
  ggsave(plot = P, 
         file = paste0("results/dG/GRB2_dG_method", method, 
                        ifelse(predict_binding0 == TRUE, "_b0", ""), 
                        "_dGpar_pll", ncol(best_model)-2, "_b_dgwt.pdf"), 
         width = 8, 
         height = 8)


  }

  # all_data[Nmut == 1, s_ddg_upper := pll_bounds[grepl("s_ddg", par)][variant %in% unique(id1), unique(upper)], id1]
  # all_data[Nmut == 1, s_ddg_lower := pll_bounds[grepl("s_ddg", par)][variant %in% unique(id1), unique(lower)], id1]
  # all_data[Nmut == 1, b_ddg_upper := pll_bounds[grepl("b_ddg", par)][variant %in% unique(id1), unique(upper)], id1]
  # all_data[Nmut == 1, b_ddg_lower := pll_bounds[grepl("b_ddg", par)][variant %in% unique(id1), unique(lower)], id1]

  # # pll_bounds[grep("ddg", par), Pos := gsub("[A-Z]","",variant), variant]
  # # table_pll_dg <- merge(pll_bounds[grep("s_ddg", par),
  # #                              .(variant, Pos,                                
  # #                                s_ddg = value, 
  # #                                s_ddg_upper = upper, 
  # #                                s_ddg_lower = lower)],
  # #                   pll_bounds[grep("b_ddg", par),
  # #                              .(variant, Pos,
  # #                                b_ddg = value, 
  # #                                b_ddg_upper = upper, 
  # #                                b_ddg_lower = lower)], all = T)
  # p1 = ggplot(table_pll_dg, aes(s_ddg, b_ddg, color = s_ddg_upper - s_ddg_lower)) +
  #   geom_point() +
  #   # geom_errorbarh(aes(xmin = s_ddg_lower, xmax = s_ddg_upper)) +
  #   labs(x = "dG folding", y = "dG binding")

  # p2 = ggplot(table_pll_dg, aes(s_ddg, b_ddg, color = b_ddg_upper - b_ddg_lower)) +
  #   geom_point() +
  #   # geom_errorbar(aes(ymin = b_ddg_lower, ymax = b_ddg_upper)) +
  #   labs(x = "dG folding", y = "dG binding")

  # p = grid.arrange(p1, p2, nrow = 1)
  # ggsave(plot = p, 
  #        file = paste0("results/dG/GRB2_dG_method", method, 
  #                       ifelse(predict_binding0 == TRUE, "_b0", ""), 
  #                       "_dGpar_pll", ncol(best_model)-2, ".pdf"), 
  #        width = 8, 
  #        height = 3.5)


  # # table_pll_dg[s_ddg > median(s_ddg), s_ddg_z := (s_ddg - median(s_ddg)) / abs(s_ddg_lower - s_ddg)]
  # # table_pll_dg[s_ddg < median(s_ddg), s_ddg_z := (s_ddg - median(s_ddg)) / abs(s_ddg_upper - s_ddg)]
  # s_ddg_median <- all_data[Nmut == 1, median(s_ddg1, na.rm = T)]
  # all_data[Nmut == 1 & s_ddg1 > s_ddg_median, 
  #   s_ddg_z := max(c(s_ddg_median, s_ddg_lower)) - 
  #     s_ddg_median, 
  #   id]
  # all_data[Nmut == 1 & s_ddg1 < s_ddg_median, 
  #   s_ddg_z := min(c(s_ddg_median, s_ddg_upper)) - 
  #     s_ddg_median, 
  #   id]
  # # table_pll_dg[b_ddg > 0, b_ddg_z := b_ddg / abs(b_ddg_lower - b_ddg)]
  # # table_pll_dg[b_ddg < 0, b_ddg_z := b_ddg / abs(b_ddg_upper - b_ddg)]
  # b_ddg_median <- all_data[Nmut == 1, median(b_ddg1, na.rm = T)]
  # all_data[Nmut == 1 & b_ddg1 > b_ddg_median, 
  #   b_ddg_z := max(c(b_ddg_median, b_ddg_lower)) - 
  #     b_ddg_median, 
  #   id]
  # all_data[Nmut == 1 & b_ddg1 < b_ddg_median, 
  #   b_ddg_z := min(c(b_ddg_median, b_ddg_upper)) - 
  #     b_ddg_median, 
  #   id]
  
  
  # p1 = ggplot(all_data[Nmut == 1], aes(s_ddg_z, b_ddg_z, color = res_type1)) +
  #   geom_density2d() + 
  #   geom_point() +
  #   labs(x = "dG folding 'zscore'", y = "dG binding 'zscore'")
  # ggsave(plot = p1, 
  #        file = paste0("results/dG/GRB2_dG_method", method, 
  #                       ifelse(predict_binding0 == TRUE, "_b0", ""), 
  #                       "_dGpar_pll", ncol(best_model)-2, "_zscore.pdf"), 
  #        width = 8, 
  #        height = 3.5)


  # all_data[,bin_RSA := findInterval(RSA_unbound,seq(0, 100, 10))]
  # all_data[,bin_dist := findInterval(HAmin_ligand,seq(2.5, 20, 2.5))]
  

  # all_data[Nmut == 1,.(b_ddg = mean(b_ddg1, na.rm = T),
  #   b_ddg_z = mean(b_ddg_z, na.rm = T),
  #   s_ddg = mean(s_ddg1, na.rm = T),
  #   s_ddg_z = mean(s_ddg_z, na.rm = T),
  #   dist = unique(HAmin_ligand),
  #   RSA = unique(RSA_unbound)), Pos1][order(s_ddg_z)]


  # ### Rsa vs dGf, dist vs dGb
  # p1 <- ggplot(all_data[Nmut == 1], aes(RSA_unbound, group = bin_RSA ,s_fitness)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "RSA unbound", y = "fitness stabilityPCA")
  # p2 <- ggplot(all_data[Nmut == 1], aes(RSA_unbound, group = bin_RSA ,s_ddg1)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "RSA unbound", y = "dG folding")
  # p3 <- ggplot(all_data[Nmut == 1], aes(RSA_unbound, group = bin_RSA ,s_ddg_z)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "RSA unbound", y = "dG folding 'zscore'")
  # p <- grid.arrange(p1, p2, p3, nrow = 2)
  # ggsave(plot = p, 
  #        file = paste0("results/dG/GRB2_dG_method", method, 
  #                       ifelse(predict_binding0 == TRUE, "_b0", ""), 
  #                       "_dGpar_pll", ncol(best_model)-2, "_dGz_RSA.pdf"), 
  #        width = 8, 
  #        height = 3.5)

  # p1 <- ggplot(all_data[Nmut == 1], aes(HAmin_ligand, group = bin_dist ,b_fitness)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "dist to ligand", y = "fitness bindingPCA")
  # p2 <- ggplot(all_data[Nmut == 1], aes(HAmin_ligand, group = bin_dist ,b_fitness - s_fitness)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "dist to ligand", y = "delta fitness")
  # p3 <- ggplot(all_data[Nmut == 1], aes(HAmin_ligand, group = bin_dist ,b_ddg1)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "dist to ligand", y = "dG binding")
  # p4 <- ggplot(all_data[Nmut == 1], aes(HAmin_ligand, group = bin_dist ,b_ddg_z)) +
  #   geom_point(aes(alpha = 0.01)) +
  #   geom_boxplot(color = "red", outlier.shape = NA, alpha = 0.1) +
  #   labs(x = "dist to ligand", y = "dG binding 'zscore'")
  # p <- grid.arrange(p1, p2, p3, p4, nrow = 2)
  # ggsave(plot = p, 
  #        file = paste0("results/dG/GRB2_dG_method", method, 
  #                       ifelse(predict_binding0 == TRUE, "_b0", ""), 
  #                       "_dGpar_pll", ncol(best_model)-2, "_dGz_distligand.pdf"), 
  #        width = 8, 
  #        height = 3.5)


# malagrino <- fread("dataset/GRB2/Malagrino2019_dG_binding.txt", skip = 2)
# dt_malagrino <- data.table(id0 = malagrino$V1,
#   KD = sapply(X = 1:nrow(malagrino),
#     FUN = function(X) {as.numeric(strsplit(malagrino[X,V4]," ")[[1]][1])}),
#   KD_sd = sapply(X = 1:nrow(malagrino),
#     FUN = function(X) {as.numeric(strsplit(malagrino[X,V4]," ")[[1]][3])}))
# dt_malagrino[, id := paste0(strsplit(id0,"")[[1]][2:nchar(id0)], collapse = ""),id0]
# all_data[Nmut == 1, KD_malagrino := dt_malagrino[id %in% id1, KD], id1]
# all_data[Nmut == 1 & !is.na(KD_malagrino),.(id1, b_ddg1, KD_malagrino, HAmin_ligand)]

}
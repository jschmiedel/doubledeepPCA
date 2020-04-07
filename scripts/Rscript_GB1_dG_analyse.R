require(data.table)
require(ggplot2)
require(GGally)
require(gridExtra)
theme_set(theme_bw(base_size = 9))

setwd("doubledeepPCA")

filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)



dataset_name_vec = c("GB1_dG_dataset","GB1_dG_dataset2","GB1_dG_dataset2reg","GB1_dG_dataset2","GB1_dG_dataset2reg")
approach_vec <- c(1, 1, 1, 2, 2)
dataset_size_vec <- c("1", "1", "1",  "1", "1", "0p33", "0p1", "0p033", "0p01","0p0033")
predict_binding0_vec <- rep(FALSE, length(dataset_size_vec))

R2_stab = c()
R2_bind_w = c()
Pf = c()
Pfb = c()

for (i in 5) {
  dataset_name <- dataset_name_vec[i]
  dataset_size <- dataset_size_vec[i]
  predict_binding0 <- predict_binding0_vec[i]
  approach <- approach_vec[i]

  initialmodels <- boot_par_models <- fread(paste0("processed_data/dG/", 
                                  dataset_name, "_approach", approach,
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_size", dataset_size,
                                  "_initialmodels.txt"))

  filename = paste0("processed_data/dG/", 
                                  dataset_name, "_approach", approach,
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_size", dataset_size,
                                  "_boot_startpar.txt")
  if (file.exists(filename)) {
    boot_par_models <- fread(filename)
    all_models <- rbind(initialmodels,
      boot_par_models)
  } else {
    all_models <- initialmodels
  }
  all_models[, objective := log10(objective)]

  best_model <- all_models[which.min(objective)]
  global_par <- best_model[, .SD, .SDcols = grep("^[sb]",names(best_model))]

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
  all_data[, s_fitness_pred := function_folding_dg2logf(
              s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
              s_dgwt = global_par$s_dgwt,
              s_f0 = global_par$b_f0, 
              s_fwt = global_par$b_fwt
            ), 
    .(s_ddg1, s_ddg2)
  ]

  # calculate bindingPCA fitness with estimated b_ddG values
  all_data[, b_fitness_pred := function_binding_dg2logf(
               b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
               s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
               b_dgwt = global_par$b_dgwt,
               s_dgwt = global_par$s_dgwt,
               b_f0 = global_par$b_f0,
               b_fwt = global_par$b_fwt
             ),
    .(b_ddg1, s_ddg1, b_ddg2, s_ddg2)
  ]

  # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[, b_fitness_pred_bdg0 := function_binding_dg2logf(
               b_ddg = 0,
               s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
               b_dgwt = global_par$b_dgwt,
               s_dgwt = global_par$s_dgwt,
               b_f0 = global_par$b_f0,
               b_fwt = global_par$b_fwt
             ), 
   .(s_ddg1, s_ddg2)
  ]

  # calculate bindingPCA fitness if s_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[, b_fitness_pred_sdg0 := function_binding_dg2logf(
               b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
               s_ddg = 0,
               b_dgwt = global_par$b_dgwt,
               s_dgwt = global_par$s_dgwt,
               b_f0 = global_par$b_f0,
               b_fwt = global_par$b_fwt
             ), 
   .(b_ddg1, b_ddg2)
  ]

  ## relationship between parameters and fitness
  # load structural property file to color residues by type
  variant_structuralproperties = fread(paste0("processed_data/GB1_structural_data.txt"))
  all_data[,res_type1 := variant_structuralproperties[Pos %in% Pos1,unique(type)],Pos1]
  all_data[,RSA_unbound := variant_structuralproperties[Pos %in% Pos1,unique(RSA_unbound)],Pos1]
  all_data[,HAmin_ligand := variant_structuralproperties[Pos %in% Pos1,unique(HAmin_ligand)],Pos1]

  all_data[, s_ddg_measured1 := 
    variant_structuralproperties[id == id1, s_ddg_measured],id1]
  all_data[, s_ddg_measured1_sd := 
    variant_structuralproperties[id == id1, s_ddg_measured_sdsmooth],id1]
  all_data[, s_ddg_measured2 := 
    variant_structuralproperties[id == id2, s_ddg_measured],id2]

  ##############################
  ### plot model performance ###
  ##############################
  all_models_globalpar <- all_models[, 
        cbind(.SD, type = "initial"), 
        .SDcols = grep("^[osb]",names(all_models))]
  P <- ggpairs(all_models_globalpar[objective < 7],
    columns = grep("^[osb]",names(all_models_globalpar)))
  ggsave(P,
    file = paste0("results/dG/",dataset_name, "_approach", approach,
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_size", dataset_size,
      "_modelperformance.pdf"),
    width = 10,
    height = 10)

  ###########################
  ### state probabilities ###
  ###########################

  rt = 1.99e-3 * 310.15
  Pf[i] = exp(-global_par$s_dgwt/rt) / 
      (1 + exp(-global_par$s_dgwt/rt) + exp(-(global_par$s_dgwt + global_par$b_dgwt)/rt))
  Pfb[i] = exp(-(global_par$s_dgwt + global_par$b_dgwt)/rt) / 
      (1 + exp(-global_par$s_dgwt/rt) + exp(-(global_par$s_dgwt + global_par$b_dgwt)/rt))
  # Pu = 1 - Pf - Pfb

  global_par_dt = data.table(experiment = "mRNA display",
                              state = c("folded", "folded & bound"),
                              probability = c(Pf[i], Pfb[i]))
  global_par_dt[, experiment := factor(experiment)]
  global_par_dt[, state := factor(state, levels = c("folded", "folded & bound"))]
  P <- ggplot(global_par_dt, aes(x = experiment, group = experiment, y = probability, fill = state)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1") + theme_classic()
  ggsave(P, 
    file = paste0("results/dG/",dataset_name, "_approach", approach,
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_size", dataset_size, "_stateprobs.pdf"),
    width = 4, height = 3)


  ######################
  ### R2 on test-set ###
  ######################

  R2_stab[i] = cor(all_data[Nmut == 1 & !is.na(s_ddg_measured1) & b_fitness > -4.5 & s_ddg_measured1 < 4,
    .(s_ddg1,s_ddg_measured1)])[2,1]
  print(R2_stab)

  # cov.wt(all_data[Nmut == 1 & !is.na(s_ddg_measured1_sd) & b_fitness > -4.5 & s_ddg_measured1 > -4,
  #   .(s_ddg1, s_ddg_measured1)],
  #   wt = all_data[Nmut == 1 & !is.na(s_ddg_measured1_sd) & b_fitness > -4.5 & s_ddg_measured1 > -4,
  #    1 / s_ddg_measured1_sd^2], cor = T)$cor[2, 1]^2
X = all_data[Nmut == 1]
X[, s_ddg_shown := s_ddg_measured1]
X[s_ddg_measured1 == 4, s_ddg_shown := rnorm(.N, 4, 0.1)]
  ps <- ggplot(X,
    aes(s_ddg_shown, s_ddg1)) +
    # geom_point(aes(color = b_fitness < -4.5)) +
    geom_point(aes(color = s_ddg_measured1 < 4 & b_fitness > -4.5 & (s_ddg1-s_ddg_measured1) > 1), alpha = 0.5) +
    geom_abline(color = "black", slope = 1) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = seq(-4, 10, 2)) +
    scale_y_continuous(breaks = seq(-4, 10, 2)) +
    scale_color_brewer(palette = "Set1") +
    geom_label(inherit.aes = F,
               data = data.table(x = -2,
                                 y = 7,
                                 label = c(paste0("R = ",round(R2_stab[i],2)))),
               aes(x,y,label=label),
               hjust = 0) +
    labs(y = "predicted ddG folding (-dG)",
         x = "measured ddG folding", 
         color = "")

  R2_bind = cor(all_data[test_set == TRUE, 
    .((b_fitness), (b_fitness_pred))], use = "pairwise.complete")[2, 1]^2
  R2_bind_w[i] = cov.wt(all_data[test_set == TRUE & b_fitness != 0 & b_fitness_pred != 0,
    .((b_fitness), (b_fitness_pred))],
    wt = all_data[test_set == TRUE, 1 / b_sigma^2], cor = T)$cor[2, 1]^2
  x = matrix(all_data[test_set == TRUE, (b_fitness) + 
    rnorm(100 * .N, mean = 0, sd = b_sigma)],
             nrow = all_data[test_set == TRUE, .N], ncol = 100)
  tmp = cor(x)
  R2_bind_nontech = mean(tmp[lower.tri(tmp)])^2
  R2_bind / R2_bind_nontech
  tmp = cov.wt(x, wt = all_data[test_set == TRUE, 1 / b_sigma^2], cor = T)
  R2_bind_nontech_w = mean(tmp$cor[lower.tri(tmp$cor)])^2
  R2_bind_w[i] / R2_bind_nontech_w
  
  pb <- ggplot(all_data[test_set == TRUE], aes(x = (b_fitness_pred), y = b_fitness)) +
    geom_point() +
    geom_abline(color = "red") +
    geom_smooth() +
    scale_x_continuous(breaks = seq(-5, 1, 1)) +
    scale_y_continuous(breaks = seq(-5, 1, 1)) +
    geom_label(inherit.aes = F,
               data = data.table(x = -6,
                 y = 1.5,
                 label = paste0("R2_weighted = ", round(R2_bind_w[i], 2), " (", round(R2_bind_w [i]/ R2_bind_nontech_w * 100, 1), "% of total)")),
               aes(x, y, label = label),
               hjust = 0, size = 2.5) +
    labs(y = "fitness",
         x = "predicted fitness")
  
  P = grid.arrange(grobs = list(ps, pb),
    widths = c(1.5, 1),
    nrow = 1)
  ggsave(P, 
    file = paste0("results/dG/",dataset_name, "_approach", approach,
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_size", dataset_size, 
      "_R2.pdf"),
    width = 8, height = 4)

  #understand which positions differ in s_ddg from measurement
  # destabilize folded state more in DMS than in vitro
  all_data[Nmut == 1 & s_ddg_measured1 < 4 & b_fitness > -4.5,
    .(N = sum((s_ddg1-s_ddg_measured1) > 1.25),
     frac = round(sum((s_ddg1-s_ddg_measured1) > 1.25)/.N,2),
     type = unique(res_type1)), Pos1][N>0][order(Pos1)]

  # destabilize folded state less in DMS than in vitro
  all_data[Nmut == 1 & b_fitness > -4.5 & (s_ddg1-s_ddg_measured1) < -1.5,
    .(.N, type = unique(res_type1)), Pos1]

    all_data[Nmut == 1 & b_fitness > -4.5 & !is.na(s_ddg_measured1),
    .(N = sum((s_ddg1-s_ddg_measured1) < -1),
     frac = round(sum((s_ddg1-s_ddg_measured1) < -1)/.N,2),
     type = unique(res_type1)), Pos1][N>0][order(Pos1)]

  ###############################################################
  ### dG value relation with each other and to fitness values ###
  ###############################################################
  # ds_fs <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s_ddg1 + global_par$s_dgwt, 
  #                                                            s_fitness, color = res_type1)) +
  #   geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
  #             aes(s_ddg1 + global_par$s_dgwt, 
  #                 s_fitness_pred),color = "black") +
  #   geom_point() +
  #   scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #   scale_x_continuous(breaks = seq(-2,6,2)) +
  #   geom_hline(yintercept = c(global_par$s_f0,global_par$s_fwt),linetype = 2) +
  #   geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
  #   labs(x = "dG folding", y = "fitness stabilityPCA", color = "") + theme(legend.position = "none")

  # unique(all_data[Nmut == 1,.( 
  #   res_type1, 
  #   b_fitness = mean(b_fitness),
  #   s_ddg = mean(s_ddg1),
  #   b_ddg = mean(b_ddg1)), Pos1])[order(Pos1)]

theme_set(theme_bw())
  ds_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(s_ddg1 + global_par$s_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = c(0.75, 0.75)) +
    labs(x = "dG folding", color = "")

  ds_fs <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(s_ddg1 + s_ddg2 + global_par$s_dgwt, 
                                                             s_fitness_pred), color = 'black') +
    # geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
    #           aes(s_ddg_measured1 + global_par$s_dgwt, 
    #               s_fitness_pred),color = "black") +
    geom_line(data = data.table(x = global_par$s_dgwt + seq(-5,10,0.1),
        y = function_folding_dg2logf(
             s_ddg = seq(-5,10,0.1),
             s_dgwt = global_par$s_dgwt,
             s_f0 = global_par$b_f0,
             s_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_point(data = all_data[Nmut == 1], aes(s_ddg1 + global_par$s_dgwt, 
                                                             s_fitness_pred, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "predicted folding fitness", color = "") + 
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set1")

  ds_fb <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(s_ddg1 + s_ddg2 + global_par$s_dgwt, 
                                                             b_fitness), color = 'black') +
    # geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
    #           aes(s_ddg_measured1 + global_par$s_dgwt, 
    #               b_fitness_pred_bdg0),color = "black") +
    geom_line(data = data.table(x = global_par$s_dgwt + seq(-5,10,0.1),
        y = function_binding_dg2logf(
             s_ddg = seq(-5,10,0.1),
             b_ddg = 0,
             s_dgwt = global_par$s_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_point(data = all_data[Nmut == 1], aes(s_ddg1 + global_par$s_dgwt, 
                                                             b_fitness, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "binding fitness", color = "") + 
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set1")

  sVSb <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s_ddg1 + global_par$s_dgwt, 
                                                            b_ddg1 + global_par$b_dgwt,
                                                            color = res_type1)) +
    geom_density2d() +
    # geom_point(alpha=0.5) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_y_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "") + theme(legend.position = "none")

  db_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    coord_flip() +
    labs(x = "dG binding") + theme(legend.position = "none")

  db_fb <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(b_ddg1 + b_ddg2 + global_par$b_dgwt, 
                                                             b_fitness), color = 'black') +
    geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
              aes(b_ddg1 + global_par$b_dgwt, 
                  b_fitness_pred_sdg0),color = "black") +
    geom_point(data = all_data[Nmut == 1], aes(b_ddg1 + global_par$b_dgwt, 
                                                             b_fitness, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG binding", y = "binding fitness", color = "") + theme(legend.position = "none")
  
  P <- grid.arrange(ds_dist, ds_fs, ds_fb, sVSb, db_dist, db_fb,
                   nrow = 2)
  ggsave(P, file = paste0("results/dG/",dataset_name, "_approach", approach,
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_size", dataset_size, "_dG_fitness.pdf"), 
         width = 10, height = 7)



####################################
#### using in vitro measured dGs####

  X = all_data[Nmut == 1]
  X[s_ddg_measured1 == 4, s_ddg_measured1 := rnorm(.N, 4, 0.251)]
  theme_set(theme_bw())
  ds_dist <- ggplot() +
    geom_density(data = X[Nmut == 1 & !is.na(b_fitness)],
      aes(s_ddg_measured1 + global_par$s_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = c(0.75, 0.75)) +
    labs(x = "dG folding", color = "")

  ds_fs <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(s_ddg_measured1 + s_ddg_measured2 + global_par$s_dgwt, 
                                                             s_fitness_pred), color = 'black') +
    # geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
    #           aes(s_ddg_measured1 + global_par$s_dgwt, 
    #               s_fitness_pred),color = "black") +
    geom_line(data = data.table(x = global_par$s_dgwt + seq(-5,10,0.1),
        y = function_folding_dg2logf(
             s_ddg = seq(-5,10,0.1),
             s_dgwt = global_par$s_dgwt,
             s_f0 = global_par$b_f0,
             s_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_point(data = X[Nmut == 1], aes(s_ddg_measured1 + global_par$s_dgwt, 
                                                             s_fitness_pred, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "predicted folding fitness", color = "") + 
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set1")

  ds_fb <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(s_ddg_measured1 + s_ddg_measured2 + global_par$s_dgwt, 
                                                             b_fitness), color = 'black') +
    # geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
    #           aes(s_ddg_measured1 + global_par$s_dgwt, 
    #               b_fitness_pred_bdg0),color = "black") +
    geom_line(data = data.table(x = global_par$s_dgwt + seq(-5,10,0.1),
        y = function_binding_dg2logf(
             s_ddg = seq(-5,10,0.1),
             b_ddg = 0,
             s_dgwt = global_par$s_dgwt,
             b_dgwt = global_par$b_dgwt,
             b_f0 = global_par$b_f0,
             b_fwt = global_par$b_fwt
           )),
      aes(x,y),color = "black") +
    geom_point(data = X[Nmut == 1], aes(s_ddg_measured1 + global_par$s_dgwt, 
                                                             b_fitness, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-10,6,2)) +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "binding fitness", color = "") + 
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set1")

  sVSb <- ggplot(all_data[Nmut==1 & !is.na(b_fitness)], aes(s_ddg_measured1 + global_par$s_dgwt, 
                                                            b_ddg1 + global_par$b_dgwt,
                                                            color = res_type1)) +
    geom_density2d() +
    # geom_point(alpha=0.5) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_y_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = global_par$s_dgwt,linetype = 2) +
    geom_hline(yintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG folding", y = "dG binding", color = "") + theme(legend.position = "none")

  db_dist <- ggplot() +
    geom_density(data = all_data[Nmut == 1 & !is.na(b_fitness)],
      aes(b_ddg1 + global_par$b_dgwt, color = res_type1)) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    coord_flip() +
    labs(x = "dG binding") + theme(legend.position = "none")

  db_fb <- ggplot() +
    geom_density2d(data = all_data[test_set == T], aes(b_ddg1 + b_ddg2 + global_par$b_dgwt, 
                                                             b_fitness), color = 'black') +
    geom_line(inherit.aes = F, data = all_data[Nmut==1 & !is.na(b_fitness)],
              aes(b_ddg1 + global_par$b_dgwt, 
                  b_fitness_pred_sdg0),color = "black") +
    geom_point(data = all_data[Nmut == 1], aes(b_ddg1 + global_par$b_dgwt, 
                                                             b_fitness, color = res_type1)) +
    scale_y_continuous(breaks = seq(-6, 1, 1)) +
    scale_x_continuous(breaks = seq(-6,6,1)) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = c(global_par$b_f0,global_par$b_fwt),linetype = 2) +
    geom_vline(xintercept = global_par$b_dgwt,linetype = 2) +
    labs(x = "dG binding", y = "binding fitness", color = "") + theme(legend.position = "none")
  
  P <- grid.arrange(ds_dist, ds_fs, ds_fb, sVSb, db_dist, db_fb,
                   nrow = 2)
  ggsave(P, file = paste0("results/dG/",dataset_name, "_approach", approach,
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_size", dataset_size, "_dGmeas_fitness.pdf"), 
         width = 10, height = 7)

#   # ##############################
#   # ### s_fitness vs b_fitness ###
#   # ##############################
  p <- ggplot() +
    geom_density2d(data = all_data,
      aes(b_fitness, s_fitness_pred), color = "black") +
    geom_point(data = all_data[is.na(id2) & !is.na(b_fitness)],
      aes(b_fitness, s_fitness_pred, color = res_type1)) +
    geom_line(data = all_data[is.na(id2) & !is.na(b_fitness)],
      aes(b_fitness_pred_bdg0, s_fitness_pred), color = "red", size = 2) +
    # geom_segment(aes(
    #   x = b_fitness_pred_bdg0, 
    #   y = s_fitness_pred, 
    #   xend = b_fitness, 
    #   yend = s_fitness), 
    #   alpha = 0.3
    # ) +
    labs(
      x = "binding fitness", 
      y = "predicted fitness stabilityPCA",
      color = "residue type"
    ) +
    scale_x_continuous() +
    scale_y_continuous() +
    scale_color_brewer(palette = "Set1")
  ggsave(plot = p, 
         file = paste0("results/dG/",dataset_name, "_approach", approach,
            ifelse(predict_binding0 == TRUE, "_b0", ""), 
            "_size", dataset_size, "_fitness_globalrelationship.pdf"), 
         width = 6, 
         height = 4.5)

}
dataset_size_vecnum = sapply(X = 1:length(dataset_size_vec), FUN = function(X) {gsub("p",".", dataset_size_vec[X])})
global_par_dt = data.table(dataset_size = dataset_size_vecnum,
                              state = rep(c("folded", "folded & bound"), each = length(Pf)),
                              probability = c(Pf, Pfb))
global_par_dt[, dataset_size := factor(dataset_size, levels = dataset_size_vecnum)]
global_par_dt[, state := factor(state, levels = c("folded", "folded & bound"))]
P <- ggplot(global_par_dt, aes(x = dataset_size, group = dataset_size, y = probability, fill = state)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1") + 
    theme_classic() +
    labs(x = "fraction of doubles")
ggsave(P, 
    file = paste0("results/dG/GB1_dG",
      ifelse(predict_binding0 == TRUE, "_b0", ""), 
      "_stateprobs_allsizes.pdf"),
    width = 4, height = 3)
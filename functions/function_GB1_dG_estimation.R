function_GB1_dG_estimation <- function(
  dataset_name = "GB1_dG_dataset",
  dataset_size = 1,
  approach = 1,
  iteration = 0,
  task,
  predict_binding0 = FALSE,
  maxit = 1e4,
  n_models = 100,
  n_cores = 1
) { 
  require(data.table)
  require(Matrix)
  require(tictoc)
  require(foreach)
  require(doMC)
  
  all_data <- fread(
    paste0("processed_data/", dataset_name, ".txt")
  )
  setkey(all_data, Nmut, id1, id2)
  
  #restrict to training set (+ downsampling)
  set.seed(1603)
  all_data = rbind(all_data[Nmut <= 1],
      all_data[Nmut == 2 & test_set == F][sample(.N,round(.N * dataset_size))])
  
  dataset_size_label = gsub("\\.", "p", dataset_size)

  # which variants are present?
  id_mut <- unique(c(
    all_data[, id1]
  ))
  id_idx <- which(id_mut %in% id_mut)
  
  idlist <- list()
  idlist$id_len <- length(id_mut)
  idlist$id_mut <- id_mut
  idlist$id_idx <- id_idx
  
  # key variants
  all_data[, id1_key := which(id_mut == id1), id1]
  all_data[, id2_key := which(id_mut == id2), id2]
  all_data[test_set == F, rank := 1:.N]
  
  #variant to mutant lookup matrix
  varXmut <- Matrix(0, nrow = all_data[test_set == F, .N], ncol = length(id_mut), sparse = T)
  varXmut[unlist(rbind(all_data[test_set == F, cbind(rank, key = id1_key)], 
                       all_data[!is.na(id2) & test_set==F, cbind(rank, key = id2_key)]))] <- 1
  
  
  set.seed(iteration)
  
  if (approach == 1) {
    #fitness and error values
    fw_list <- list(
      b_fitness = all_data[test_set == F, b_fitness],
      b_sigma = all_data[test_set == F, b_sigma]
    )
  } else if (approach == 2) {
   ## use in vitro s_ddg values from Nishtal 2019 to constrain s_ddg values
    structural_data = fread("processed_data/GB1_structural_data.txt")

    all_data[Nmut == 1, s_ddg_measured := 
      structural_data[id == id1, s_ddg_measured],id1]
    all_data[Nmut == 1, s_ddg_measured_sdsmooth := 
      structural_data[id == id1, s_ddg_measured_sdsmooth],id1]

    #fitness and error values + in vitro s_ddg values
    fw_list <- list(
      b_fitness = all_data[test_set == F, b_fitness],
      b_sigma = all_data[test_set == F, b_sigma],
      s_ddg_invitro = all_data[Nmut <= 1][idlist$id_idx, s_ddg_measured],
      s_ddg_invitro_sd = all_data[Nmut <= 1][idlist$id_idx, s_ddg_measured_sdsmooth]
    )    
  }

  # upper/lower bounds for parameters
  lower_bounds <- c(
    rep(-15, idlist$id_len *2),
    c(-15, -1, -8), # b_f0 (third parameter here) is on log-scale
    -15
  )
  upper_bounds <- c(
    rep(15, idlist$id_len * 2),
    c(15, 1, -3), # b_f0 (third parameter here) is on log-scale
    15
  )

  # fixed parameters
  fixed_par <- c(0, 0)
  names(fixed_par) <- c("0X_s_ddg", "0X_b_ddg")
  
  if (task == 1) {
    ###### calculate initial set of models (100)

    tic() 
    # sample start parameters
    if (approach == 1) {
      start_par <- c(
        rep(0, idlist$id_len * 2), # s_ddg and b_ddg values + wildtype
        rnorm(3, mean = c(0, 0, -5.5), sd = c(1, 0.05, 1)), # binding dgwt, fwt and f0
        rnorm(1, mean = c(0), sd = c(1)) # stability dgwt
      )
    } else if (approach == 2) {
      start_par <- c(
        fw_list$s_ddg_invitro,
        rep(0, idlist$id_len), # s_ddg and b_ddg values + wildtype
        rnorm(3, mean = c(0, 0, -5.5), sd = c(1, 0.05, 1)), # binding dgwt, fwt and f0
        rnorm(1, mean = c(0), sd = c(1)) # stability dgwt
      )
      start_par[is.na(start_par)] = 0
    }
    
    names(start_par) <- c(
      paste0(idlist$id_mut, "_s_ddg"),
      paste0(idlist$id_mut, "_b_ddg"),
      c("b_dgwt", "b_fwt", "b_f0"), # b_fwt & b_f0 are on log-scale
      c("s_dgwt")
    )
    
    # fix parameters and enforce bounds
    for (i in seq_along(fixed_par)) {
      start_par[names(start_par) == names(fixed_par)[i]] <- fixed_par[i]
    }
    start_par[start_par < lower_bounds] = lower_bounds[start_par < lower_bounds]
    start_par[start_par > upper_bounds] = upper_bounds[start_par > upper_bounds]
    
    global_model <- optim(
      par = start_par,
      fn = function_GB1_dG_fitting,
      gr = function_GB1_dG_gradient,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = maxit),
      idlist = idlist,
      fw_list = fw_list,
      varXmut = varXmut,
      predict_binding0 = predict_binding0,
      fixed_par = fixed_par
    )

    toc()
    
    # convert to dt
    x <- data.table(t(global_model$par), global_model$value, global_model$convergence)
    names(x) <- c(names(global_model$par), "objective", "convergence")
    
    filename = paste0("processed_data/dG/", 
                            dataset_name, "_approach", approach,
                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
                            "_size", dataset_size_label, 
                            "_initialmodels.txt")
    if (file.exists(filename)) {
      write.table(x, 
                file = filename,
                quote = F,
                row.names = F,
                col.names = F,
                append = T
      )
    } else {
      write.table(x, 
                file = filename,
                quote = F,
                row.names = F
      )
    }

    models_dt <- fread(filename)
    if (nrow(models_dt) >= n_models) {


      ########## boostrap best model to improve ############
      ####### by varying best parameters randomly
      models_dt <- fread(paste0("processed_data/dG/", 
                              dataset_name, "_approach", approach,
                              ifelse(predict_binding0 == TRUE, "_b0", ""),
                              "_size", dataset_size_label, 
                              "_initialmodels.txt"))
      best_model <- models_dt[which.min(objective)]
      start_par <- best_model[, 
                              unlist(.SD), 
                              .SDcols = !(names(models_dt) %in% c("objective", "convergence"))]
      
      while (nrow(best_model) < 10) {
        
        start_par0 <- start_par + 
          rnorm(length(start_par), mean = 0, sd = 0.01 * abs(start_par))
        # fix parameters and enforce bounds
        for (i in seq_along(fixed_par)) {
          start_par0[names(start_par0) == names(fixed_par)[i]] <- fixed_par[i]
        }
        start_par0[start_par0 < lower_bounds] = lower_bounds[start_par0 < lower_bounds]
        start_par0[start_par0 > upper_bounds] = upper_bounds[start_par0 > upper_bounds]
        
        
        boot_model <- optim(
          par = start_par0,
          fn = function_GB1_dG_fitting,
          gr = function_GB1_dG_gradient,
          method = "L-BFGS-B",
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(maxit = maxit),
          idlist = idlist,
          fw_list = fw_list,
          varXmut = varXmut,
          predict_binding0 = predict_binding0,
          fixed_par = fixed_par
        )
        
        boot_model_dt <- data.table(
          t(boot_model$par), 
          boot_model$value, 
          boot_model$convergence
        )
        names(boot_model_dt) <- c(
          names(boot_model$par),
          "objective",
          "convergence"
        )
        obj_diff <- (boot_model_dt$objective / best_model[1, objective]) - 1
        if (obj_diff < -1e-6) {
          best_model <- boot_model_dt
          start_par <- start_par <- best_model[, 
                       unlist(.SD), 
                       .SDcols = !(names(models_dt) %in% c("objective", "convergence"))]
        } else {
          best_model <- rbind(
            best_model,
            boot_model_dt
          )
        }
        
        print(paste0("n = ", nrow(best_model), ", diff:", round(obj_diff, 6)))
      }
      write.table(best_model, 
                  file = paste0("processed_data/dG/", 
                                dataset_name, "_approach", approach,
                                ifelse(predict_binding0 == TRUE, "_b0", ""),
                                "_size", dataset_size_label,
                                "_boot_startpar.txt"),
                  quote = F,
                  row.names = F
      )

    }
  } else if (task == 2) {
    #######################################################
    ########## parameter uncertainty estimates ############
    #######################################################
  

    ############ bootstrap models by varying fitness values 
    boot_par_models <- fread(paste0("processed_data/dG/", 
                                    dataset_name, "_approach", approach,
                                    ifelse(predict_binding0 == TRUE, "_b0", ""),
                                    "_size", dataset_size_label, 
                                    "_boot_startpar.txt"))
    best_model <- boot_par_models[which.min(objective)]
    start_par <- best_model[, 
                                 unlist(.SD), 
                                 .SDcols = !(names(boot_par_models) %in% c("objective", "convergence"))]
    
    fw_boot <- fw_list
    fw_boot$b_fitness <- fw_boot$b_fitness + 
      rnorm(n = length(fw_boot$b_fitness), mean = 0, sd = fw_boot$b_sigma)
    
    boot_model <- optim(
      par = start_par,
      fn = function_GB1_dG_fitting,
      gr = function_GB1_dG_gradient,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = maxit),
      idlist = idlist,
      fw_list = fw_boot,
      varXmut = varXmut,
      predict_binding0 = predict_binding0,
      fixed_par = fixed_par
    )
    
    x <- data.table(t(c(boot_model$par, boot_model$value, boot_model$convergence)))
    names(x) <- c(names(boot_model$par), "objective", "convergence")

    
    filename = paste0("processed_data/dG/", 
                            dataset_name, "_approach", approach,
                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
                            "_size", dataset_size_label,
                            "_boot_fitness.txt")

    if (file.exists(filename)) {
      write.table(x, 
                file = filename,
                quote = F,
                row.names = F,
                col.names = F,
                append = T
      )
    } else {
      write.table(x, 
                file = filename,
                quote = F,
                row.names = F
      )
    }

  } else if (task == 3) {
    #######################################################
    ########## parameter uncertainty estimates ############
    #######################################################
    ########## perform profile likelihood estimates #######

    boot_par_models <- fread(paste0("processed_data/dG/", 
                                    dataset_name, "_approach", approach,
                                    ifelse(predict_binding0 == TRUE, "_b0", ""),
                                    "_size", dataset_size_label, 
                                    "_boot_startpar.txt"))
    best_model <- boot_par_models[which.min(objective)]
    start_par <- best_model[, 
                                 unlist(.SD), 
                                 .SDcols = !(names(boot_par_models) %in% c("objective", "convergence"))]

    pll_par <- data.table(
      par = names(start_par), 
      value = start_par, 
      lower = 0, 
      upper = 0
    )
    alpha <- 0.05
    df <- nrow(pll_par)
    rel_stepsize <- 0.05
    
    for (df_idx in df) {
      pll_bounds <- foreach(
        p = 1:nrow(pll_par),
        .combine = rbind
      ) %dopar% {


        record <- TRUE
        # step size for variation
        if (grepl("dd[gG]", pll_par[p, par])) { 
          #it's a ddg value
          ddg_type <- pll_par[p, paste0(strsplit(par, "_")[[1]][2:3], 
                                        collapse = "_")]
          range <- pll_par[grep(ddg_type, par), quantile(value, c(0.01, 0.99))]
          step_size <- abs(diff(range)) * rel_stepsize
          record <- FALSE
        } else if (grepl("s_dgwt", pll_par[p, par])) { 
          #it's a stability dgwt value
          range <- pll_par[grep("s_ddg", par), quantile(value, c(0.01, 0.99))]
          step_size <- abs(diff(range)) * rel_stepsize
        } else if (grepl("b_dgwt", pll_par[p, par])) { 
          #it's binding dgwt value
          range <- pll_par[grep("b_ddg", par), quantile(value, c(0.01, 0.99))]
          step_size <- abs(diff(range)) * rel_stepsize
        } else { 
          #it's a fitness value
          step_size <- rel_stepsize
        }
        
        tic()
        pll_models <- data.table(t(start_par))
        if (!(pll_par[p, par] %in% names(fixed_par))) { 
          for (dir in c(-1, 1)) {
            non_sig <- TRUE
            i <- 1
            e <- 0
            start_par0 <- start_par
            
            while (non_sig & i < 20 & between(pll_par[p, value],lower_bounds[p],upper_bounds[p])) {
              fixed_par0 <- c(fixed_par, pll_par[p, value] + step_size * dir * i)
              names(fixed_par0)[length(fixed_par0)] <- pll_par[p, par]
              
              for (idx in seq_along(fixed_par0)) {
                start_par0[names(start_par0) == names(fixed_par0)[idx]] <- 
                  fixed_par0[idx]
              }
              
              pll_model <- optim(
                par = start_par0,
                fn = function_GB1_dG_fitting,
                gr = function_GB1_dG_gradient,
                method = "L-BFGS-B",
                lower = lower_bounds,
                upper = upper_bounds,
                control = list(maxit = maxit),
                idlist = idlist,
                fw_list = fw_list,
                varXmut = varXmut,
                predict_binding0 = predict_binding0,
                fixed_par = fixed_par0
              )
              
              #evaluate objective function diff to best model
              obj_diff <- (pll_model$value - best_model[1, objective])
              if (obj_diff > qchisq(1 - alpha, df)) {
                if (e >= -2) {
                  e <- e - 1
                  i <- i - 2^e
                } else {
                  non_sig <- FALSE
                }
              } else {
                i <- i + 2^e
                if (record == TRUE) {
                  pll_models <- rbind(
                    pll_models, 
                    data.table(t(pll_model$par))
                  )
                }
                
              }
              start_par0 <- pll_model$par
              
              #save results
              if (dir == -1 & (non_sig == FALSE | i >= 20)) {
                pll_par[p, lower := value + step_size * dir * (i - 2^e)]
              } else if (dir == 1 & (non_sig == FALSE | i >= 20)) {
                pll_par[p, upper := value + step_size * dir * (i - 2^e)]
              }
            }
            
          }
          if (record == TRUE) {
            write.table(pll_models, 
                        file = paste0("processed_data/dG/", 
                                      dataset_name,"_approach", approach,
                                      ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                      "_size", dataset_size_label,
                                      "_parameter_pll_df", df, "_", pll_par[p, par], ".txt"),
                        quote = F,
                        row.names = F
            )
          }
          print(paste0("parameter ", p, ", df = ", df))
          toc()
          return(pll_par[p, ])
        }
      }



      write.table(pll_bounds, 
                  file = paste0("processed_data/dG/", 
                                dataset_name,"_approach", approach,
                                ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                "_size", dataset_size_label,
                                "_parameter_pll_df", df, ".txt"),
                  quote = F,
                  row.names = F
      )
    }
  }
}
    

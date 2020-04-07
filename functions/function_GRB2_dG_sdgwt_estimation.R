function_GRB2_dG_sdgwt_estimation <- function(
  dataset_name = "GRB2_dG_dataset",
  predict_binding0 = FALSE,
  method,
  maxit = 1e6,
  n_models = 100,
  n_cores = 15
) { 

### allows for s_dgwt to be different between stability and binding PCA assays

  require(data.table)
  require(Matrix)
  require(tictoc)
  require(foreach)
  require(doMC)
  registerDoMC(cores = n_cores)
  
  all_data <- fread(
    paste0("processed_data/", dataset_name, ".txt")
  )
  setkey(all_data, Nmut, id1, id2)
  
  #restrict variant set
  if (method == 1) {
    all_data = all_data[Nmut  <= 1 & !is.na(b_fitness)] #model 1: only use single mutants
  } else if (method == 2) {
    all_data = all_data[test_set == F]
  } else if (method == 3) {
    all_data = all_data[test_set == F & !is.na(b_fitness)]
    all_data[, s_fitness := NA] #model 3: delete all stability fitness values
  } else if (method == 4) {
    all_data = all_data[test_set == F & Nmut %in% c(0,2)] #model 4: only use double mutants
  }
  
  # which variants are present?
  id_mut <- unique(c(
    all_data[!is.na(s_fitness) | !is.na(b_fitness), id1], 
    all_data[(!is.na(s_fitness) | !is.na(b_fitness)) & !is.na(id2), id2]
  ))
  
  s_id_mut <- unique(c(
    all_data[!is.na(s_fitness) | !is.na(b_fitness), id1], 
    all_data[(!is.na(s_fitness) | !is.na(b_fitness)) & !is.na(id2), id2]
  ))
  s_id_idx <- which(id_mut %in% s_id_mut)
  
  b_id_mut <- unique(c(
    all_data[!is.na(b_fitness), id1], 
    all_data[!is.na(b_fitness) & !is.na(id2), id2]
  ))
  b_id_idx <- which(id_mut %in% b_id_mut)
  
  idlist <- list()
  idlist$id_len <- length(id_mut)
  idlist$s_id_mut <- s_id_mut
  idlist$s_id_idx <- s_id_idx
  idlist$s_id_len <- length(s_id_idx)
  idlist$b_id_mut <- b_id_mut
  idlist$b_id_idx <- b_id_idx
  idlist$b_id_len <- length(b_id_idx)
  
  # key variants
  all_data[, id1_key := which(id_mut == id1), id1]
  all_data[, id2_key := which(id_mut == id2), id2]
  all_data[, rank := 1:.N]
  
  #variant to mutant lookup matrix
  varXmut <- Matrix(0, nrow = all_data[, .N], ncol = length(id_mut), sparse = T)
  varXmut[unlist(rbind(all_data[, cbind(rank, key = id1_key)], 
                       all_data[!is.na(id2), cbind(rank, key = id2_key)]))] <- 1
  
  #fitness and error values
  fw_list <- list(
    s_fitness = all_data[, s_fitness],
    s_sigma = all_data[, s_sigma],
    b_fitness = all_data[, b_fitness],
    b_sigma = all_data[, b_sigma]
  )
  
  # upper/lower bounds for parameters
  lower_bounds <- c(
    rep(-15, idlist$s_id_len + idlist$b_id_len),
    c(-15, 0.7, -0.2),
    c(-15, -15, 0.7, -0.2)
  )
  upper_bounds <- c(
    rep(15, idlist$s_id_len + idlist$b_id_len),
    c(15, 1.3, 0.6),
    c(15, 15, 1.3, 0.6)
  )
  if (method == 3) {
    lower_bounds = lower_bounds[1 : (length(lower_bounds)-2)]
    upper_bounds = upper_bounds[1 : (length(upper_bounds)-2)]
  }

  # fixed parameters
  fixed_par <- c(0, 0)#,
                  # 0.15, 0.275)
  names(fixed_par) <- c("0X_s_ddg", "0X_b_ddg")#,
                          # "s_f0","b_f0")
  ### calculate initial set of models
   
  models <- foreach(m = 1:n_models) %dopar% {
    tic() 
    # sample start parameters
    start_par <- c(
      rep(0, idlist$s_id_len + idlist$b_id_len), # s_ddg and b_ddg values + wildtype
      rnorm(3, mean = c(0, 1, 0.275), sd = c(1, 0.05, 0.05)), # binding dgwt, fwt and f0
      rnorm(4, mean = c(0, 0, 1, 0.15), sd = c(1, 1, 0.05, 0.05)) # stability dgwt, fwt and f0
    )
    names(start_par) <- c(
      paste0(idlist$s_id_mut, "_s_ddg"),
      paste0(idlist$b_id_mut, "_b_ddg"),
      c("b_dgwt", "b_fwt", "b_f0"),
      c("sb_dgwt", "ss_dgwt", "s_fwt", "s_f0")
    )
    if (method == 3) {
      start_par = start_par[1 : (length(start_par)-2)]
    }
    
    # fix parameters and enforce bounds
    for (i in seq_along(fixed_par)) {
      start_par[names(start_par) == names(fixed_par)[i]] <- fixed_par[i]
    }
    start_par[start_par < lower_bounds] = lower_bounds[start_par < lower_bounds]
    start_par[start_par > upper_bounds] = upper_bounds[start_par > upper_bounds]
    
    global_model <- optim(
      par = start_par,
      fn = function_GRB2_dG_sdgwt_fitting,
      gr = function_GRB2_dG_sdgwt_gradient,
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
    return(global_model)
    print(m)
    toc()
  }
  
  
  models_dt <- data.table(t(sapply(X = seq_along(models), FUN = function(X) {
    x <- c(models[[X]]$par, models[[X]]$value, models[[X]]$convergence)
    names(x) <- c(names(models[[X]]$par), "objective", "convergence")
    return(x)
  })))
  write.table(models_dt, 
              file = paste0("processed_data/dG/", 
                            dataset_name, "_sdgwt_method", method, 
                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
                            "_initialmodels.txt"),
              quote = F,
              row.names = F
  )
  
  ########## boostrap best model to improve ############
  ####### by varying best parameters randomly
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
      fn = function_GRB2_dG_sdgwt_fitting,
      gr = function_GRB2_dG_sdgwt_gradient,
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
                            dataset_name, "_sdgwt_method", method, 
                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
                            "_boot_startpar.txt"),
              quote = F,
              row.names = F
  )
  
  
  #######################################################
  ########## parameter uncertainty estimates ############
  #######################################################
  
  boot_par_models <- fread(paste0("processed_data/dG/", 
                                  dataset_name, "_sdgwt_method", method, 
                                  ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                  "_boot_startpar.txt"))
  best_model <- boot_par_models[which.min(objective)]
  start_par <- best_model[, 
                               unlist(.SD), 
                               .SDcols = !(names(boot_par_models) %in% c("objective", "convergence"))]
  
  ############ bootstrap models by varying fitness values 
  boot_models <- foreach(m = 1:n_models) %dopar% {
    fw_boot <- fw_list
    if (method != 3) {
      fw_boot$s_fitness <- fw_boot$s_fitness + 
        rnorm(n = length(fw_boot$s_fitness), mean = 0, sd = fw_boot$s_sigma)  
    }
    fw_boot$b_fitness <- fw_boot$b_fitness + 
      rnorm(n = length(fw_boot$b_fitness), mean = 0, sd = fw_boot$b_sigma)
    
    boot_model <- optim(
      par = start_par,
      fn = function_GRB2_dG_sdgwt_fitting,
      gr = function_GRB2_dG_sdgwt_gradient,
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
    
    return(boot_model)
  }
  boot_models_dt <- data.table(t(sapply(X = seq_along(boot_models), FUN = function(X) {
    x <- c(boot_models[[X]]$par, boot_models[[X]]$value, boot_models[[X]]$convergence)
    names(x) <- c(names(boot_models[[X]]$par), "objective", "convergence")
    return(x)
  })))
  
  write.table(boot_models_dt, 
              file = paste0("processed_data/dG/", 
                            dataset_name, "_sdgwt_method", method,
                            ifelse(predict_binding0 == TRUE, "_b0", ""), 
                            "_boot_fitness.txt"),
              quote = F,
              row.names = F
  )
  
  ############################################################
  ########## perform profile likelihood estimates ############
  
  pll_par <- data.table(
    par = names(start_par), 
    value = start_par, 
    lower = 0, 
    upper = 0
  )
  alpha <- 0.05
  df <- c(1 , nrow(pll_par))
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
      } else if (grepl("s[sb]_dgwt", pll_par[p, par])) { 
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
              fn = function_GRB2_dG_sdgwt_fitting,
              gr = function_GRB2_dG_sdgwt_gradient,
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
            if (obj_diff > qchisq(1 - alpha, df_idx)) {
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
                                    dataset_name, "_sdgwt_method", method,
                                    ifelse(predict_binding0 == TRUE, "_b0", ""), 
                                    "_parameter_pll_df", df_idx, "_", pll_par[p, par], ".txt"),
                      quote = F,
                      row.names = F
          )
        }
        print(paste0("parameter ", p, ", df = ", df_idx))
        toc()
        return(pll_par[p, ])
      }
    }
    write.table(pll_bounds, 
                file = paste0("processed_data/dG/", 
                              dataset_name, "_sdgwt_method", method,
                              ifelse(predict_binding0 == TRUE, "_b0", ""), 
                              "_parameter_pll_df", df_idx, ".txt"),
                quote = F,
                row.names = F
    )
  }
}
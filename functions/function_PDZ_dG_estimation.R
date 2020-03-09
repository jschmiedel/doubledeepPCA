function_PDZ_dG_estimation <- function(
  dataset_name = "PDZ_dG_dataset",
  dataset_suffix,
  maxit = 1e6,
  n_models = 100,
  n_cores = 15
) {

  require(data.table)
  require(tictoc)
  require(foreach)
  require(doMC)
  registerDoMC(cores = n_cores)

  # read dataset
  all_data <- fread(
    paste0("processed_data/", dataset_name, ".txt")
  )
  setkey(all_data, id1)
  
  
  # which variants are present?
  id_var <- unique(c(
    all_data[!is.na(s_fitness) | !is.na(b1_fitness) | !is.na(b2_fitness) | !is.na(b3_fitness), id1]
  ))

  # how many variants?
  idlist <- list()
  idlist$s <- length(id_var)
  idlist$b12 <- length(id_var)
  idlist$b3 <- length(id_var)
  
  #list of fitness and sigma values
  fw_list <- list(
    s_fitness = all_data[, s_fitness],
    s_sigma = all_data[, s_sigma],
    b1_fitness = all_data[, b1_fitness],
    b1_sigma = all_data[, b1_sigma],
    b2_fitness = all_data[, b2_fitness],
    b2_sigma = all_data[, b2_sigma],
    b3_fitness = all_data[, b3_fitness],
    b3_sigma = all_data[, b3_sigma]
  )

  #upper/lower bounds for parameters
  lower <- c(
    rep(-15, idlist$s + idlist$b12 + idlist$b3 + 5),
    # rep(-15, idlist$s + idlist$b12 + 5),
    rep(c(0.8, 0), 4)
  )
  upper <- c(
    rep(15, idlist$s + idlist$b12 + idlist$b3 + 5),
    # rep(15, idlist$s + idlist$b12 + 5),
    rep(c(1.2, 0.6), 4)
  )

  fixed_par <- c(
    0, 0, 0,
    0.05,
    0.1,
    0.25,
    0.425
  )
  names(fixed_par) <- c(
    "0X_s_ddg", "0X_b12_ddg", "0X_b3_ddg",
    "s_f0",
    "b1_f0",
    "b2_f0",
    "b3_f0"
  )

  ### calculate a set of models
  tic()  
  models <- foreach(m = 1:n_models) %dopar% {
    
    start_par <- c(
      rep(0, idlist$s), # s_ddg values
      rep(0, idlist$b12), # b_ddg values to CRIPT
      rep(0, idlist$b3), # b_ddg values to mutated peptide
      rnorm(5, 0), # wildtype dG values (s1, s23, b1, b2, b3),
      rnorm(1, 0.95, 0.05), # s_fwt parameter
      rnorm(1, 0.05, 0.05), # s_f0 parameter
      rnorm(1, 0.95, 0.05), # b1_fwt parameter
      rnorm(1, 0.1, 0.05), # b1_f0 parameter
      rnorm(1, 1.025, 0.05), # b2_fwt parameter
      rnorm(1, 0.25, 0.05), # b2_f0 parameter
      rnorm(1, 0.95, 0.05), # b3_fwt parameter
      rnorm(1, 0.425, 0.05) # b3_f0 parameter
    )
    names(start_par) <- c(
      paste0(id_var, "_s_ddg"),
      paste0(id_var, "_b12_ddg"),
      paste0(id_var, "_b3_ddg"),
      "s1_dgwt",
      "s23_dgwt",
      "b1_dgwt",
      "b2_dgwt",
      "b3_dgwt",
      "s_fwt",
      "s_f0",
      "b1_fwt",
      "b1_f0",
      "b2_fwt",
      "b2_f0",
      "b3_fwt",
      "b3_f0"
    )

    for (i in seq_along(fixed_par)) {
        start_par[names(start_par) == names(fixed_par)[i]] <- fixed_par[i]
    }
    
    global_model <- optim(
      par = start_par,
      fn = function_PDZ_dG_fitting2,
      gr = function_PDZ_dG_gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit),
      idlist = idlist,
      fw_list = fw_list,
      fixed_par = fixed_par
    )
    
   return(global_model)
  }
  toc()

  models_dt <- data.table(t(sapply(X = seq_along(models), FUN = function(X) {
    x <- c(models[[X]]$par, models[[X]]$value, models[[X]]$convergence)
    names(x) <- c(names(models[[X]]$par), "objective", "convergence")
    return(x)
  })))
  write.table(models_dt, 
              file = paste0("processed_data/", 
                dataset_name, "_initialmodels.txt"),
              quote = F,
              row.names = F
  )
    
  ########## boostrap best model to improve ############
  
  ####### varying best parameters randomly
  best_model <- models_dt[which.min(objective)]
  start_par <- best_model[, 
      unlist(.SD), 
      .SDcols = !(names(models_dt) %in% c("objective", "convergence"))]

  while (nrow(best_model) < 100) {
    
    start_par0 <- start_par + 
      rnorm(length(start_par), mean = 0, sd = 0.01 * abs(start_par))
    for (i in seq_along(fixed_par)) {
        start_par0[names(start_par0) == names(fixed_par)[i]] <- fixed_par[i]
    }

    boot_model <- optim(
      par = start_par0,
      fn = function_PDZ_dG_fitting2,
      gr = function_PDZ_dG_gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit),
      idlist = idlist,
      fw_list = fw_list,
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
              file = paste0("processed_data/", 
                dataset_name, "_boot_startpar.txt"),
              quote = F,
              row.names = F
  )

  ############ bootstrapping fitness values
  best_model <- best_model[which.min(objective)]
  start_par <- best_model[, 
      unlist(.SD), 
      .SDcols = !(names(best_model) %in% c("objective", "convergence"))]

  while (nrow(best_model) < 100) {
    fw_boot <- fw_list
    fw_boot$s_fitness <- fw_boot$s_fitness + 
      rnorm(n = length(fw_boot$s_fitness), mean = 0, sd = fw_boot$s_sigma)
    fw_boot$b1_fitness <- fw_boot$b1_fitness + 
      rnorm(n = length(fw_boot$b1_fitness), mean = 0, sd = fw_boot$b1_sigma)
    fw_boot$b2_fitness <- fw_boot$b2_fitness + 
      rnorm(n = length(fw_boot$b2_fitness), mean = 0, sd = fw_boot$b2_sigma)
    fw_boot$b3_fitness <- fw_boot$b3_fitness + 
      rnorm(n = length(fw_boot$b3_fitness), mean = 0, sd = fw_boot$b3_sigma)

    boot_model <- optim(
      par = start_par0,
      fn = function_PDZ_dG_fitting2,
      gr = function_PDZ_dG_gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit),
      idlist = idlist,
      fw_list = fw_boot,
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
              file = paste0("processed_data/", 
                dataset_name, "_boot_fitness.txt"),
              quote = F,
              row.names = F
  )

  

  ############################################################
  ########## perform profile likelihood estimates ############
  ############################################################

  pll_par <- data.table(
    par = names(start_par), 
    value = start_par, 
    lower = 0, 
    upper = 0
  )
  alpha <- 0.05
  df <- 1#nrow(pll_par)
  rel_stepsize <- 0.05

  pll_bounds <- foreach(
    # p = 1:nrow(pll_par),
    p = 4498:4510,
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
    } else if (grepl("s[0-9]+_dgwt", pll_par[p, par])) { 
    #it's a stability dgwt value
      range <- pll_par[grep("s_ddg", par), quantile(value, c(0.01, 0.99))]
      step_size <- abs(diff(range)) * rel_stepsize
    } else if (grepl("b[12]+_dgwt", pll_par[p, par])) { 
    #it's binding 1 or 2 dgwt value
      range <- pll_par[grep("b12_ddg", par), quantile(value, c(0.01, 0.99))]
      step_size <- abs(diff(range)) * rel_stepsize
    } else if (grepl("b3_dgwt", pll_par[p, par])) { 
    #it's a binding 3 dgwt value
      range <- pll_par[grep("b3_dd[gG]", par), quantile(value, c(0.01, 0.99))]
      step_size <- abs(diff(range)) * rel_stepsize
    } else { 
    #it's a fitness value
      step_size <- rel_stepsize
    }

    tic()
    pll_models <- data.table(t(start_par))
    for (dir in c(-1, 1)) {
      non_sig <- TRUE
      i <- 1
      e <- 0
      start_par0 <- start_par
    
      while (non_sig & i < 20) {
        fixed_par0 <- c(fixed_par, pll_par[p, value] + step_size * dir * i)
        names(fixed_par0)[length(fixed_par0)] <- pll_par[p, par]
        
        for (idx in seq_along(fixed_par0)) {
          start_par0[names(start_par0) == names(fixed_par0)[idx]] <- 
            fixed_par0[idx]
        }
        
        pll_model <- optim(
          par = start_par0,
          fn = function_PDZ_dG_fitting2,
          gr = function_PDZ_dG_gradient,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = list(maxit = maxit),
          idlist = idlist,
          fw_list = fw_list,
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
              file = paste0("processed_data/tmp/", 
                dataset_name, "_parameter_pll_", pll_par[p, par], ".txt"),
              quote = F,
              row.names = F
      )
    }
    print(p)
    toc()
    return(pll_par[p, ])
  }
  write.table(pll_bounds, 
              file = paste0("processed_data/", 
                dataset_name, "_parameter_pll_point.txt"),
              quote = F,
              row.names = F
  )
  
}

function_dG_estimation <- function(
  dataset_name,
  dataset_suffix,
  method,
  maxit = 100,
  execute = TRUE
) {
  if (!execute) {
    return()
  }
  require(data.table)
  require(tictoc)
  
  print(paste0("method: ", method))
  
  # read dataset
  all_data <- fread(
    paste0("processed_data/", dataset_name, ".txt")
  )
  setkey(all_data, Nmut, id1, id2)
  
  
  if (method == 1) {
    all_data = all_data[Nmut == 1]
  } else if (method == 3) {
    all_data[, s_fitness := NA]
  } else if (method == 4) {
    all_data = all_data[Nmut == 2]
  }
  # which variants are present?
  s_id_var <- unique(c(
    all_data[!is.na(s_fitness) | !is.na(b_fitness), id1], 
    all_data[(!is.na(s_fitness) | !is.na(b_fitness)) & !is.na(id2), id2]
  ))
  
  b_id_var <- unique(c(
    all_data[!is.na(b_fitness), id1], 
    all_data[!is.na(b_fitness) & !is.na(id2), id2]
  ))
  
  # how many variants?
  s_idl <- length(s_id_var)
  b_idl <- length(b_id_var)
  
  # key variants
  all_data[, s_id1_key := which(s_id_var == id1), id1]
  all_data[, s_id2_key := which(s_id_var == id2), id2]
  all_data[, b_id1_key := which(b_id_var == id1), id1]
  all_data[, b_id2_key := which(b_id_var == id2), id2]
  
  all_data_long = rbind(
    all_data[!is.na(s_fitness) & !is.na(s_sigma), 
             .(id1_key = s_id1_key,
               id2_key = s_id2_key,
               fitness = s_fitness,
               sigma = s_sigma,
               type = "stability"
             )],
    all_data[!is.na(b_fitness) & !is.na(b_sigma), 
             .(id1_key = b_id1_key,
               id2_key = b_id2_key,
               fitness = b_fitness,
               sigma = b_sigma,
               type = "binding"
             )]
  )
  
  i <- 1
  best_model <- c()
  while (max(table(best_model)) <= 10) {
    
    print(paste0("model: ", i))
    tic()
    
    # create fitness & error vectors using bootstrapped fitness values
    if (max(table(best_model)) < 10) {
      all_data_iteration <- copy(all_data_long[, .(
        id1_key,
        id2_key,
        type,
        fitness = fitness + rnorm(.N, mean = 0, sd = sigma),
        sigma
      )])
    } else { #last iteration with best estimates of fitness
       all_data_iteration <- copy(all_data_long[, .(
        id1_key,
        id2_key,
        type,
        fitness,
        sigma
      )]) 
    }
    
    #upper/lower bounds for parameters
    lower <- c(
      rep(-15, s_idl + b_idl + 2),
      rep(0, 4)
    )
    upper <- c(
      rep(15, s_idl + b_idl + 2),
      rep(1.2, 4)
    )

    #starting parameters for model fit
    if (i == 1) {
      if (method == 1) {
        
        # initialize global parameters
        start_par <- c(
          rep(0, s_idl), # s_ddg values
          rnorm(2, 0), # wildtype dG values (first s then b),
          rnorm(1, 1, 0.1), # b_fwt parameter
          rnorm(1, 0, 0.1), # b_f0 parameter
          rnorm(1, 1, 0.1), # s_fwt parameter
          rnorm(1, 0, 0.1) # s_f0 parameter
        )

        #fit global relationship assuming b_ddg = 0
        global_model <- optim(
          par = start_par,
          fn = function_dG_fitting,
          method = "L-BFGS-B",
          lower = lower[c(1:s_idl, s_idl + 1:4)],
          upper = upper[c(1:s_idl, s_idl + 1:4)],
          control = list(maxit = 500),
          s_idl = s_idl,
          b_idl = b_idl,
          dt = copy(all_data_iteration),
          m1_mod = 0
        )
      
        #update starting parameter
        start_par = c(
          global_model$par[1 : s_idl],
          rep(0, b_idl), #b_ddg values
          global_model$par[(s_idl + 1) : (s_idl + 6)]
        )
        
      }
      else if (method %in% c(2,4)) {
        start_par <- c(
          rep(0, s_idl), # s_ddg values
          rep(0, b_idl), # b_ddg values
          rnorm(2, 0), # wildtype dG values,
          rnorm(1, 1, 0.1), # b_fwt parameter
          rnorm(1, 0, 0.1), # b_f0 parameter
          rnorm(1, 1, 0.1), # s_fwt parameter
          rnorm(1, 0, 0.1) # s_f0 parameter
        )
      } else if (method == 3) {
        start_par <- c(
          rep(0, s_idl), # s_ddg values
          rep(0, b_idl), # b_ddg values
          rnorm(2, 0), # wildtype dG values,
          rnorm(1, 1, 0.1), # b_fwt parameter
          rnorm(1, 0, 0.1) # b_f0 parameter
        )
        lower <- lower[1 : (s_idl + b_idl + 4)]
        upper <- upper[1 : (s_idl + b_idl + 4)]
      }
    }
    
    global_model <- optim(
      par = start_par,
      fn = function_dG_fitting,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit),
      s_idl = s_idl, 
      b_idl = b_idl, 
      dt = copy(all_data_iteration)
    )
    
    model_i <- data.table(global_model$value,t(global_model$par))  
    if (method != 3) {
      names(model_i) <- c(
        "objective",
        paste0(s_id_var,"_s_ddG"),
        paste0(b_id_var,"_b_ddG"),
        "s_dgwt",
        "b_dgwt",
        "b_fwt",
        "b_f0",
        "s_fwt",
        "s_f0"
      )
    } else {
      names(model_i) <- c(
        "objective",
        paste0(s_id_var,"_s_ddG"),
        paste0(b_id_var,"_b_ddG"),
        "s_dgwt",
        "b_dgwt",
        "b_fwt",
        "b_f0"
      )
    }
    
    # add model to model data.table
    if (i == 1) {
      dO <- 0
      models <- copy(model_i)
    } else {
      dO <- log10(abs(model_i$objective - min(models$objective)) /
        min(models$objective))
      models <- rbind(models,model_i)
    }
    
    if (max(table(best_model)) < 10) {
      best_model[i] <- models[, which.min(objective)]
    } else {
      best_model[i] <- best_model[i - 1]
    }
    start_par <- models[which.min(objective), 
      unlist(.SD), 
      .SDcols = !grepl("objective",names(models))]
    
    print(paste0("best model: " ,best_model[i], ", n = ", sum(best_model == best_model[i]), ", log10(rel_diff) = ", round(dO,2)))
    
    write.table(models[, 
                  cbind(objective,round(.SD,3)),
                  .SDcols = !grepl("objective",names(models))], 
                file = paste0("processed_data/tmp/method", method, "_", dataset_name, "_", dataset_suffix, ".txt"),
                quote = F,
                row.names = F
    )
    
    toc()
    i = i + 1
    
  }  
}

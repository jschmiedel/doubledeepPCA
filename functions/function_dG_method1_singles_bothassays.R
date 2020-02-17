######## dG estimation method 1: calculate free energies of folding and binding from stabilityPCA and deepPCA data of single mutants ########
#### version: 3.0
#### created: 2019/11/27 by J??rn
#### last modified: 2020/01/15 by J??rn
# v3: parallelized v2

function_dG_method1_singles_bothassays <- function(
                                     dataset_name,
                                     dataset_suffix,
                                     n_cores,
                                     n_bootstraps,
                                     execute = TRUE
) {

  if (!execute) {
    return()
  }

  require(data.table)
  require(foreach)
  require(doMC)
  registerDoMC(cores = n_cores)

  # read dataset
  all_data <- 
    fread(paste0("processed_data/", dataset_name, "_dG_dataset", dataset_suffix, ".txt"))
  setkey(all_data, Nmut, id1)

  # which/how many variants?
  id_L <- all_data[Nmut == 1, .N]
  id_var <- all_data[Nmut == 1, id]
  
  for (i in 1:n_bootstraps) {
    print(i)

    if (i == 1) { #calculate 50 models with starting parameters
      n_par <- 50
    } else { #then optimize best model
      n_par <- n_cores
    } 

    models0 <- foreach(m = 1:n_par) %dopar% {

      #starting parameters for model fit
      if (i == 1) {
        start_par <- c(
          rep(0, id_L * 2), #s_dG & b_dG values
          1e-5 + (1 - 2e-5) * runif(2), # scale parameters
          0.1 * runif(2), # bgr parameters
          rep(0, 2) # wildtype dG values
        )
      }

      # bootstrap fitness data
      list_fs <- list(
        s_fitness = all_data[Nmut == 1, 
          s_fitness + rnorm(.N, mean = 0, sd = s_sigma)],
        s_sigma = all_data[Nmut == 1, s_sigma],
        b_fitness = all_data[Nmut == 1, 
          b_fitness + rnorm(.N, mean = 0, sd = b_sigma)],
        b_sigma = all_data[Nmut == 1, b_sigma]
      )

      ####### first, fit global relationship by assuming binding dG are all zero
      fit_global_relationship <- optim(
        par = start_par[c(1:id_L, (id_L * 2 + 1):(id_L * 2 + 6))],
        fn = function_dG_method1_fitting,
        method = "L-BFGS-B",
        lower = c(
          rep(-10, id_L),
          rep(1e-5, 2),
          rep(0, 2),
          rep(-10, 2)
        ),
        upper = c(
          rep(10, id_L),
          rep(1 - 1e-5, 2),
          rep(1, 2),
          rep(10, 2)
        ),
        control = list(maxit = 500),
        id_L = id_L, 
        global_par = global_par, 
        list_fs = list_fs
      )

      ####### second, estimate binding dGs (and again stability dG)
      fit_dGs <- optim(
        par = c(
          fit_global_relationship$par[1:id_L],
          start_par[(id_L + 1):(2 * id_L)]
        ),
        fn = function_dG_method1_fitting,
        method = "L-BFGS-B",
        lower = c(rep(-10, 2 * id_L)),
        upper = c(rep(10, 2 * id_L)),
        control = list(maxit = 500),
        id_L = id_L, 
        global_par = fit_global_relationship$par[(id_L + 1):(id_L + 6)], 
        list_fs = list_fs
      )

      # gather model parameters
      fit_dGs$par <- c(
        fit_dGs$par,
        fit_global_relationship$par[(id_L + 1):(id_L + 6)]
      )
      names(fit_dGs$par) <- c(
        id_var,
        id_var,
        "s_scale",
        "b_scale",
        "s_bgr",
        "b_bgr",
        "s_dGwt",
        "b_dGwt"
      )

      return(fit_dGs)
    }

    if (i == 1) {
        models <- models0
    } else {
        models <- c(models, models0)
    }

    parameters <- t(sapply(X = seq_along(models), FUN = function(X) {
        models[[X]]$par
    }))
    objective <- sapply(X = seq_along(models), FUN = function(X) {
        models[[X]]$value
    })

    start_par <- parameters[which.min(objective), ]
  }


  # save
  save(models, 
    file = paste0("processed_data/dG_method1_",dataset_name, dataset_suffix, "_",      
      n_bootstraps * n_cores, "models_", n_bootstraps, "iterations.Rdata"
    ))

  # plot models
  function_dG_plot_models(
    models_name = paste0("dG_method1_", dataset_name, dataset_suffix, "_",
      n_bootstraps * n_cores, "models_", n_bootstraps, "iterations"),
    dataset_name = dataset_name,
  )
}


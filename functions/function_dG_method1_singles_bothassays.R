######## dG estimation method 1: calculate free energies of folding and binding from stabilityPCA and deepPCA data of single mutants ########
#### version: 3.0
#### created: 2019/11/27 by J??rn
#### last modified: 2020/01/15 by J??rn
# v3: parallelized v2

function_dG_method1_singles_bothassays <- function(
                                                   dataset_name,
                                                   n_cores = 15,
                                                   n_bootstraps = 100,
                                                   bgr_set = c(), # c(s_bgr,b_bgr)
                                                   execute = TRUE) {
  if (!execute) {
    return()
  }

  require(data.table)
  require(foreach)
  require(doMC)

  registerDoMC(cores = n_cores)

  # read dataset
  all_data <- fread(paste0("processed_data/", dataset_name, "_dG_dataset.txt"))
  setkey(all_data, Nmut, id1)

  # use only singles
  list_fs <- list(
    s_fitness = all_data[Nmut == 1, s_fitness],
    s_sigma = all_data[Nmut == 1, s_sigma],
    b_fitness = all_data[Nmut == 1, b_fitness],
    b_sigma = all_data[Nmut == 1, b_sigma]
  )

  # which/how many variants?
  id_L <- all_data[Nmut == 1, .N]
  id_var <- all_data[Nmut == 1, id]

  models <- foreach(m = 1:n_bootstraps) %dopar% {
    ####### first, fit global relationship by assuming binding dG are all zero
    if (is.null(bgr_set)) {
      start_par <- c(
        # rep(0, id_L), #s_dGs
        1e-5 + (1 - 2e-5) * runif(2), #scale parameters
        0.1 * runif(2), #bgr parameters
        rep(1, 2) #wildtype fitness values
      )
      global_par <- c()
    } else {
      start_par <- c(
        # rep(0, id_L), #s_dGs
        1e-5 + (1 - 2e-5) * runif(2), #scale parameters
        rep(1, 2) #wildtype fitness values
      )
      global_par <- bgr_set
    }

    fit_global_relationship <- optim( ###### WHY DOESN"T THIS WORK ANYMORE? GIVES ERROR 52
      par = start_par,
      fn = function_dG_method1_fitting,
      method = "L-BFGS-B",
      lower = c(
        # rep(-10, id_L), 
        rep(1e-5, 2), 
        rep(0, ifelse(is.null(bgr_set), 2, 0)), 
        rep(0.8, 2) #adjust this to apparent wild-type fitness values
      ),
      upper = c(
        # rep(10, id_L), 
        rep(1 - 1e-5, 2), 
        rep(1, ifelse(is.null(bgr_set), 2, 0)),
        rep(1.1, 2) #adjust this to apparent wild-type fitness values
      ),
      control = list(maxit = 5000),
      id_L = id_L, global_par = global_par, list_fs = list_fs
    )

    # global_par <- c(fit_global_relationship$par[(id_L + 1):
    #   (id_L + ifelse(is.null(bgr_set), 6, 4))], bgr_set)
    global_par <- c(fit_global_relationship$par[1:ifelse(is.null(bgr_set), 6, 4)],
      bgr_set)

    ####### second, estimate binding dGs (and again stability dG)
    fit_dGs <- optim(
      par = c(
        fit_global_relationship$par[1:id_L],
        rep(0, id_L)
      ),
      fn = function_dG_method1_fitting,
      method = "L-BFGS-B",
      lower = c(rep(-10, 2 * id_L)),
      upper = c(rep(10, 2 * id_L)),
      control = list(maxit = 5000),
      id_L = id_L, global_par = global_par, list_fs = list_fs
    )

    # gather model parameters
    fit_dGs$par <- c(
      fit_dGs$par,
      global_par,
      function_folding_F2dG(
        s_fitness = global_par[5], 
        s_bgr = global_par[3], 
        s_scale = global_par[1]
      ),
      function_binding_F2dG(
        b_fitness = global_par[6], 
        s_dG = function_folding_F2dG(
          s_fitness = global_par[5], 
          s_bgr = global_par[3], 
          s_scale = global_par[1]
        ),
        b_bgr = global_par[4], 
        b_scale = global_par[2]
      )
    )
    names(fit_dGs$par) <- c(
      id_var, 
      id_var, 
      "s_scale", 
      "b_scale", 
      "s_bgr", 
      "b_bgr", 
      "s_Ft", 
      "b_Fwt",
      "s_dGwt",
      "b_dGwt"
    )

    return(fit_dGs)
  }
  # save
  save(models, file = paste0("processed_data/dG_method1_", dataset_name, "_", n_bootstraps, "models.Rdata"))

  # plot models
  function_dG_plot_models(
    models_name = paste0("dG_method1_", dataset_name, "_", n_bootstraps, "models"),
    dataset_name = dataset_name
  )
}


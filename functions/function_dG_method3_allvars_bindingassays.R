function_dG_method3_allvars_bindingassay <- function(
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

    # this model uses only binding assay data
    # it doesn't see stability fitness data and doesn't fit scale and 
    # background parameters for those
    # the wild-type stability dG value becomes a fit parameter

    # read dataset
    all_data <- fread(paste0("processed_data/", 
        dataset_name, "_dG_dataset", dataset_suffix, ".txt"))
    setkey(all_data, Nmut, id1, id2)

    # which variants are present?
    id_var <- unique(all_data[!is.na(s_fitness) & !is.na(b_fitness), id1])
    all_data[, id1_key := which(id_var == id1), id1]
    all_data[, id2_key := which(id_var == id2), id2]
    id_L <- length(id_var)

    # grab singles and doubles
    singles <- all_data[Nmut == 1,
        .(id1_key, s_fitness, b_fitness, s_sigma, b_sigma)]
    doubles <- all_data[Nmut == 2 , 
        .(id1_key, id2_key, s_fitness, b_fitness, s_sigma, b_sigma)]
    # key them
    list_keys <- list(
        singles_key = singles[, id1_key],
        doubles_key1 = doubles[, id1_key],
        doubles_key2 = doubles[, id2_key]
    )
    
    for (i in seq_len(n_bootstraps)) {
        print(i)

        if (i == 1) { #calculate 50 models with starting parameters
          n_par <- 50
        } else { #then optimize best model
          n_par <- n_cores
        } 

        # fit parameters
        models0 <- foreach(m = 1:n_par) %dopar% {
            
            #starting parameters for model fit
            if (i == 1) {
                start_par <- c(
                    rep(0, id_L),
                    rep(0, id_L),
                    1e-5 + (1 - 2e-5) * runif(1),
                    0.1 * runif(1),
                    rep(0, 2)
                )
            }

            # create fitness & error vectors
            list_fs <- list(
                b_fitness = c(
                    singles[, b_fitness + rnorm(.N, mean = 0, sd = b_sigma)], 
                    doubles[, b_fitness + rnorm(.N, mean = 0, sd = b_sigma)]),
                b_sigma = c(
                    singles[, b_sigma], 
                    doubles[, b_sigma])
            )
            global_model <- optim(
                par = start_par,
                fn = function_dG_method3_fitting,
                method = "L-BFGS-B",
                lower = c(rep(-10, 2 * id_L), 1e-5, 0, rep(-10, 2)),
                upper = c(rep(10, 2 * id_L), 1 - 1e-5, 1, rep(10, 2)),
                control = list(maxit = 500),
                id_L = id_L, 
                list_fs = list_fs, 
                list_keys = list_keys
            )
            names(global_model$par) <- c(
                id_var, 
                id_var, 
                "b_scale", 
                "b_bgr", 
                "s_dGwt", 
                "b_dGwt")

            global_model$iteration <- i

            return(global_model)
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
    save(
        models, 
        file = 
            paste0("processed_data/dG_method3_", dataset_name, dataset_suffix, "_", 
                n_bootstraps * n_cores, "models_", n_bootstraps, "iterations.Rdata")
    )

    # plot models
    function_dG_plot_models(
        models_name = 
            paste0("dG_method3_", dataset_name, dataset_suffix, "_", 
            n_bootstraps * n_cores, "models_", n_bootstraps, "iterations"),
        dataset_name = dataset_name
    ) 
}

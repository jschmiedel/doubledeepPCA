function_GRB2_dG_gradient <- function(
  parameters,
  fixed_par,
  predict_binding0 = FALSE,
  idlist,
  fw_list,
  varXmut
 ) {

    for (i in seq_along(fixed_par)) {
        parameters[names(parameters) == names(fixed_par)[i]] <- fixed_par[i]
    }

    s_ddg_mut <- rep(NA, idlist$id_len)
    s_ddg_mut[idlist$s_id_idx] <- parameters[1:idlist$s_id_len]
    s_ddg_var <- varXmut %*% s_ddg_mut

    b_ddg_mut <- rep(NA, idlist$id_len)
    b_ddg_mut[idlist$b_id_idx] <- parameters[(idlist$s_id_len + 1) : 
                                            (idlist$s_id_len + idlist$b_id_len)]
    b_ddg_var <- varXmut %*% b_ddg_mut
   
    all_idl <- idlist$s_id_len + idlist$b_id_len

    mutXvar <- t(varXmut)

    ## stability phenotype
    if (sum(!is.na(fw_list$s_fitness)) > 0) { #is zero for method 3
        gradient_s <- function_stability_dG2F_gradient(
          s_ddg = s_ddg_var, 
          s_dgwt = parameters[all_idl + 4], 
          s_fwt = parameters[all_idl + 5],
          s_f0 = parameters[all_idl + 6],  
          f = fw_list$s_fitness,
          w = fw_list$s_sigma,
          mutXvar = mutXvar
        )
    }

    ## binding phenotype
    gradient_b <- function_binding_dG2F_gradient( 
        s_ddg = s_ddg_var,
        b_ddg = b_ddg_var,
        b_dgwt = parameters[all_idl + 1],
        b_fwt = parameters[all_idl + 2],
        b_f0 = parameters[all_idl + 3], 
        s_dgwt = parameters[all_idl + 4],
        f = fw_list$b_fitness,
        w = fw_list$b_sigma,
        mutXvar = mutXvar
    )

    #binding phenotype if ddg = 0
    if (predict_binding0 == TRUE) {
        gradient_b0 <- function_binding_dG2F_gradient(
            s_ddg = s_ddg_var,
            b_ddg = 0,   
            b_dgwt = parameters[all_idl + 1],
            b_fwt = parameters[all_idl + 2],
            b_f0 = parameters[all_idl + 3], 
            s_dgwt = parameters[all_idl + 4],
            f = fw_list$b_fitness,
            w = fw_list$b_sigma,
            mutXvar = mutXvar
        )
    }

    #list and sum gradients from individual contributions
    if (sum(!is.na(fw_list$s_fitness)) > 0) {
        if (predict_binding0 == TRUE) {
            gradient <- c(
                gradient_s$s_ddg[idlist$s_id_idx] + 
                    gradient_b$s_ddg[idlist$s_id_idx] + 
                    gradient_b0$s_ddg[idlist$s_id_idx], #s_dgg
                gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
                gradient_b$b_dgwt + gradient_b0$b_dgwt, #b_dgwt
                gradient_b$b_fwt + gradient_b0$b_fwt, #b_fwt
                gradient_b$b_f0 + gradient_b0$b_f0, #b_f0
                gradient_s$s_dgwt + 
                    gradient_b$s_dgwt + gradient_b0$s_dgwt, #s_dgwt
                gradient_s$s_fwt, #s_fwt
                gradient_s$s_f0 #s_f0
            )
        } else {
            gradient <- c(
                gradient_s$s_ddg[idlist$s_id_idx] + 
                    gradient_b$s_ddg[idlist$s_id_idx], #s_dgg
                gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
                gradient_b$b_dgwt, #b_dgwt
                gradient_b$b_fwt, #b_fwt
                gradient_b$b_f0, #b_f0
                gradient_s$s_dgwt + gradient_b$s_dgwt, #s_dgwt
                gradient_s$s_fwt, #s_fwt
                gradient_s$s_f0 #s_f0
            )
        }
    } else {
        if (predict_binding0 == TRUE) {
            gradient <- c(
                gradient_b$s_ddg[idlist$s_id_idx] + 
                    gradient_b0$s_ddg[idlist$s_id_idx], #s_dgg
                gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
                gradient_b$b_dgwt + gradient_b0$b_dgwt, #b_dgwt
                gradient_b$b_fwt + gradient_b0$b_fwt, #b_fwt
                gradient_b$b_f0 + gradient_b0$b_f0, #b_f0
                gradient_b$s_dgwt + gradient_b0$s_dgwt #s_dgwt
            )
        } else {
            gradient <- c(
                gradient_b$s_ddg[idlist$s_id_idx], #s_dgg
                gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
                gradient_b$b_dgwt, #b_dgwt
                gradient_b$b_fwt, #b_fwt
                gradient_b$b_f0, #b_f0
                gradient_b$s_dgwt #s_dgwt
            )
        }
    }
    names(gradient) <- names(parameters)
    
    #set gradients from fixed parameters to 0
    for (i in seq_along(fixed_par)) {
        gradient[names(gradient) == names(fixed_par)[i]] <- 0
    }
    return(gradient)
}

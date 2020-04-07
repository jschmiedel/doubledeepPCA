function_GRB2_dG_4state_gradient <- function(
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

    s1_ddg_mut <- rep(NA, idlist$id_len)
    s1_ddg_mut[idlist$s_id_idx] <- parameters[1 : idlist$s_id_len]
    s1_ddg_var <- varXmut %*% s1_ddg_mut

    s2_ddg_mut <- rep(NA, idlist$id_len)
    s2_ddg_mut[idlist$s_id_idx] <- parameters[(idlist$s_id_len + 1) : (idlist$s_id_len * 2)]
    s2_ddg_var <- varXmut %*% s2_ddg_mut

    b_ddg_mut <- rep(NA, idlist$id_len)
    b_ddg_mut[idlist$b_id_idx] <- parameters[(idlist$s_id_len * 2 + 1) : 
                                          (idlist$s_id_len * 2 + idlist$b_id_len)]
    b_ddg_var <- varXmut %*% b_ddg_mut

    all_idl <- idlist$s_id_len * 2 + idlist$b_id_len

    mutXvar <- t(varXmut)

    ## stability phenotype
    if (sum(!is.na(fw_list$s_fitness)) > 0) { #is zero for method 3
        gradient_s <- function_dG_4state_folding_dg2fitness_gradient(
          s1_ddg = s1_ddg_var, 
          s2_ddg = s2_ddg_var, 
          s1_dgwt = parameters[all_idl + 4], 
          s2_dgwt = parameters[all_idl + 5], 
          s_fwt = parameters[all_idl + 6],
          s_f0 = parameters[all_idl + 7],  
          f = fw_list$s_fitness,
          w = fw_list$s_sigma,
          mutXvar = mutXvar
        )
    }

    ## binding phenotype
    gradient_b <- function_dG_4state_binding_dg2fitness_gradient( 
        b_ddg = b_ddg_var,
        s1_ddg = s1_ddg_var,
        s2_ddg = s2_ddg_var,
        b_dgwt = parameters[all_idl + 1],
        s1_dgwt = parameters[all_idl + 4],
        s2_dgwt = parameters[all_idl + 5],
        b_fwt = parameters[all_idl + 2],
        b_f0 = parameters[all_idl + 3], 
        f = fw_list$b_fitness,
        w = fw_list$b_sigma,
        mutXvar = mutXvar
    )

    #binding phenotype if ddg = 0
    if (predict_binding0 == TRUE) {
        gradient_b0 <- function_dG_4state_binding_dg2fitness_gradient( 
          b_ddg = 0,
          s1_ddg = s1_ddg_var,
          s2_ddg = s2_ddg_var,
          b_dgwt = parameters[all_idl + 1],
          s1_dgwt = parameters[all_idl + 4],
          s2_dgwt = parameters[all_idl + 5],
          b_fwt = parameters[all_idl + 2],
          b_f0 = parameters[all_idl + 3], 
          f = fw_list$b_fitness,
          w = fw_list$b_sigma,
          mutXvar = mutXvar
        )
    }

    #list and sum gradients from individual contributions
    
    if (predict_binding0 == TRUE) {
        gradient <- c(
            gradient_s$s1_ddg[idlist$s_id_idx] + 
                gradient_b$s1_ddg[idlist$s_id_idx] + 
                gradient_b0$s1_ddg[idlist$s_id_idx], #s1_dgg
            gradient_s$s2_ddg[idlist$s_id_idx] + 
                gradient_b$s2_ddg[idlist$s_id_idx] + 
                gradient_b0$s2_ddg[idlist$s_id_idx], #s2_dgg
            gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
            gradient_b$b_dgwt + gradient_b0$b_dgwt, #b_dgwt
            gradient_b$b_fwt + gradient_b0$b_fwt, #b_fwt
            gradient_b$b_f0 + gradient_b0$b_f0, #b_f0
            gradient_s$s1_dgwt + 
                gradient_b$s1_dgwt + gradient_b0$s1_dgwt, #s1_dgwt
            gradient_s$s2_dgwt + 
                gradient_b$s2_dgwt + gradient_b0$s2_dgwt, #s2_dgwt
            gradient_s$s_fwt, #s_fwt
            gradient_s$s_f0 #s_f0
        )
    } else {
        gradient <- c(
            gradient_s$s1_ddg[idlist$s_id_idx] + 
                gradient_b$s1_ddg[idlist$s_id_idx], #s1_dgg
            gradient_s$s2_ddg[idlist$s_id_idx] + 
                gradient_b$s2_ddg[idlist$s_id_idx], #s2_dgg
            gradient_b$b_ddg[idlist$b_id_idx], #b_dgg
            gradient_b$b_dgwt, #b_dgwt
            gradient_b$b_fwt, #b_fwt
            gradient_b$b_f0, #b_f0
            gradient_s$s1_dgwt + gradient_b$s1_dgwt, #s1_dgwt
            gradient_s$s2_dgwt + gradient_b$s2_dgwt, #s2_dgwt
            gradient_s$s_fwt, #s_fwt
            gradient_s$s_f0 #s_f0
        )
    }
    
    names(gradient) <- names(parameters)
    
    #set gradients from fixed parameters to 0
    for (i in seq_along(fixed_par)) {
        gradient[names(gradient) == names(fixed_par)[i]] <- 0
    }
    return(gradient)
}

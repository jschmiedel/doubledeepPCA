function_GB1_dG_gradient <- function(
  parameters,
  fixed_par,
  predict_binding0 = FALSE,
  approach = 1,
  idlist,
  fw_list,
  varXmut
 ) {

    for (i in seq_along(fixed_par)) {
        parameters[names(parameters) == names(fixed_par)[i]] <- fixed_par[i]
    }

    s_ddg_mut <- rep(NA, idlist$id_len)
    s_ddg_mut[idlist$id_idx] <- parameters[1 : idlist$id_len]
    s_ddg_var <- varXmut %*% s_ddg_mut

    b_ddg_mut <- rep(NA, idlist$id_len)
    b_ddg_mut[idlist$id_idx] <- parameters[(idlist$id_len + 1) : 
                                            (idlist$id_len * 2)]
    b_ddg_var <- varXmut %*% b_ddg_mut
   
    all_idl <- idlist$id_len * 2

    mutXvar <- t(varXmut)

    ## binding phenotype
    gradient_b <- function_binding_dg2logf_gradient( 
        s_ddg = s_ddg_var,
        b_ddg = b_ddg_var,
        b_dgwt = parameters[all_idl + 1],
        b_fwt = parameters[all_idl + 2],
        b_f0 = parameters[all_idl + 3], 
        s_dgwt = parameters[all_idl + 4],
        f = fw_list$b_fitness,
        w = fw_list$b_sigma,
        mutXvar = mutXvar,
        logF = TRUE
    )

    #binding phenotype if ddg = 0
    if (predict_binding0 == TRUE) {
        gradient_b0 <- function_binding_dg2logf_gradient(
            s_ddg = s_ddg_var,
            b_ddg = 0,   
            b_dgwt = parameters[all_idl + 1],
            b_fwt = parameters[all_idl + 2],
            b_f0 = parameters[all_idl + 3], 
            s_dgwt = parameters[all_idl + 4],
            f = fw_list$b_fitness,
            w = fw_list$b_sigma,
            mutXvar = mutXvar,
            logF = TRUE
        )
    }

    #list and sum gradients from individual contributions
    if (predict_binding0 == TRUE & approach == 1) {
        gradient <- c(
            gradient_b$s_ddg[idlist$id_idx] + 
                gradient_b0$s_ddg[idlist$id_idx], #s_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt + gradient_b0$b_dgwt, #b_dgwt
            gradient_b$b_fwt + gradient_b0$b_fwt, #b_fwt
            gradient_b$b_f0 + gradient_b0$b_f0, #b_f0
            gradient_b$s_dgwt + gradient_b0$s_dgwt #s_dgwt
        )
    } else if (predict_binding0 == FALSE & approach == 1) {
        gradient <- c(
            gradient_b$s_ddg[idlist$id_idx], #s_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt, #b_dgwt
            gradient_b$b_fwt, #b_fwt
            gradient_b$b_f0, #b_f0
            gradient_b$s_dgwt #s_dgwt
        )
    } else if (predict_binding0 == FALSE & approach == 2) {
        gradient_sddg_invitro <- - 2 * (fw_list$s_ddg_invitro - s_ddg_mut) / fw_list$s_ddg_invitro_sd^2
        gradient_sddg_invitro[s_ddg_mut > 4 & fw_list$s_ddg_invitro == 4] = 0 
        gradient <- c(
            rowSums(cbind(gradient_b$s_ddg[idlist$id_idx],gradient_sddg_invitro), na.rm = T), #s_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt, #b_dgwt
            gradient_b$b_fwt, #b_fwt
            gradient_b$b_f0, #b_f0
            gradient_b$s_dgwt #s_dgwt
        )
    }
    
    names(gradient) <- names(parameters)
    
    #set gradients from fixed parameters to 0
    for (i in seq_along(fixed_par)) {
        gradient[names(gradient) == names(fixed_par)[i]] <- 0
    }
    return(gradient)
}

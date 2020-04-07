function_GB1_dG_4state_gradient <- function(
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

    s1_ddg_mut <- rep(NA, idlist$id_len)
    s1_ddg_mut[idlist$id_idx] <- parameters[1:idlist$id_len]
    s1_ddg_var <- varXmut %*% s1_ddg_mut

    s2_ddg_mut <- rep(NA, idlist$id_len)
    s2_ddg_mut[idlist$id_idx] <- parameters[(idlist$id_len + 1): (idlist$id_len * 2)]
    s2_ddg_var <- varXmut %*% s2_ddg_mut

    b_ddg_mut <- rep(NA, idlist$id_len)
    b_ddg_mut[idlist$id_idx] <- parameters[(idlist$id_len * 2 + 1) : 
                                        (idlist$id_len * 3)]
    b_ddg_var <- varXmut %*% b_ddg_mut

    all_idl <- idlist$id_len * 3

    mutXvar <- t(varXmut)

    ## binding phenotype
    gradient_b <- function_dG_4state_binding_dg2logf_gradient( 
        s1_ddg = s1_ddg_var,
        s2_ddg = s2_ddg_var,
        b_ddg = b_ddg_var,
        b_dgwt = parameters[all_idl + 1],
        b_fwt = parameters[all_idl + 2],
        b_f0 = parameters[all_idl + 3], 
        s1_dgwt = parameters[all_idl + 4],
        s2_dgwt = parameters[all_idl + 5],
        f = fw_list$b_fitness,
        w = fw_list$b_sigma,
        mutXvar = mutXvar
    )

    #binding phenotype if ddg = 0
    if (predict_binding0 == TRUE) {
        gradient_b0 <- function_dG_4state_binding_dg2logf_gradient(
            s1_ddg = s1_ddg_var,
            s2_ddg = s2_ddg_var,
            b_ddg = 0,
            b_dgwt = parameters[all_idl + 1],
            b_fwt = parameters[all_idl + 2],
            b_f0 = parameters[all_idl + 3], 
            s1_dgwt = parameters[all_idl + 4],
            s2_dgwt = parameters[all_idl + 5],
            f = fw_list$b_fitness,
            w = fw_list$b_sigma,
            mutXvar = mutXvar
        )
    }

    #list and sum gradients from individual contributions
    if (predict_binding0 == TRUE & approach == 1) {
        gradient <- c(
            gradient_b$s1_ddg[idlist$id_idx] + 
                gradient_b0$s1_ddg[idlist$id_idx], #s1_dgg
            gradient_b$s2_ddg[idlist$id_idx] + 
                gradient_b0$s2_ddg[idlist$id_idx], #s2_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt + gradient_b0$b_dgwt, #b_dgwt
            gradient_b$b_fwt + gradient_b0$b_fwt, #b_fwt
            gradient_b$b_f0 + gradient_b0$b_f0, #b_f0
            gradient_b$s1_dgwt + gradient_b0$s1_dgwt, #s1_dgwt
            gradient_b$s2_dgwt + gradient_b0$s2_dgwt #s2_dgwt
        )
    } else if (predict_binding0 == FALSE & approach == 1) {
        gradient <- c(
            gradient_b$s1_ddg[idlist$id_idx], #s1_dgg
            gradient_b$s2_ddg[idlist$id_idx], #s2_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt, #b_dgwt
            gradient_b$b_fwt, #b_fwt
            gradient_b$b_f0, #b_f0
            gradient_b$s1_dgwt, #s1_dgwt
            gradient_b$s2_dgwt #s2_dgwt
        )
    } else if (predict_binding0 == FALSE & approach == 2) {

        #gradient for comparision to in vitro folding ddgs
        #assuming 
        rt = 1.99e-3 * 310.15
        e1wt <- exp(-parameters[all_idl + 4]/rt)
        e2wt <- exp(-parameters[all_idl + 5]/rt)
        e1iwt <- e1wt * exp(-s1_ddg_mut / rt)
        e2iwt <- e2wt * exp(-s2_ddg_mut / rt)
        s_ddg_wt = -log(e1wt + e2wt) * rt
        s_ddg_mut = -log(e1iwt + e2iwt) * rt - s_ddg_wt

        gradient_sddg_invitro <- - 2 * (fw_list$s_ddg_invitro - s_ddg_mut) / fw_list$s_ddg_invitro_sd^2
        gradient_sddg_invitro[s_ddg_mut > 4 & fw_list$s_ddg_invitro == 4] = 0 

        gradient_sddg_invitro_s1_ddg <- gradient_sddg_invitro * e1iwt / (e1iwt + e2iwt)
        gradient_sddg_invitro_s2_ddg <- gradient_sddg_invitro * e2iwt / (e1iwt + e2iwt)
        gradient_sddg_invitro_s1_dgwt <- sum(gradient_sddg_invitro * (e1iwt / (e1iwt + e2iwt) + e1wt / (e1wt + e2wt)), na.rm = T)
        gradient_sddg_invitro_s2_dgwt <- sum(gradient_sddg_invitro * (e2iwt / (e1iwt + e2iwt) + e2wt / (e1wt + e2wt)), na.rm = T)

        gradient <- c(
            rowSums(cbind(gradient_b$s1_ddg[idlist$id_idx],gradient_sddg_invitro_s1_ddg), na.rm = T), #s1_dgg
            rowSums(cbind(gradient_b$s2_ddg[idlist$id_idx],gradient_sddg_invitro_s2_ddg), na.rm = T), #s2_dgg
            gradient_b$b_ddg[idlist$id_idx], #b_dgg
            gradient_b$b_dgwt, #b_dgwt
            gradient_b$b_fwt, #b_fwt
            gradient_b$b_f0, #b_f0
            gradient_b$s1_dgwt + gradient_sddg_invitro_s1_dgwt, #s1_dgwt
            gradient_b$s2_dgwt + gradient_sddg_invitro_s2_dgwt #s2_dgwt
        )
    }
    
    names(gradient) <- names(parameters)
    
    #set gradients from fixed parameters to 0
    for (i in seq_along(fixed_par)) {
        gradient[names(gradient) == names(fixed_par)[i]] <- 0
    }
    return(gradient)
}

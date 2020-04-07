function_PDZ_dG_gradient <- function(
  parameters,
  fixed_par, 
  predict_binding0,
  idlist,
  fw_list,
  varXmut
 ) {

    for (i in seq_along(fixed_par)) {
        parameters[names(parameters) == names(fixed_par)[i]] <- fixed_par[i]
    }


    s_ddg_mut <- rep(NA, idlist$id_len)
    s_ddg_mut[idlist$id_idx] <- parameters[1:idlist$id_len]
    s_ddg_var <- varXmut %*% s_ddg_mut

    b12_ddg_mut <- rep(NA, idlist$id_len)
    b12_ddg_mut[idlist$id_idx] <- parameters[(idlist$id_len + 1):(idlist$id_len * 2)]
    b12_ddg_var <- varXmut %*% b12_ddg_mut

    b3_ddg_mut <- rep(NA, idlist$id_len)
    b3_ddg_mut[idlist$id_idx] <- parameters[(idlist$id_len * 2 + 1):(idlist$id_len * 3)]
    b3_ddg_var <- varXmut %*% b3_ddg_mut

    all_idl <- idlist$id_len * 3

    mutXvar <- t(varXmut)

    ## stability phenotype
    gradient_s <- function_stability_dG2F_gradient(
      s_ddg = s_ddg_var, 
      s_dgwt = parameters[all_idl + 1], 
      s_fwt = parameters[all_idl + 6],
      s_f0 = parameters[all_idl + 7],  
      f = fw_list$s_fitness,
      w = fw_list$s_sigma,
      mutXvar = mutXvar
    ) 

    ## binding1 phenotype (ddPCA measurement)
    gradient_b1 <- function_binding_dG2F_gradient(
        s_ddg = s_ddg_var,
        b_ddg = b12_ddg_var,
        s_dgwt = parameters[all_idl + 1],   
        b_dgwt = parameters[all_idl + 3],
        b_fwt = parameters[all_idl + 8],
        b_f0 = parameters[all_idl + 9], 
        f = fw_list$b1_fitness,
        w = fw_list$b1_sigma,
        mutXvar = mutXvar
    )
    #if b12_ddG are 0
    if (predict_binding0 == TRUE) {
        gradient_b1_0 <- function_binding_dG2F_gradient(
            s_ddg = s_ddg_var,
            b_ddg = 0,
            s_dgwt = parameters[all_idl + 1],   
            b_dgwt = parameters[all_idl + 3],
            b_fwt = parameters[all_idl + 8],
            b_f0 = parameters[all_idl + 9], 
            f = fw_list$b1_fitness,
            w = fw_list$b1_sigma,
            mutXvar = mutXvar
        )
    }
    
    ## binding2 phenotype (McLaughlin CRIPT measurement)
    gradient_b2 <- function_binding_dG2F_gradient(
        s_ddg = s_ddg_var,
        b_ddg = b12_ddg_var,
        s_dgwt = parameters[all_idl + 2],   
        b_dgwt = parameters[all_idl + 4],
        b_fwt = parameters[all_idl + 10],
        b_f0 = parameters[all_idl + 11], 
        f = fw_list$b2_fitness,
        w = fw_list$b2_sigma,
        mutXvar = mutXvar
    )
    #if b12_ddG are 0
    if (predict_binding0 == TRUE) {
        gradient_b2_0 <- function_binding_dG2F_gradient(
            s_ddg = s_ddg_var,
            b_ddg = 0,
            s_dgwt = parameters[all_idl + 2],   
            b_dgwt = parameters[all_idl + 4],
            b_fwt = parameters[all_idl + 10],
            b_f0 = parameters[all_idl + 11], 
            f = fw_list$b2_fitness,
            w = fw_list$b2_sigma,
            mutXvar = mutXvar
        )
    }

    ## binding3 phenotype (McLaughlin Tm2F measurement)
    gradient_b3 <- function_binding_dG2F_gradient(
        s_ddg = s_ddg_var,
        b_ddg = b3_ddg_var,
        s_dgwt = parameters[all_idl + 2],   
        b_dgwt = parameters[all_idl + 5],
        b_fwt = parameters[all_idl + 12],
        b_f0 = parameters[all_idl + 13], 
        f = fw_list$b3_fitness,
        w = fw_list$b3_sigma,
        mutXvar = mutXvar
    )
    #if only b3_ddG are 0
    if (predict_binding0 == TRUE) {
        gradient_b3_0 <- function_binding_dG2F_gradient(
            s_ddg = s_ddg_var,
            b_ddg = 0,
            s_dgwt = parameters[all_idl + 2],   
            b_dgwt = parameters[all_idl + 5],
            b_fwt = parameters[all_idl + 12],
            b_f0 = parameters[all_idl + 13], 
            f = fw_list$b3_fitness,
            w = fw_list$b3_sigma,
            mutXvar = mutXvar
        )
    }

    #list and sum gradients from individual contributions
    if (predict_binding0 == TRUE) {
        gradient <- c(
            gradient_s$s_ddg[idlist$id_idx] + 
                gradient_b1$s_ddg[idlist$id_idx] + gradient_b1_0$s_ddg[idlist$id_idx] + 
                gradient_b2$s_ddg[idlist$id_idx] + gradient_b2_0$s_ddg[idlist$id_idx] +
                gradient_b3$s_ddg[idlist$id_idx] + gradient_b3_0$s_ddg[idlist$id_idx], #s_dgg
            gradient_b1$b_ddg[idlist$id_idx] + gradient_b2$b_ddg[idlist$id_idx], #b12_dgg
            gradient_b3$b_ddg[idlist$id_idx], #b3_ddg
            gradient_s$s_dgwt + gradient_b1$s_dgwt + gradient_b1_0$s_dgwt, #s1_dgwt
            gradient_b2$s_dgwt + gradient_b2_0$s_dgwt + 
                gradient_b3$s_dgwt + gradient_b3_0$s_dgwt, #s23_dgwt
            gradient_b1$b_dgwt + gradient_b1_0$b_dgwt, #b1_dgwt
            gradient_b2$b_dgwt + gradient_b2_0$b_dgwt, #b2_dgwt
            gradient_b3$b_dgwt + gradient_b3_0$b_dgwt, #b3_dgwt
            gradient_s$s_fwt, #s_fwt
            gradient_s$s_f0, #s_f0
            gradient_b1$b_fwt + gradient_b1_0$b_fwt, #b1_fwt
            gradient_b1$b_f0 + gradient_b1_0$b_f0, #b1_f0
            gradient_b2$b_fwt + gradient_b2_0$b_fwt, #b2_fwt
            gradient_b2$b_f0 + gradient_b2_0$b_f0, #b2_f0
            gradient_b3$b_fwt + gradient_b3_0$b_fwt, #b3_fwt
            gradient_b3$b_f0 + gradient_b3_0$b_f0 #b3_f0
        )
    } else {
        gradient <- c(
            gradient_s$s_ddg[idlist$id_idx] + 
                gradient_b1$s_ddg[idlist$id_idx] +  
                gradient_b2$s_ddg[idlist$id_idx] + 
                gradient_b3$s_ddg[idlist$id_idx], #s_dgg
            gradient_b1$b_ddg[idlist$id_idx] + gradient_b2$b_ddg[idlist$id_idx], #b12_dgg
            gradient_b3$b_ddg[idlist$id_idx], #b3_ddg
            gradient_s$s_dgwt + gradient_b1$s_dgwt, #s1_dgwt
            gradient_b2$s_dgwt + 
                gradient_b3$s_dgwt, #s23_dgwt
            gradient_b1$b_dgwt, #b1_dgwt
            gradient_b2$b_dgwt, #b2_dgwt
            gradient_b3$b_dgwt, #b3_dgwt
            gradient_s$s_fwt, #s_fwt
            gradient_s$s_f0, #s_f0
            gradient_b1$b_fwt, #b1_fwt
            gradient_b1$b_f0, #b1_f0
            gradient_b2$b_fwt, #b2_fwt
            gradient_b2$b_f0, #b2_f0
            gradient_b3$b_fwt, #b3_fwt
            gradient_b3$b_f0 #b3_f0
        )
    }

    names(gradient) <- names(parameters)

    #set gradients from fixed parameters to 0
    for (i in seq_along(fixed_par)) {
        gradient[names(gradient) == names(fixed_par)[i]] <- 0
    }
    return(gradient)
}

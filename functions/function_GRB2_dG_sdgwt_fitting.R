function_GRB2_dG_sdgwt_fitting <- function(
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
  
  ## stability phenotype
  if (sum(!is.na(fw_list$s_fitness)) > 0) { #is zero for method 3
    s_fpred <- function_folding_dG2F(
      s_ddg = s_ddg_var, 
      s_dgwt = parameters[all_idl + 5], 
      s_fwt = parameters[all_idl + 6],
      s_f0 = parameters[all_idl + 7]
    )
    msd_stab <- sum(((fw_list$s_fitness - s_fpred) / fw_list$s_sigma)^2,
      na.rm = T)
  } else {
    msd_stab <- 0
  }
  
  ## binding phenotype
  b_fpred <- function_binding_dG2F(
    s_ddg = s_ddg_var,
    b_ddg = b_ddg_var,
    b_dgwt = parameters[all_idl + 1], 
    b_fwt = parameters[all_idl + 2],
    b_f0 = parameters[all_idl + 3],
    s_dgwt = parameters[all_idl + 4]
  )
  msd_bind <- sum(((fw_list$b_fitness - b_fpred) / fw_list$b_sigma)^2,
      na.rm = T)

  #if binding phenotype when ddg = 0
  if (predict_binding0 == TRUE) {
    b0_fpred <- function_binding_dG2F(
      s_ddg = s_ddg_var, 
      b_ddg = 0,
      b_dgwt = parameters[all_idl + 1], 
      b_fwt = parameters[all_idl + 2],
      b_f0 = parameters[all_idl + 3],
      s_dgwt = parameters[all_idl + 4]
    )
    msd_bind0 <- sum(((fw_list$b_fitness - b0_fpred) / fw_list$b_sigma)^2,
      na.rm = T)
  } else {
    msd_bind0 <- 0
  }
  
  
  ## mean square deviation
  msd <- msd_stab + msd_bind + msd_bind0
  return(msd)
}
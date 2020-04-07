function_PDZ_dG_fitting <- function(
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
  
  ## stability phenotype
  s_fpred <- function_folding_dG2F(
    s_ddg = s_ddg_var, 
    s_dgwt = parameters[all_idl + 1], 
    s_fwt = parameters[all_idl + 6],
    s_f0 = parameters[all_idl + 7]
  )
  msd_s <- sum(((fw_list$s_fitness - s_fpred) / fw_list$s_sigma)^2,
      na.rm = T)
  
  ## binding1 phenotype (ddPCA measurement)
  b1_fpred <- function_binding_dG2F(
    s_ddg = s_ddg_var, 
    b_ddg = b12_ddg_var,
    s_dgwt = parameters[all_idl + 1],
    b_dgwt = parameters[all_idl + 3], 
    b_fwt = parameters[all_idl + 8],
    b_f0 = parameters[all_idl + 9]
  )
  msd_b1 <- sum(((fw_list$b1_fitness - b1_fpred) / fw_list$b1_sigma)^2,
      na.rm = T)

  if (predict_binding0 == TRUE) {
    #if b12_ddG are 0
    b1_0_fpred <- function_binding_dG2F(
      s_ddg = s_ddg_var, 
      b_ddg = 0,
      s_dgwt = parameters[all_idl + 1],
      b_dgwt = parameters[all_idl + 3], 
      b_fwt = parameters[all_idl + 8],
      b_f0 = parameters[all_idl + 9]
    )
    msd_b10 <- sum(((fw_list$b1_fitness - b1_0_fpred) / fw_list$b1_sigma)^2,
      na.rm = T)

  }

  ## binding2 phenotype (McLaughlin CRIPT measurement)
  b2_fpred <- function_binding_dG2F(
    s_ddg = s_ddg_var, 
    b_ddg = b12_ddg_var,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 4], 
    b_fwt = parameters[all_idl + 10],
    b_f0 = parameters[all_idl + 11]
  )
  msd_b2 <- sum(((fw_list$b2_fitness - b2_fpred) / fw_list$b2_sigma)^2,
      na.rm = T)


  if (predict_binding0 == TRUE) {
    #if b12_ddG are 0
    b2_0_fpred <- function_binding_dG2F(
      s_ddg = s_ddg_var, 
      b_ddg = 0,
      s_dgwt = parameters[all_idl + 2],
      b_dgwt = parameters[all_idl + 4], 
      b_fwt = parameters[all_idl + 10],
      b_f0 = parameters[all_idl + 11]
    )
    msd_b20 <- sum(((fw_list$b2_fitness - b2_0_fpred) / fw_list$b2_sigma)^2,
      na.rm = T)
  }

  ## binding3 phenotype (McLaughlin Tm2F measurement), assuming b3 = 0
  b3_fpred <- function_binding_dG2F(
    s_ddg = s_ddg_var, 
    b_ddg = b3_ddg_var,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 5], 
    b_fwt = parameters[all_idl + 12],
    b_f0 = parameters[all_idl + 13]
  )
  msd_b3 <- sum(((fw_list$b3_fitness - b3_fpred) / fw_list$b3_sigma)^2,
      na.rm = T)

  if (predict_binding0 == TRUE) {
    #if b_12ddG are 0
    b3_0_fpred <- function_binding_dG2F(
      s_ddg = s_ddg_var, 
      b_ddg = 0,
      s_dgwt = parameters[all_idl + 2],
      b_dgwt = parameters[all_idl + 5], 
      b_fwt = parameters[all_idl + 12],
      b_f0 = parameters[all_idl + 13]
    )
    msd_b30 <- sum(((fw_list$b3_fitness - b3_0_fpred) / fw_list$b3_sigma)^2,
      na.rm = T)
  }

  ## mean square deviation

  if (predict_binding0 == TRUE) {
    msd <- msd_s + msd_b1 + msd_b2 + msd_b3 +
      msd_b10 + msd_b20 + msd_b30
  } else {
    msd <- msd <- msd_s + msd_b1 + msd_b2 + msd_b3
  }

  return(msd)
}
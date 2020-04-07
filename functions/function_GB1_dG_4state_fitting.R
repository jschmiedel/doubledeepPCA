function_GB1_dG_4state_fitting <- function(
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
  
  ## binding phenotype
  b_fpred <- function_dG_4state_binding_dg2logf(
    s1_ddg = s1_ddg_var,
    s2_ddg = s2_ddg_var,
    b_ddg = b_ddg_var,
    b_dgwt = parameters[all_idl + 1], 
    b_fwt = parameters[all_idl + 2],
    b_f0 = parameters[all_idl + 3],
    s1_dgwt = parameters[all_idl + 4],
    s2_dgwt = parameters[all_idl + 5]
  )
  msd_bind <- sum(((fw_list$b_fitness - b_fpred) / fw_list$b_sigma)^2,
      na.rm = T)

  #if binding phenotype when ddg = 0
  if (predict_binding0 == TRUE) {
    b0_fpred <- function_dG_4state_binding_dg2logf(
      s1_ddg = s1_ddg_var,
      s2_ddg = s2_ddg_var,
      b_ddg = 0,
      b_dgwt = parameters[all_idl + 1], 
      b_fwt = parameters[all_idl + 2],
      b_f0 = parameters[all_idl + 3],
      s1_dgwt = parameters[all_idl + 4],
      s2_dgwt = parameters[all_idl + 5]
    )
    msd_bind0 <- sum(((fw_list$b_fitness - b0_fpred) / fw_list$b_sigma)^2,
      na.rm = T)
  } else {
    msd_bind0 <- 0
  }
  
  if (approach == 2) {
    #comparision to in vitro folding ddgs
    #assuming ddgs measured in vitro correspond to sum of folded sets 1&2
    rt = 1.99e-3 * 310.15
    s_ddg_wt = -log(exp(-parameters[all_idl + 4]/rt) + exp(-parameters[all_idl + 5]/rt)) * rt
    s_ddg_mut = -log(exp(-(s1_ddg_mut + parameters[all_idl + 4])/rt) + 
                      exp(-(s2_ddg_mut + parameters[all_idl + 5])/rt)) * rt - s_ddg_wt

    sd_sddg = ((s_ddg_mut - fw_list$s_ddg_invitro) / fw_list$s_ddg_invitro_sd)^2
    sd_sddg[s_ddg_mut > 4 & fw_list$s_ddg_invitro == 4] = 0 #inequality for ddG == 4; because not known if larger than that
    msd_sddg = sum(sd_sddg, na.rm = T)
  } else {
    msd_sddg = 0
  }
  
  ## mean square deviation
  msd <- msd_bind + msd_bind0 + msd_sddg
  return(msd)
}
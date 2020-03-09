function_PDZ_dG_fitting <- function(
  parameters, 
  idlist,
  dt
) {

  s_dG_par <- parameters[1:idlist$s]
  b12_dG_par <- parameters[(idlist$s + 1):(idlist$s + idlist$b12)] 
  # b3_dG_par <- parameters[(idlist$s + idlist$b12 + 1):(idlist$s + idlist$b12 + idlist$b3)] 
  
  # all_idl <- idlist$s + idlist$b12 + idlist$b3
  all_idl <- idlist$s + idlist$b12

  ## stability phenotype
  dt[type == "stability", fitness_pred := function_folding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    s_dgwt = parameters[all_idl + 1], 
    s_fwt = parameters[all_idl + 6],
    s_f0 = parameters[all_idl + 7] 
  )]
  
  ## binding1 phenotype (ddPCA measurement)
  dt[type == "binding1", fitness_pred := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = b12_dG_par[id1_key],
    s_dgwt = parameters[all_idl + 1],
    b_dgwt = parameters[all_idl + 3], 
    b_fwt = parameters[all_idl + 8],
    b_f0 = parameters[all_idl + 9]
    # b_fwt = global_par[3],
    # b_f0 = global_par[4]
  )]
#if b12_ddG are 0
  dt[type == "binding1", fitness_pred_b0 := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 1],
    b_dgwt = parameters[all_idl + 3], 
    b_fwt = parameters[all_idl + 8],
    b_f0 = parameters[all_idl + 9]
    # b_fwt = global_par[3],
    # b_f0 = global_par[4]
  )]
  
  ## binding2 phenotype (McLaughlin CRIPT measurement)
  dt[type == "binding2", fitness_pred := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = b12_dG_par[id1_key],
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 4], 
    b_fwt = parameters[all_idl + 10],
    b_f0 = parameters[all_idl + 11]
    # b_fwt = global_par[5],
    # b_f0 = global_par[6]
  )]
#if b12_ddG are 0
  dt[type == "binding2", fitness_pred_b0 := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 4], 
    b_fwt = parameters[all_idl + 10],
    b_f0 = parameters[all_idl + 11]
    # b_fwt = global_par[5],
    # b_f0 = global_par[6]
  )]

  ## binding3 phenotype (McLaughlin Tm2F measurement), assuming b3 = 0
  dt[type == "binding3", fitness_pred := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = b12_dG_par[id1_key],
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 5], 
    b_fwt = parameters[all_idl + 12],
    b_f0 = parameters[all_idl + 13]
    # b_fwt = global_par[7],
    # b_f0 = global_par[8]
  )]

  #if b_12ddG are 0
  dt[type == "binding3", fitness_pred_b0 := function_binding_dG2F(
    s_ddg = s_dG_par[id1_key], 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 5], 
    b_fwt = parameters[all_idl + 12],
    b_f0 = parameters[all_idl + 13]
    # b_fwt = global_par[7],
    # b_f0 = global_par[8]
  )]

  ## mean square deviation
  msd <- dt[, sum(((fitness - fitness_pred) / sigma)^2, na.rm = T)] +
    dt[grep("binding", type), sum(((fitness - fitness_pred_b0) / sigma)^2,
     na.rm = T)]
  return(msd)
}












function_PDZ_dG_fitting2 <- function(
  parameters,
  fixed_par,
  idlist,
  fw_list
) {

  for (i in seq_along(fixed_par)) {
    parameters[names(parameters) == names(fixed_par)[i]] <- fixed_par[i]
  }

  s_ddg <- parameters[1:idlist$s]
  b12_ddg <- parameters[(idlist$s + 1):(idlist$s + idlist$b12)]
  b3_ddg <- parameters[(idlist$s + idlist$b12 + 1):(idlist$s + idlist$b12 + idlist$b3)]
   
  all_idl <- idlist$s + idlist$b12 + idlist$b3
  # all_idl <- idlist$s + idlist$b12
  
  ## stability phenotype
  s_fpred <- function_folding_dG2F(
    s_ddg = s_ddg, 
    s_dgwt = parameters[all_idl + 1], 
    s_fwt = parameters[all_idl + 6],
    s_f0 = parameters[all_idl + 7]
  )
  
  ## binding1 phenotype (ddPCA measurement)
  b1_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = b12_ddg,
    s_dgwt = parameters[all_idl + 1],
    b_dgwt = parameters[all_idl + 3], 
    b_fwt = parameters[all_idl + 8],
    b_f0 = parameters[all_idl + 9]
  )

#if b12_ddG are 0
  b1_0_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 1],
    b_dgwt = parameters[all_idl + 3], 
    b_fwt = parameters[all_idl + 8],
    b_f0 = parameters[all_idl + 9]
  )
  
  ## binding2 phenotype (McLaughlin CRIPT measurement)
  b2_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = b12_ddg,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 4], 
    b_fwt = parameters[all_idl + 10],
    b_f0 = parameters[all_idl + 11]
  )

  #if b12_ddG are 0
  b2_0_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 4], 
    b_fwt = parameters[all_idl + 10],
    b_f0 = parameters[all_idl + 11]
  )

  ## binding3 phenotype (McLaughlin Tm2F measurement), assuming b3 = 0
  b3_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = b3_ddg,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 5], 
    b_fwt = parameters[all_idl + 12],
    b_f0 = parameters[all_idl + 13]
  )

  #if b_12ddG are 0
  b3_0_fpred <- function_binding_dG2F(
    s_ddg = s_ddg, 
    b_ddg = 0,
    s_dgwt = parameters[all_idl + 2],
    b_dgwt = parameters[all_idl + 5], 
    b_fwt = parameters[all_idl + 12],
    b_f0 = parameters[all_idl + 13]
  )

  ## mean square deviation
  msd <- sum(
    c(
      ((fw_list$s_fitness - s_fpred) / fw_list$s_sigma) ^ 2,
      ((fw_list$b1_fitness - b1_fpred) / fw_list$b1_sigma) ^ 2,
      ((fw_list$b1_fitness - b1_0_fpred) / fw_list$b1_sigma) ^ 2,
      ((fw_list$b2_fitness - b2_fpred) / fw_list$b2_sigma) ^ 2,
      ((fw_list$b2_fitness - b2_0_fpred) / fw_list$b2_sigma) ^ 2,
      ((fw_list$b3_fitness - b3_fpred) / fw_list$b3_sigma) ^ 2,
      ((fw_list$b3_fitness - b3_0_fpred) / fw_list$b3_sigma) ^ 2
    ),
    na.rm = T
  )
  return(msd)
}
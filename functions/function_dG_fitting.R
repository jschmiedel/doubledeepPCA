function_dG_fitting <- function(
  parameters, 
  s_idl = s_idl,
  b_idl = b_idl,
  dt,
  m1_mod = 1
) {

  s_dG_par <- parameters[1:s_idl]
  if (m1_mod != 0) {
    b_dG_par <- parameters[(s_idl + 1):(s_idl + b_idl)] 
  }
  
  ## stability phenotype
  if (dt[type == "stability", .N] > 0) {
    dt[type == "stability", fitness_pred := function_folding_dG2F(
      s_ddg = ifelse(
        is.na(id2_key),
        s_dG_par[id1_key],
        s_dG_par[id1_key] + s_dG_par[id2_key]
      ), 
      s_dgwt = parameters[s_idl + m1_mod * b_idl + 1], 
      s_fwt = parameters[s_idl + m1_mod * b_idl + 5],
      s_f0 = parameters[s_idl + m1_mod * b_idl + 6] 
    )]
  }
  ## binding phenotype
  if (m1_mod == 0) { # first iteration of method 1 fitting global rel.
    dt[type == "binding", b_ddg := 0]
  } else {
    dt[type == "binding", b_ddg := ifelse(
      is.na(id2_key),
      b_dG_par[id1_key],
      b_dG_par[id1_key] + b_dG_par[id2_key]
    )]
  }
    
  dt[type == "binding", fitness_pred := function_binding_dG2F(
    s_ddg = ifelse(
      is.na(id2_key),
      s_dG_par[id1_key],
      s_dG_par[id1_key] + s_dG_par[id2_key]
    ), 
    b_ddg = b_ddg,
    s_dgwt = parameters[s_idl + m1_mod * b_idl + 1],
    b_dgwt = parameters[s_idl + m1_mod * b_idl + 2], 
    b_fwt = parameters[s_idl + m1_mod * b_idl + 3],
    b_f0 = parameters[s_idl + m1_mod * b_idl + 4]
  )]
  
  ## mean square deviation; divide by weights to correct for variants being NA
  MSD <- dt[,sum((fitness - fitness_pred)^2 / sigma^2, na.rm = T) / 
    sum(sigma^-2, na.rm = T)]
  return(MSD)
}
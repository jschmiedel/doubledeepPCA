function_dG_method3_fitting <- function(
  parameters, 
  id_L, 
  list_fs, 
  list_keys
) {
  s_dG_par <- parameters[1:id_L]
  b_dG_par <- parameters[(id_L + 1):(2 * id_L)]
  # binding phenotype
  bf <- function_binding_dG2F(
    b_dG = c(
      b_dG_par[list_keys$singles_key],
      b_dG_par[list_keys$doubles_key1] + b_dG_par[list_keys$doubles_key2]
    ), 
    s_dG = c(
      s_dG_par[list_keys$singles_key],
      s_dG_par[list_keys$doubles_key1] + s_dG_par[list_keys$doubles_key2]
    ), 
    b_dGwt = parameters[2 * id_L + 4], 
    s_dGwt = parameters[2 * id_L + 3], 
    b_bgr = parameters[2 * id_L + 2], 
    b_scale = parameters[2 * id_L + 1]
  )
  # mean square deviation
  MSD <- sum((list_fs$b_fitness - bf)^2 / list_fs$b_sigma^2, na.rm = T) / 
    sum(list_fs$b_sigma^-2, na.rm = T)
  return(MSD)
}
function_dG_method2_fitting <- function(
  parameters, 
  id_L, 
  list_fs, 
  list_keys
) {
  s_dG_par <- parameters[1:id_L]
  b_dG_par <- parameters[(id_L + 1):(2 * id_L)]
  
  ## stability phenotype
  sf <- function_folding_dG2F(
    s_dG = c(
      s_dG_par[list_keys$singles_key],
      s_dG_par[list_keys$doubles_key1] + s_dG_par[list_keys$doubles_key2]
    ), 
    s_dGwt = parameters[2 * id_L + 5], 
    s_bgr = parameters[2 * id_L + 3], 
    s_scale = parameters[2 * id_L + 1]
  )
  ## binding phenotype
  bf <- function_binding_dG2F(
    b_dG = c(
      b_dG_par[list_keys$singles_key],
      b_dG_par[list_keys$doubles_key1] + b_dG_par[list_keys$doubles_key2]
    ), 
    s_dG = c(
      s_dG_par[list_keys$singles_key],
      s_dG_par[list_keys$doubles_key1] + s_dG_par[list_keys$doubles_key2]
    ), 
    b_dGwt = parameters[2 * id_L + 6], 
    s_dGwt = parameters[2 * id_L + 5], 
    b_bgr = parameters[2 * id_L + 4], 
    b_scale = parameters[2 * id_L + 2]
  )
  ## mean square deviation; divide by weights to correct for variants being NA
  MSD <- sum((list_fs$b_fitness - bf)^2 / list_fs$b_sigma^2 + 
    (list_fs$s_fitness - sf)^2 / list_fs$s_sigma^2, na.rm = T) / 
    sum(list_fs$b_sigma^-2 + list_fs$s_sigma^-2, na.rm = T)
  return(MSD)
}



function_dG_method2_fitting = function(parameters,id_L,list_fs,list_keys) {
  
  s_dG_par = parameters[1:id_L]
  b_dG_par = parameters[(id_L+1):(2*id_L)]
  if (length(parameters) > (2*id_L)) {
    b_scale = parameters[2*id_L+1]
    s_scale = parameters[2*id_L+2]
    b_bgr = parameters[2*id_L+3]
    s_bgr = parameters[2*id_L+4]
  }
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness being 1 !!! change if using fitness simply from log-ratios of counts
  s_dGwt = function_folding_F2dG(1,s_bgr,s_scale)
  b_dGwt = function_binding_F2dG(1,s_dGwt,b_bgr,b_scale)
  
  #mix and match mutants and variants
  #stability phenotype
  s_dG = c(s_dG_par[list_keys[["singles_key"]]],s_dG_par[list_keys[["doubles_key1"]]]+s_dG_par[list_keys[["doubles_key2"]]])
  sf = function_folding_dG2F(s_dG,s_dGwt,s_bgr,s_scale)
  #binding phenotype
  b_dG = c(b_dG_par[list_keys[["singles_key"]]],b_dG_par[list_keys[["doubles_key1"]]]+b_dG_par[list_keys[["doubles_key2"]]])
  bf = function_binding_dG2F(b_dG,s_dG,b_dGwt,s_dGwt,b_bgr,b_scale)
  
  #some bf/sf values might be NA because they are below background growth, set fitness and error NA to include in deviation calcuation
  s_sigma = list_fs[["s_sigma"]]
  s_sigma[is.na(sf)] = NA
  s_sigma[is.infinite(sf)] = NA
  
  b_sigma = list_fs[["b_sigma"]]
  b_sigma[is.na(bf)] = NA
  b_sigma[is.infinite(bf)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  MSD = sum((list_fs[["b_fitness"]] - bf)^2 / b_sigma^2 + (list_fs[["s_fitness"]] - sf)^2 / s_sigma^2,na.rm=T) / sum(b_sigma^-2 + s_sigma^-2,na.rm=T)
  return(MSD)
}
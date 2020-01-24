function_dGestimation_method2_fitting = function(parameters) {
  
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
  s_dG_for_f = c(s_dG_par[s_singles_key],s_dG_par[s_doubles_key1]+s_dG_par[s_doubles_key2])
  sf = function_folding_dG2F(s_dG_for_f,s_dGwt,s_bgr,s_scale)
  #binding phenotype
  s_dG_for_b = c(s_dG_par[b_singles_key],s_dG_par[b_doubles_key1]+s_dG_par[b_doubles_key2])
  b_dG_for_b = c(b_dG_par[b_singles_key],b_dG_par[b_doubles_key1]+b_dG_par[b_doubles_key2])
  bf = function_binding_dG2F(b_dG_for_b,s_dG_for_b,b_dGwt,s_dGwt,b_bgr,b_scale)
  
  #some bf/sf values might be NA because they are below background growth, set fitness and error NA to include in deviation calcuation
  binding_fitness2 = matrix(binding_fitness) #get rid of this
  stability_fitness2 = matrix(stability_fitness)
  binding_error2 = matrix(binding_error)
  binding_error2[is.na(bf)] = NA
  binding_error2[is.infinite(bf)] = NA
  stability_error2 = matrix(stability_error)
  stability_error2[is.na(sf)] = NA
  stability_error2[is.infinite(sf)] = NA
  
  #mean square deviation; for both phenotypes separately
  # MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2,na.rm=T) / sum(binding_error2^-2,na.rm=T) +
    # sum((stability_fitness2 - sf)^2 / stability_error2^2,na.rm=T) / sum(stability_error2^-2,na.rm=T)
  #for both phenotypes together
  MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - sf)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
  return(MSD)
}
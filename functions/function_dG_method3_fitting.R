function_dG_method3_fitting = function(parameters,id_L,list_fs,list_keys) {
  
  s_dG_par = parameters[1:id_L]
  b_dG_par = parameters[(id_L+1):(2*id_L)]
  b_scale = parameters[2*id_L+1]
  b_bgr = parameters[2*id_L+2]
  s_dGwt = parameters[2*id_L+3]
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness being 1 !!! change if using fitness simply from log-ratios of counts
  b_dGwt = function_binding_F2dG(b_fitness=1,s_dG=s_dGwt,b_bgr=b_bgr,b_scale=b_scale)
  
  #mix and match mutants and variants
  #binding phenotype
  s_dG = c(s_dG_par[list_keys$singles_key],
            s_dG_par[list_keys$doubles_key1]+s_dG_par[list_keys$doubles_key2])
  b_dG = c(b_dG_par[list_keys$singles_key],
            b_dG_par[list_keys$doubles_key1]+b_dG_par[list_keys$doubles_key2])
  bf = function_binding_dG2F(b_dG=b_dG,s_dG=s_dG,b_dGwt=b_dGwt,s_dGwt=s_dGwt,b_bgr=b_bgr,b_scale=b_scale)
  
  #some bf values might be NA because they are below background growth, set fitness and error NA to include in deviation calcuation
  b_sigma = list_fs$b_sigma
  b_sigma[is.na(bf)] = NA
  b_sigma[is.infinite(bf)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  MSD = sum((list_fs$b_fitness - bf)^2 / b_sigma^2,na.rm=T) / sum(b_sigma^-2,na.rm=T)
  return(MSD)
}
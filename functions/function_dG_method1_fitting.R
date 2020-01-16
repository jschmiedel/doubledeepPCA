############################################################
#### fit global relationship between deltaGs and fitess from binding and stabilityPCA data assuming deltadeltaG binding = 0 for all used variants
############################################################
function_dG_method1_fitting = function(parameters,id_L,global_par,list_fs) {
  
  s_dG = parameters[1:id_L] #fit stability delta Gs
  if (length(parameters) == (id_L+2)) { # + scale parameters [>> determine global relationship; background parameters set manually]
    s_scale = parameters[id_L+1]
    b_scale = parameters[id_L+2]
    s_bgr = global_par[1]
    b_bgr = global_par[2]
    b_dG = 0
  } else if (length(parameters) == (id_L+4)) { # + scale and background parameters [>> determine global relationship]
    s_scale = parameters[id_L+1]
    b_scale = parameters[id_L+2]
    s_bgr = parameters[id_L+3]
    b_bgr = parameters[id_L+4]
    b_dG = 0
  } else if (length(parameters) == 2*id_L) { # + binding delta Gs [>> determine dGs after global relationship was already determined]
    b_dG = parameters[(id_L + 1) : (2*id_L)]
    s_scale = global_par[1]
    b_scale = global_par[2]
    s_bgr = global_par[3]
    b_bgr = global_par[4]
  }
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness == 1!
  s_dGwt = function_folding_F2dG(1,s_bgr,s_scale)
  b_dGwt = function_binding_F2dG(1,s_dGwt,b_bgr,b_scale)
  
  bf = function_binding_dG2F(b_dG,s_dG,b_dGwt,s_dGwt,b_bgr,b_scale)
  sf = function_folding_dG2F(s_dG,s_dGwt,s_bgr,s_scale)
  
  #some bf/sf values might be NA because they are below background growth, set fitness and error NA to include in deviation calcuation
  s_sigma = list_fs$s_sigma
  s_sigma[is.na(sf)] = NA
  s_sigma[is.infinite(sf)] = NA
  
  b_sigma = list_fs$b_sigma
  b_sigma[is.na(bf)] = NA
  b_sigma[is.infinite(bf)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  MSD = sum((list_fs$b_fitness - bf)^2 / b_sigma^2 + (list_fs$s_fitness - sf)^2 / s_sigma^2,na.rm=T) / sum(b_sigma^-2 + s_sigma^-2,na.rm=T)
  return(MSD)
}
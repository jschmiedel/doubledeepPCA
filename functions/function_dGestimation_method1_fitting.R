############################################################
#### fit global relationship between deltaGs and fitess from binding and stabilityPCA data assuming deltadeltaG binding = 0 for all used variants
############################################################
function_dGestimation_method1_fitting = function(parameters) {
  
  s_dG = parameters[1:id_L] #fit stability delta Gs
  if (length(parameters) == (id_L+2)) { # + scale parameters [>> determine global relationship; background parameters set manually]
    b_scale = parameters[id_L+1]
    s_scale = parameters[id_L+2]
    b_dG = 0
  } else if (length(parameters) == (id_L+4)) { # + scale and background parameters [>> determine global relationship]
    b_scale = parameters[id_L+1]
    s_scale = parameters[id_L+2]
    b_bgr = parameters[id_L+3]
    s_bgr = parameters[id_L+4]
    b_dG = 0
  } else if (length(parameters) == 2*id_L) { # + binding delta Gs [>> determine dGs after global relationship was already determined]
    b_dG = parameters[(id_L + 1) : (2*id_L)]
  } 
  # else if (length(parameters) == 2*id_L+2) { # + binding delta Gs and scale parameters
  #   b_dG = parameters[(id_L + 1) : (2*id_L)]
  #   b_scale = parameters[2*id_L+1]
  #   s_scale = parameters[2*id_L+2]
  # } else if (length(parameters) == 2*id_L+4) { # + binding delta Gs, scale and background parameters
  #   b_dG = parameters[(id_L + 1) : (2*id_L)]
  #   b_scale = parameters[2*id_L+1]
  #   s_scale = parameters[2*id_L+2]
  #   b_bgr = parameters[2*id_L+3]
  #   s_bgr = parameters[2*id_L+4]
  # }
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness == 1!
  s_dGwt = function_folding_F2dG(1,s_bgr,s_scale)
  b_dGwt = function_binding_F2dG(1,s_dGwt,b_bgr,b_scale)
  
  bf = function_binding_dG2F(b_dG,s_dG,b_dGwt,s_dGwt,b_bgr,b_scale)
  ff = function_folding_dG2F(s_dG,s_dGwt,s_bgr,s_scale)
  
  #some bf/ff values might be NA because they are below background growth, set fitness and error NA to include in deviation calcuation
  binding_fitness2 = binding_fitness
  stability_fitness2 = stability_fitness
  binding_error2 = binding_error
  binding_error2[is.na(bf)] = NA
  binding_error2[is.infinite(bf)] = NA
  stability_error2 = stability_error
  stability_error2[is.na(ff)] = NA
  stability_error2[is.infinite(ff)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
  return(MSD)
}
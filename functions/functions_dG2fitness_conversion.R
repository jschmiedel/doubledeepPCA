#functions to convert dG to fitness and vice versa for both stability and binding phenotypes

function_folding_dG2F = function(s_dG,s_dGwt,s_bgr,s_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  sf = (1-s_bgr)/s_scale/(1+exp((s_dGwt+s_dG)/R/Temp)) + s_bgr
}
function_binding_dG2F = function(b_dG,s_dG,b_dGwt,s_dGwt,b_bgr,b_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  bf = (1-b_bgr)/b_scale/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((s_dGwt+s_dG)/R/Temp))) + b_bgr
}
function_folding_F2dG = function(s_fitness,s_bgr,s_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  s_dG = R*Temp*log((1-s_bgr)/s_scale/(s_fitness-s_bgr)-1)
}
function_binding_F2dG = function(b_fitness,s_dG,b_bgr,b_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  b_dG = R*Temp*(log((1-b_bgr)/b_scale/(b_fitness-b_bgr)-1) - log(1+exp(s_dG/R/Temp)))
}
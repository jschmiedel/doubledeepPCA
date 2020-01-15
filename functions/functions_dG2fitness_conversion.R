#functions to convert dG to fitness and vice versa for both stability and binding phenotypes

function_folding_dG2F = function(f_dG,f_dGwt,f_bgr,f_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  ff = (1-f_bgr)/f_scale/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr
}
function_binding_dG2F = function(b_dG,f_dG,b_dGwt,f_dGwt,b_bgr,b_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  bf = (1-b_bgr)/b_scale/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr
}
function_folding_F2dG = function(f_fitness,f_bgr,f_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  f_dG = R*Temp*log((1-f_bgr)/f_scale/(f_fitness-f_bgr)-1)
}
function_binding_F2dG = function(b_fitness,f_dG,b_bgr,b_scale) {
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  b_dG = R*Temp*(log((1-b_bgr)/b_scale/(b_fitness-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
}
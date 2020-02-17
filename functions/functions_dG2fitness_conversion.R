## functions to convert dG to fitness and vice versa
## for both stability and binding phenotypes

function_folding_dG2F <- function(
  # s_dG,
  # s_dGwt,
  # s_bgr,
  # s_scale,
  s_ddg, 
  s_dgwt, 
  s_fwt,
  s_f0,  
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  # sf <- (1 - s_bgr) / s_scale / (1 + exp((s_dGwt + s_dG) / r / temp)) + s_bgr

  scale <- (1 + exp(s_dgwt / r / temp)) * (s_fwt - s_f0)
  sf <- s_f0 + scale / (1 + exp(s_dgwt + s_ddg))
}

function_binding_dG2F <- function(
  # b_dG, 
  # s_dG, 
  # b_dGwt, 
  # s_dGwt, 
  # b_bgr, 
  # b_scale,
  b_ddg, 
  s_ddg, 
  b_dgwt, 
  s_dgwt, 
  b_fwt,
  b_f0, 
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  # bf <- (1 - b_bgr) / b_scale / (1 + exp((b_dGwt + b_dG) / r / temp) * 
  #   (1 + exp((s_dGwt + s_dG) / r / temp))) + b_bgr
  scale <- (1 + exp(b_dgwt / r / temp) * 
    (1 + exp(s_dgwt / r / temp))) * (b_fwt - b_f0)
  bf <- b_f0 + scale / (1 + exp((b_dgwt + b_ddg) / r / temp) * 
    (1 + exp((s_dgwt + s_ddg) / r / temp)))
}

function_folding_F2dG <- function(
  # s_fitness, 
  # s_bgr, 
  # s_scale,
  s_fitness, 
  s_fwt, 
  s_f0,
  s_dgwt,
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  # s_dG <- r * temp * log((1 - s_bgr) / s_scale / (s_fitness - s_bgr) - 1)
  
  scale <- (1 + exp(s_dgwt / r / temp)) * (s_fwt - s_f0)
  s_dG <- r * temp * log(scale / (s_fitness - s_f0) - 1)
}

function_binding_F2dG <- function(
  b_fitness, 
  # s_dG, 
  # b_bgr, 
  # b_scale,
  s_ddg, 
  b_dgwt,
  s_dgwt,
  b_f0, 
  b_fwt,
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  # b_dG <- r * temp * 
  #   (log((1 - b_bgr) / b_scale / (b_fitness - b_bgr) - 1) - 
  #     log(1 + exp(s_dG / r / temp)))

  scale <- (1 + exp(b_dgwt / r / temp) * 
    (1 + exp(s_dgwt / r / temp))) * (b_fwt - b_f0)
  b_dG <- r * temp * 
    (log(scale / (b_fitness - b_f0) - 1) - 
      log(1 + exp(s_ddg / r / temp)))
}
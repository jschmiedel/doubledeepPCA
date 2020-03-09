## functions to convert dG to fitness and vice versa
## for both stability and binding phenotypes

function_folding_dG2F <- function(
  s_ddg, 
  s_dgwt, 
  s_fwt,
  s_f0,
  rt = 1.98e-3 * 310.15
) {
  sf <- s_f0 + 
    (s_fwt - s_f0) * (1 + exp(s_dgwt / rt)) /
      (1 + exp((s_dgwt + s_ddg) / rt))
}

function_stability_dG2F_gradient <- function(
  s_ddg, 
  s_dgwt, 
  s_fwt,
  s_f0,  
  f,
  w,
  n = 2,
  rt = 1.98e-3 * 310.15
) {
  cf <- (s_fwt - s_f0) / rt
  ei <- exp(s_ddg / rt)
  ewt <- exp(s_dgwt / rt)
  eiwt <- ei * ewt
  df_ds_ddg <- -cf * eiwt * (1 + ewt) / (1 + eiwt)^2
  df_ds_dgwt <- cf * ewt * (1 - ei) / (1 + eiwt)^2
  df_ds_fwt <- (1 + ewt) / (1 + eiwt)
  df_ds_f0 <- ewt * (ei - 1) / (1 + eiwt)
  f_pred <- function_folding_dG2F(
    s_ddg = s_ddg, 
    s_dgwt = s_dgwt, 
    s_fwt = s_fwt,
    s_f0 = s_f0
  )
  derr_df <- -n * (f - f_pred)^ (n - 1) / w^n
  gradient_s_ddg <- derr_df * df_ds_ddg
  gradient_s_dgwt <- sum(derr_df * df_ds_dgwt)
  gradient_s_fwt <- sum(derr_df * df_ds_fwt)
  gradient_s_f0 <- sum(derr_df * df_ds_f0)
  
  return(list(
    s_ddg = gradient_s_ddg,
    s_dgwt = gradient_s_dgwt,
    s_fwt = gradient_s_fwt,
    s_f0 = gradient_s_f0
  ))
}

function_binding_dG2F <- function(
  b_ddg, 
  s_ddg, 
  b_dgwt, 
  s_dgwt, 
  b_fwt,
  b_f0,
  rt = 1.98e-3 * 310.15
) {
  bf <- b_f0 + 
    (b_fwt - b_f0) * (1 + exp(b_dgwt / rt) * (1 + exp(s_dgwt / rt))) /
     (1 + exp((b_dgwt + b_ddg) / rt) * (1 + exp((s_dgwt + s_ddg) / rt)))
}

function_binding_dG2F_gradient <- function(
  s_ddg,  
  b_ddg,
  s_dgwt,
  b_dgwt,
  b_fwt,
  b_f0, 
  f,
  w,
  n = 2,
  rt = 1.98e-3 * 310.15
) {
  cf <- (b_fwt - b_f0) / rt
  ebi <- exp(b_ddg / rt)
  ebwt <- exp(b_dgwt / rt)
  esi <- exp(s_ddg / rt)
  eswt <- exp(s_dgwt / rt)

  ebiwt <- ebi * ebwt
  esiwt <- esi * eswt
  eiwt <- ebiwt * (1 + esiwt)
  ewt <- ebwt * (1 + eswt)

  df_ds_ddg <- -cf * ebiwt * esiwt * (1 + ewt) /
    (1 + eiwt) ^ 2
  df_db_ddg <- -cf * (1 + ewt) * eiwt /
    (1 + eiwt) ^ 2
  df_ds_dgwt <- cf * ebwt * eswt * (1 + eiwt - ebi * esi * (1 + ewt)) / 
    (1 + eiwt) ^ 2
  df_db_dgwt <- cf * (ewt + eiwt) /
    (1 + eiwt) ^ 2
  df_db_fwt <- (1 + ewt) / (1 + eiwt)
  df_db_f0 <- (eiwt - ewt) / (1 + eiwt)
  f_pred <- function_binding_dG2F(
    s_ddg = s_ddg,
    b_ddg = b_ddg,
    s_dgwt = s_dgwt,  
    b_dgwt = b_dgwt, 
    b_fwt = b_fwt,
    b_f0 = b_f0
  )

  derr_df <- -n * (f - f_pred)^ (n - 1) / w^n
  gradient_s_ddg <- derr_df * df_ds_ddg
  gradient_b_ddg <- derr_df * df_db_ddg
  gradient_s_dgwt <- sum(derr_df * df_ds_dgwt)
  gradient_b_dgwt <- sum(derr_df * df_db_dgwt)
  gradient_b_fwt <- sum(derr_df * df_db_fwt)
  gradient_b_f0 <- sum(derr_df * df_db_f0)
  
  return(list(
    s_ddg = gradient_s_ddg, 
    b_ddg = gradient_b_ddg, 
    s_dgwt = gradient_s_dgwt,
    b_dgwt = gradient_b_dgwt, 
    b_fwt = gradient_b_fwt, 
    b_f0 = gradient_b_f0
  ))
}

function_folding_F2dG <- function(
  s_fitness, 
  s_fwt, 
  s_f0,
  s_dgwt,
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  scale <- (1 + exp(s_dgwt / r / temp)) * (s_fwt - s_f0)
  s_dG <- r * temp * log(scale / (s_fitness - s_f0) - 1)
}

function_binding_F2dG <- function(
  b_fitness,
  s_ddg, 
  b_dgwt,
  s_dgwt,
  b_f0, 
  b_fwt,
  r = 1.98 * 10^(-3), # kcal/mol
  temp = 310.15
) {
  scale <- (1 + exp(b_dgwt / r / temp) * 
    (1 + exp(s_dgwt / r / temp))) * (b_fwt - b_f0)
  b_dG <- r * temp * 
    (log(scale / (b_fitness - b_f0) - 1) - 
      log(1 + exp(s_ddg / r / temp)))
}
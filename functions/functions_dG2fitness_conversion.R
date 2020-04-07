## functions to convert dG to fitness and vice versa
## for both stability and binding phenotypes

function_folding_dG2F <- function(
  s_ddg, 
  s_dgwt, 
  s_fwt,
  s_f0,
  rt = 1.99e-3 * 310.15
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
  mutXvar,
  rt = 1.99e-3 * 310.15,
  logF = FALSE
) {
  cf <- (s_fwt - s_f0) / rt
  ei <- exp(s_ddg / rt)
  ewt <- exp(s_dgwt / rt)
  eiwt <- ei * ewt
  df_ds_ddg <- -cf * eiwt * (1 + ewt) / (1 + eiwt)^2 #
  df_ds_ddg[is.na(df_ds_ddg)] <- 0

  df_ds_dgwt <- cf * (ewt - eiwt) / (1 + eiwt)^2 #
  df_ds_fwt <- (1 + ewt) / (1 + eiwt) #
  df_ds_f0 <- (eiwt - ewt) / (1 + eiwt) #
  f_pred <- function_folding_dG2F(
    s_ddg = s_ddg, 
    s_dgwt = s_dgwt, 
    s_fwt = s_fwt,
    s_f0 = s_f0
  )
  if (logF == FALSE) {
    derr_df <- -2 * (f - f_pred) / w^2
  } else {
    derr_df <- -2 * (f - log(f_pred)) / w^2 / f_pred
  }
  derr_df[is.na(derr_df)] <- 0
  gradient_s_ddg <- mutXvar %*% (derr_df * df_ds_ddg)
  gradient_s_dgwt <- sum(derr_df * df_ds_dgwt, na.rm = T)
  gradient_s_fwt <- sum(derr_df * df_ds_fwt, na.rm = T)
  gradient_s_f0 <- sum(derr_df * df_ds_f0, na.rm = T)
  
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
  rt = 1.99e-3 * 310.15
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
  mutXvar,
  rt = 1.99e-3 * 310.15,
  logF = FALSE
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
    (1 + eiwt) ^ 2 #
  df_ds_ddg[is.na(df_ds_ddg)] <- 0
  
  df_db_ddg <- -cf * (1 + ewt) * eiwt /
    (1 + eiwt) ^ 2 #
  df_db_ddg[is.na(df_db_ddg)] <- 0

  df_ds_dgwt <- cf * ebwt * eswt * (1 + eiwt - ebi * esi * (1 + ewt)) / 
    (1 + eiwt) ^ 2 #
  df_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!
  df_db_fwt <- (1 + ewt) / (1 + eiwt) #
  df_db_f0 <- (eiwt - ewt) / (1 + eiwt) #
  f_pred <- function_binding_dG2F(
    s_ddg = s_ddg,
    b_ddg = b_ddg,
    s_dgwt = s_dgwt,  
    b_dgwt = b_dgwt, 
    b_fwt = b_fwt,
    b_f0 = b_f0
  )
  if (logF == FALSE) {
    derr_df <- -2 * (f - f_pred) / w^2
  } else {
    derr_df <- -2 * (f - log(f_pred)) / w^2 / f_pred
  }
  derr_df[is.na(derr_df)] <- 0
  gradient_s_ddg <- mutXvar %*% (derr_df * df_ds_ddg)
  gradient_b_ddg <- mutXvar %*% (derr_df * df_db_ddg)
  gradient_s_dgwt <- sum(derr_df * df_ds_dgwt, na.rm = T)
  gradient_b_dgwt <- sum(derr_df * df_db_dgwt, na.rm = T)
  gradient_b_fwt <- sum(derr_df * df_db_fwt, na.rm = T)
  gradient_b_f0 <- sum(derr_df * df_db_f0, na.rm = T)
  
  return(list(
    s_ddg = gradient_s_ddg, 
    b_ddg = gradient_b_ddg, 
    s_dgwt = gradient_s_dgwt,
    b_dgwt = gradient_b_dgwt, 
    b_fwt = gradient_b_fwt, 
    b_f0 = gradient_b_f0
  ))
}

# function_folding_F2dG <- function(
#   s_fitness, 
#   s_fwt, 
#   s_f0,
#   s_dgwt,
#   r = 1.99 * 10^(-3), # kcal/mol
#   temp = 310.15
# ) {
#   scale <- (1 + exp(s_dgwt / r / temp)) * (s_fwt - s_f0)
#   s_dG <- r * temp * log(scale / (s_fitness - s_f0) - 1)
# }

# function_binding_F2dG <- function(
#   b_fitness,
#   s_ddg, 
#   b_dgwt,
#   s_dgwt,
#   b_f0, 
#   b_fwt,
#   r = 1.99 * 10^(-3), # kcal/mol
#   temp = 310.15
# ) {
#   scale <- (1 + exp(b_dgwt / r / temp) * 
#     (1 + exp(s_dgwt / r / temp))) * (b_fwt - b_f0)
#   b_dG <- r * temp * 
#     (log(scale / (b_fitness - b_f0) - 1) - 
#       log(1 + exp(s_ddg / r / temp)))
# }


function_dG_4state_folding_dg2fitness <- function(
  s1_ddg,
  s2_ddg,
  s1_dgwt, 
  s2_dgwt, 
  s_fwt,
  s_f0,
  rt = 1.99e-3 * 310.15
) {
  sf <- s_f0 + 
    (s_fwt - s_f0) * 
    (exp(-(s1_dgwt + s1_ddg) / rt) + exp(-(s2_dgwt + s2_ddg) / rt)) * 
    (1 + exp(-s1_dgwt / rt) + exp(-s2_dgwt / rt)) /
      ((1 + exp(-(s1_dgwt + s1_ddg) / rt) + exp(-(s2_dgwt + s2_ddg) / rt)) *
      (exp(-s1_dgwt / rt) + exp(-s2_dgwt / rt))) 
}

function_dG_4state_folding_dg2fitness_gradient <- function(
  s1_ddg, 
  s2_ddg,
  s1_dgwt, 
  s2_dgwt,
  s_fwt,
  s_f0,  
  f,
  w,
  mutXvar,
  rt = 1.99e-3 * 310.15
) {
  cf <- (s_fwt - s_f0) / rt

  e1i <- exp(-s1_ddg / rt)
  e2i <- exp(-s2_ddg / rt)

  e1wt <- exp(-s1_dgwt / rt)
  e2wt <- exp(-s2_dgwt / rt)

  e1iwt <- e1i * e1wt
  e2iwt <- e2i * e2wt

  ewt <- e1wt + e2wt
  eiwt <- e1iwt + e2iwt

  df_ds1_ddg <- -cf * e1iwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2 #
  df_ds1_ddg[is.na(df_ds1_ddg)] <- 0

  df_ds2_ddg <- -cf * e2iwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2 #
  df_ds2_ddg[is.na(df_ds2_ddg)] <- 0

  df_ds1_dgwt <- cf * (e1wt * eiwt * (1 + eiwt)  - e1iwt * ewt * (1 + ewt)) / 
    ((1 + eiwt) * ewt)^2 #

  df_ds2_dgwt <- cf * (e2wt * eiwt * (1 + eiwt)  - e2iwt * ewt * (1 + ewt)) /
    ((1 + eiwt) * ewt)^2 #

  df_ds_fwt <- eiwt * (1 + ewt) / ((1 + eiwt) * ewt) #

  df_ds_f0 <- (ewt - eiwt) / ((1 + eiwt) * ewt) #

  f_pred <- function_dG_4state_folding_dg2fitness(
    s1_ddg = s1_ddg, 
    s2_ddg = s2_ddg,
    s1_dgwt = s1_dgwt, 
    s2_dgwt = s2_dgwt, 
    s_fwt = s_fwt,
    s_f0 = s_f0
  )

  derr_df <- -2 * (f - f_pred) / w^2
  derr_df[is.na(derr_df)] <- 0

  gradient_s1_ddg <- mutXvar %*% (derr_df * df_ds1_ddg)
  gradient_s2_ddg <- mutXvar %*% (derr_df * df_ds2_ddg)
  gradient_s1_dgwt <- sum(derr_df * df_ds1_dgwt, na.rm = T)
  gradient_s2_dgwt <- sum(derr_df * df_ds2_dgwt, na.rm = T)
  gradient_s_fwt <- sum(derr_df * df_ds_fwt, na.rm = T)
  gradient_s_f0 <- sum(derr_df * df_ds_f0, na.rm = T)
  
  return(list(
    s1_ddg = gradient_s1_ddg,
    s2_ddg = gradient_s2_ddg,
    s1_dgwt = gradient_s1_dgwt,
    s2_dgwt = gradient_s2_dgwt,
    s_fwt = gradient_s_fwt,
    s_f0 = gradient_s_f0
  ))
}

function_dG_4state_binding_dg2fitness <- function(
  b_ddg, 
  s1_ddg, 
  s2_ddg, 
  b_dgwt, 
  s1_dgwt, 
  s2_dgwt, 
  b_fwt,
  b_f0,
  rt = 1.99e-3 * 310.15
) {
  bf <- b_f0 + 
    (b_fwt - b_f0) * (1 + exp(b_dgwt / rt) * (1 + exp(s2_dgwt / rt) * (1 + exp(-s1_dgwt / rt)))) /
     (1 + exp((b_dgwt + b_ddg) / rt) * 
        (1 + exp((s2_dgwt + s2_ddg) / rt) * (1 + exp(-(s1_dgwt + s1_ddg) / rt))))
}

function_dG_4state_binding_dg2fitness_gradient <- function(
  s1_ddg,  
  s2_ddg,  
  b_ddg,
  s1_dgwt,
  s2_dgwt,
  b_dgwt,
  b_fwt,
  b_f0, 
  f,
  w,
  mutXvar,
  rt = 1.99e-3 * 310.15
) {
  cf <- (b_fwt - b_f0) / rt

  ebi <- exp(b_ddg / rt)
  ebwt <- exp(b_dgwt / rt)
  ebiwt <- ebi * ebwt

  e1si <- exp(-s1_ddg / rt)
  e1swt <- exp(-s1_dgwt / rt)
  e1siwt <- e1si * e1swt
  
  e2si <- exp(s2_ddg / rt)
  e2swt <- exp(s2_dgwt / rt)
  e2siwt <- e2si * e2swt
  
  eiwt <- ebiwt * (1 + e2siwt * (1 + e1siwt))
  ewt <- ebwt * (1 + e2swt * (1 + e1swt))

  df_db_ddg <- -cf * (1 + ewt) * eiwt /
    (1 + eiwt) ^ 2 #
  df_db_ddg[is.na(df_db_ddg)] <- 0

  df_ds1_ddg <- cf * ebiwt * e2siwt * e1siwt * (1 + ewt) /
    (1 + eiwt) ^ 2 #
  df_ds1_ddg[is.na(df_ds1_ddg)] <- 0
  
  df_ds2_ddg <- -cf * ebiwt * e2siwt * (1 + e1siwt) * (1 + ewt) /
    (1 + eiwt) ^ 2 #
  df_ds2_ddg[is.na(df_ds2_ddg)] <- 0

  df_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!

  df_ds1_dgwt <- cf * (-ebwt * e2swt * e1swt * (1 + eiwt) + (1 + ewt) * ebiwt * e2siwt * e1siwt) / 
    (1 + eiwt) ^ 2 #

  df_ds2_dgwt <- cf * (ebwt * e2swt * (1 + e1swt) * (1 + eiwt) - 
    (1 + ewt) * ebiwt * e2siwt * (1 + e1siwt)) / 
    (1 + eiwt) ^ 2 #
    
  df_db_fwt <- (1 + ewt) / (1 + eiwt) #
  df_db_f0 <- 1 - (1 + ewt) / (1 + eiwt) #

  f_pred <- function_dG_4state_binding_dg2fitness(
    b_ddg = b_ddg,
    s1_ddg = s1_ddg,
    s2_ddg = s2_ddg,
    b_dgwt = b_dgwt,
    s1_dgwt = s1_dgwt,  
    s2_dgwt = s2_dgwt,  
    b_fwt = b_fwt,
    b_f0 = b_f0
  )

  derr_df <- -2 * (f - f_pred) / w^2
  derr_df[is.na(derr_df)] <- 0
  gradient_b_ddg <- mutXvar %*% (derr_df * df_db_ddg)
  gradient_s1_ddg <- mutXvar %*% (derr_df * df_ds1_ddg)
  gradient_s2_ddg <- mutXvar %*% (derr_df * df_ds2_ddg)
  
  gradient_b_dgwt <- sum(derr_df * df_db_dgwt, na.rm = T)
  gradient_s1_dgwt <- sum(derr_df * df_ds1_dgwt, na.rm = T)
  gradient_s2_dgwt <- sum(derr_df * df_ds2_dgwt, na.rm = T)
  
  gradient_b_fwt <- sum(derr_df * df_db_fwt, na.rm = T)
  gradient_b_f0 <- sum(derr_df * df_db_f0, na.rm = T)
  
  return(list(
    b_ddg = gradient_b_ddg, 
    s1_ddg = gradient_s1_ddg, 
    s2_ddg = gradient_s2_ddg,
    b_dgwt = gradient_b_dgwt, 
    s1_dgwt = gradient_s1_dgwt,
    s2_dgwt = gradient_s2_dgwt,
    b_fwt = gradient_b_fwt, 
    b_f0 = gradient_b_f0
  ))
}
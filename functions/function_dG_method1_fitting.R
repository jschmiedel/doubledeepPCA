############################################################
#### fit global relationship between deltaGs and fitess
#### from binding and stabilityPCA data assuming
#### deltadeltaG binding = 0 for all used variants
############################################################
function_dG_method1_fitting <- function(
  parameters, 
  id_L, 
  global_par, 
  list_fs
) {

  # # fit stability delta Gs
  s_dG <- parameters[1:id_L]
  if (length(parameters) == (id_L + 6)) {
    s_scale <- parameters[id_L + 1]
    b_scale <- parameters[id_L + 2]
    s_bgr <- parameters[id_L + 3]
    b_bgr <- parameters[id_L + 4]
    s_dGwt <- parameters[id_L + 5]
    b_dGwt <- parameters[id_L + 6]
    b_dG <- 0
  } else if (length(parameters) == 2 * id_L) {
    # + binding delta Gs [>> determine dGs after global relationship was already determined]
    b_dG <- parameters[(id_L + 1):(2 * id_L)]
    s_scale <- global_par[1]
    b_scale <- global_par[2]
    s_bgr <- global_par[3]
    b_bgr <- global_par[4]
    s_dGwt <- global_par[5]
    b_dGwt <- global_par[6]
  }


  sf <- function_folding_dG2F(
    s_dG = s_dG,
    s_dGwt = s_dGwt,
    s_bgr = s_bgr,
    s_scale = s_scale
  )
  bf <- function_binding_dG2F(
    b_dG = b_dG,
    s_dG = s_dG,
    b_dGwt = b_dGwt,
    s_dGwt = s_dGwt,
    b_bgr = b_bgr,
    b_scale = b_scale
  )

  # mean square deviation; divide by weights to correct for variants being NA
  MSD <- sum((list_fs$b_fitness - bf)^2 / list_fs$b_sigma^2 +
    (list_fs$s_fitness - sf)^2 / list_fs$s_sigma^2, na.rm = T) /
    sum(list_fs$b_sigma^-2 + list_fs$s_sigma^-2, na.rm = T)

  return(MSD)
}
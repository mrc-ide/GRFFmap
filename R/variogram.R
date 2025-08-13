
#------------------------------------------------
#' @title TODO
#'   
#' @description Get the correlation from a spatial-temporal variogram model.
#' TODO - more text.
#'   
#' @param x TODO
#'
#' @export

get_vario_cor <- function(dist_space,
                          dist_time,
                          length_space,
                          length_time,
                          kernel_space_power = 1,
                          kernel_time_power = 2) {
  
  exp(-(dist_space / length_space)^kernel_space_power) * exp(-(dist_time / length_time)^kernel_time_power)
}

#------------------------------------------------
#' @title TODO
#'   
#' @description Get the likelihood from a spatial-temporal variogram model.
#' TODO - more text.
#'   
#' @param x TODO
#'
#' @export

get_vario_loglike <- function(dist_space,
                              dist_time,
                              dist_z,
                              length_space,
                              length_time,
                              nugget,
                              partial_sill,
                              kernel_space_power = 1,
                              kernel_time_power = 2) {
  
  rho <- get_vario_cor(dist_space = dist_space,
                       dist_time = dist_time,
                       length_space = length_space,
                       length_time = length_time,
                       kernel_space_power = kernel_space_power,
                       kernel_time_power = kernel_time_power)
  c <- nugget + partial_sill*(1 - rho)
  ll <- dnorm(dist_z, sd = sqrt(c), log = TRUE)
  sum(ll)
}

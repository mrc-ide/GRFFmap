
#------------------------------------------------
#' @title Empirical logit function
#'   
#' @description The empirical logit is a version of the logit transform that is
#'   shrunk towards a value, in our case 0.5. This avoids issues with infinite
#'   numbers when prevalence is 0% or 100%.
#'   
#' @param x TODO
#'
#' @export

emplogit <- function(x, N) {
  log(x + 0.5) - log(N - x + 0.5)
}

#------------------------------------------------
#' @title TODO
#'   
#' @description TODO
#'   
#' @param x a vector.
#' @param n the number of values to expand in each direction.
#'
#' @export

buffer_vector <- function(x, n) {
  
  # check inputs
  assert_vector_numeric(x)
  assert_increasing(x)
  assert_single_pos_int(n, zero_allowed = TRUE)
  
  dx <- x[2] - x[1]
  c(x[1] - rev(seq_len(n))*dx, x, x[length(x)] + seq_len(n)*dx)
}

#------------------------------------------------
#' @title Simulate values on a space-time grid by drawing from a GRF using spectral methods
#'   
#' @description simulate values on a space-time grid. The range in each
#'   dimension is set by the user along with the number of grid cells.
#'   Simulation is via a spectral approximation, which is efficient for large
#'   grid sizes.
#'   
#' @param x TODO.
#'
#' @export

sim_spacetime_spectral <- function(x_range,
                                   y_range,
                                   t_range,
                                   nx,
                                   ny,
                                   x_buffer = 10,
                                   y_buffer = 10,
                                   t_buffer = 10,
                                   length_space,
                                   length_time,
                                   kernel_space_type = "RBF",
                                   kernel_time_type = "RBF",
                                   sill = 1,
                                   mu = 0) {
  
  # check inputs
  assert_limit(x_range)
  assert_limit(y_range)
  assert_limit(t_range)
  assert_int(t_range)
  assert_single_pos_int(nx, zero_allowed = FALSE)
  assert_single_pos_int(ny, zero_allowed = FALSE)
  assert_single_pos_int(x_buffer, zero_allowed = TRUE)
  assert_single_pos_int(y_buffer, zero_allowed = TRUE)
  assert_single_pos_int(t_buffer, zero_allowed = TRUE)
  assert_single_pos(length_space, zero_allowed = FALSE)
  assert_single_pos(length_time, zero_allowed = FALSE)
  assert_in(kernel_space_type, c("RBF", "Exp"))
  assert_in(kernel_time_type, c("RBF", "Exp"))
  assert_single_pos(sill)
  assert_single_numeric(mu)
  
  # make simulation grid in space and time. The values in the spatial grid
  # represent the midpoints of cells, and are chosen so that the extreme edges
  # of the cells span the desired range. The values in the time grid are
  # integers spanning the stated range (inclusive)
  dx <- diff(x_range) / nx
  dy <- diff(y_range) / ny
  x_vec <- seq(x_range[1] + dx/2, x_range[2] - dx/2, by = dx)
  y_vec <- seq(y_range[1] + dy/2, y_range[2] - dy/2, by = dy)
  t_vec <- t_range[1]:t_range[2]
  
  # buffer the vectors. This is needed because the spectral method returns a
  # periodic solution, and so to avoid this we run with a buffer and then chop
  # it off later
  x_vec_buffer <- buffer_vector(x_vec, x_buffer)
  y_vec_buffer <- buffer_vector(y_vec, y_buffer)
  t_vec_buffer <- buffer_vector(t_vec, t_buffer)
  nx_buffer <- nx + 2*x_buffer
  ny_buffer <- ny + 2*y_buffer
  nt <- diff(t_range) + 1
  nt_buffer <- nt + 2*t_buffer
  
  # define white noise array. All array dimensions are chosen to go x then y then t
  m1_buffer <- sqrt(sill) * array(data = rnorm(nx_buffer * ny_buffer * nt_buffer),
                     dim = c(nx_buffer, ny_buffer, nt_buffer))
  
  # define kernel arrays
  m2_x_buffer <- array(x_vec_buffer - mean(x_range), dim = c(nx_buffer, ny_buffer, nt_buffer))
  m2_y_buffer <- aperm(array(y_vec_buffer - mean(y_range), dim = c(ny_buffer, nx_buffer, nt_buffer)), perm = c(2, 1, 3))
  m2_t_buffer <- aperm(array(t_vec_buffer - mean(t_range), dim = c(nt_buffer, ny_buffer, nx_buffer)), perm = c(3, 2, 1))
  
  # calculate distance in space and time and compute kernel
  m2_dist_space <- sqrt(m2_x_buffer^2 + m2_y_buffer^2)
  m2_dist_time <- m2_t_buffer
  if (kernel_space_type == "RBF") {
    m2_z_buffer <- exp(-0.5*(m2_dist_space / length_space)^2)
  } else if (kernel_space_type == "Exp") {
    m2_z_buffer <- exp(-m2_dist_space / length_space)
  }
  if (kernel_time_type == "RBF") {
    m2_z_buffer <- m2_z_buffer * exp(-0.5*(m2_dist_time / length_time)^2)
  } else if (kernel_time_type == "Exp") {
    m2_z_buffer <- m2_z_buffer * exp(-m2_dist_time / length_time)
  }
  
  # perform convolution
  m3_buffer <- Re(fft(sqrt(fft(m2_z_buffer)) * fft(m1_buffer), inverse = TRUE)) / length(m1_buffer)
  
  # discard buffer
  m3 <- m3_buffer[1:nx + x_buffer, 1:ny + y_buffer, 1:nt + t_buffer, drop = FALSE]
  
  # apply mean
  m3 <- m3 + mu
  
  # return list
  list(x = x_vec,
       y = y_vec,
       t = t_vec,
       grid = m3)
}

#------------------------------------------------
#' @title Simulate data by sampling without replacement from a space-time grid
#'   
#' @description Given a space-time grid produced by
#'   \code{sim_spacetime_spectral()}, draw sampling locations by replacement.
#'   Then draw positive sample numbers given a defined sample size.
#'   
#' @param x TODO
#'
#' @export

array_draw_data <- function(sim_array, n_data, N) {
  
  # check inputs
  assert_list(sim_array)
  assert_in(c("x", "y", "t", "grid"), names(sim_array))
  assert_single_pos_int(n_data, zero_allowed = FALSE)
  assert_in(length(N), c(1, n_data), message = "N must be either a single value, or a vector of length n_data")
  assert_pos_int(N, zero_allowed = FALSE)
  if (length(N) == 1) {
    N <- rep(N, n_data)
  }
  
  # sample indices at random from the grid
  dims <- dim(sim_array$grid)
  indices <- expand_grid(1:dims[1], 1:dims[2], 1:dims[3]) |>
    sample_n(n_data) |>
    as.matrix()
  
  # extract values at these indices, and draw from the binomial distribution
  df_dat <- data.frame(x = sim_array$x[indices[,1]],
                       y = sim_array$y[indices[,2]],
                       t = sim_array$t[indices[,3]],
                       z = sim_array$grid[as.matrix(indices)]) |>
    mutate(p = plogis(z),
           N = N,
           k = rbinom(n(), size = N, prob = p),
           p_est = k / N,
           p_emplogit = emplogit(k, N))
  
  return(df_dat)
}

#------------------------------------------------
#' @title TODO
#'   
#' @description TODO
#'   
#' @param x TODO
#'
#' @export

get_data_dist <- function(x, y, t, z) {
  
  # check inputs
  assert_vector_numeric(x)
  assert_vector_numeric(y)
  assert_vector_pos_int(t)
  assert_vector_numeric(z)
  assert_same_length_multiple(x, y, t, z)
  n <- length(x)
  
  # calculate distances
  dist_space <- dist(cbind(x, y))
  dist_time <- dist(t)
  dist_z_mat <- outer(z, z, FUN = "-")
  dist_z <- dist_z_mat[lower.tri(dist_z_mat)]
  
  # make data.frame
  ret <- which(lower.tri(dist_space), arr.ind = TRUE) |>
    as.data.frame() |>
    setNames(c("ind1", "ind2")) |>
    mutate(dist_space = as.vector(dist_space),
           dist_time = as.vector(dist_time),
           dist_z = as.vector(dist_z))
  
  ret
}

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

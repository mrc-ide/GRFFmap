
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
                                   nt,
                                   x_buffer = 10,
                                   y_buffer = 10,
                                   t_buffer = 10,
                                   length_space = 1,
                                   length_time = 1,
                                   kernel_space_type = "RBF",
                                   kernel_time_type = "RBF",
                                   mu = 0,
                                   sigma = 1) {
  
  # check inputs
  assert_limit(x_range)
  assert_limit(y_range)
  assert_limit(t_range)
  assert_single_pos_int(nx, zero_allowed = FALSE)
  assert_single_pos_int(ny, zero_allowed = FALSE)
  assert_single_pos_int(nt, zero_allowed = FALSE)
  assert_single_pos_int(x_buffer, zero_allowed = TRUE)
  assert_single_pos_int(y_buffer, zero_allowed = TRUE)
  assert_single_pos_int(t_buffer, zero_allowed = TRUE)
  assert_single_pos(length_space, zero_allowed = FALSE)
  assert_single_pos(length_time, zero_allowed = FALSE)
  assert_in(kernel_space_type, c("RBF", "Exp"))
  assert_in(kernel_time_type, c("RBF", "Exp"))
  assert_single_numeric(mu)
  assert_single_pos(sigma)
  
  # Grid spacings
  dx <- diff(x_range) / nx
  dy <- diff(y_range) / ny
  dt <- diff(t_range) / nt
  
  # Total number of cells in each dimension including buffer
  px <- nx + x_buffer
  py <- ny + y_buffer
  pt <- nt + t_buffer
  
  # Embedding vectors: lags increase then wrap to negatives
  hx <- c(0:(px / 2 - 1), (-px / 2):-1) * dx
  hy <- c(0:(py / 2 - 1), (-py / 2):-1) * dy
  ht <- c(0:(pt / 2 - 1), (-pt / 2):-1) * dt
  
  # Create 3D coordinate grids
  Hx <- array(hx, dim = c(px, py, pt))
  Hy <- aperm(array(hy, dim = c(py, px, pt)), perm = c(2, 1, 3))
  Ht <- aperm(array(ht, dim = c(pt, py, px)), perm = c(3, 2, 1))
  
  # Apply kernel
  if (kernel_space_type == "RBF") {
    kernel_array <- sigma^2 * exp(- (Hx^2 + Hy^2) / (2 * length_space^2))
  } else if (kernel_space_type == "Exp") {
    kernel_array <- sigma^2 * exp(- sqrt(Hx^2 + Hy^2) / length_space)
  }
  if (kernel_time_type == "RBF") {
    kernel_array <- kernel_array * exp(- Ht^2 / (2 * length_time^2))
  } else if (kernel_time_type == "Exp") {
    kernel_array <- kernel_array * exp(- sqrt(Ht^2) / length_time)
  }
  
  # FFT of the kernel
  fft_kernel <- fft(kernel_array)
  
  # Simulate white noise
  noise <- array(rnorm(px * py * pt), dim = c(px, py, pt))
  
  # Multiply in Fourier space
  sqrt_fft_kernel <- sqrt(pmax(Re(fft_kernel), 0))
  fft_noise <- fft(noise)
  fft_field <- sqrt_fft_kernel * fft_noise
  
  # Inverse FFT to real space
  field <- Re(fft(fft_field, inverse = TRUE)) / (px * py * pt)
  
  # Crop to desired size and add mean
  field_out <- field[1:nx, 1:ny, 1:nt, drop = FALSE] + mu
  
  # Coordinate vectors (cell centres)
  x_vals <- seq(from = x_range[1] + dx / 2, to = x_range[2] - dx / 2, length.out = nx)
  y_vals <- seq(from = y_range[1] + dy / 2, to = y_range[2] - dy / 2, length.out = ny)
  t_vals <- seq(from = t_range[1] + dt / 2, to = t_range[2] - dt / 2, length.out = nt)
  
  return(list(
    x = x_vals,
    y = y_vals,
    t = t_vals,
    grid = field_out
  ))
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
#' @importFrom dplyr sample_n
#' @export

array_draw_data <- function(sim_array, n_data, n_samp) {
  
  # check inputs
  assert_list(sim_array)
  assert_in(c("x", "y", "t", "grid"), names(sim_array))
  assert_single_pos_int(n_data, zero_allowed = FALSE)
  assert_in(length(n_samp), c(1, n_data), message = "n_samp must be either a single value, or a vector of length n_data")
  assert_pos_int(n_samp, zero_allowed = FALSE)
  if (length(n_samp) == 1) {
    n_samp <- rep(n_samp, n_data)
  }
  
  # sample indices at random from the grid
  dims <- dim(sim_array$grid)
  indices <- data.frame(x = sample.int(dims[1], n_data, replace = TRUE),
                        y = sample.int(dims[2], n_data, replace = TRUE),
                        t = sample.int(dims[3], n_data, replace = TRUE))
  
  # extract values at these indices, and draw from the binomial distribution
  df_dat <- data.frame(x = sim_array$x[indices[,1]],
                       y = sim_array$y[indices[,2]],
                       t = sim_array$t[indices[,3]],
                       z = sim_array$grid[as.matrix(indices)]) |>
    mutate(p = plogis(z),
           n_samp = n_samp,
           n_pos = rbinom(n(), size = n_samp, prob = p),
           p_est = n_pos / n_samp,
           p_emplogit = emplogit(n_pos, n_samp))
  
  return(df_dat)
}

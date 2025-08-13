# 17.circulant_fix.R
#
# Author: Bob Verity
# Date: 2025-07-30
#
# Purpose:
# Finally found and fixed a bug that was making it look like RFF variance
# estimation was way off. Turned out to be largely caused by a bug in my
# simulator, meaning the simulated variance was not what I thought it was. This
# was caused by the way the "circulant embedding" spectral method generated the
# kernel matrix, along with discretization correction factors and other fiddly
# things. These are now fixed, variances make sense, estimation of sigma via the
# true likelihood (no approximation) gives the correct values, and RFFs are
# looking good again.
#
# ------------------------------------------------------------------

library(tidyverse)
library(plotly)
library(rstan)

get_unscaled_rbf_kernel <- function(x1, y1, t1, x2, y2, t2, length_space, length_time) {
  D2_space <- (outer(x1, x2, "-"))^2 + (outer(y1, y2, "-"))^2
  D2_time <- (outer(t1, t2, "-"))^2
  K <- exp(-D2_space / (2 * length_space^2) - D2_time / (2 * length_time^2))
  return(K)
}

sim_spacetime_spectral_circulant <- function(x_range, y_range, t_range,
                                             nx, ny, nt,
                                             x_buffer = 10, y_buffer = 10, t_buffer = 10,
                                             length_space, length_time,
                                             mu = 0, sigma = 1) {
  
  # Grid spacings
  dx <- diff(x_range) / nx
  dy <- diff(y_range) / ny
  dt <- diff(t_range) / nt
  
  # Total size including buffer
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
  
  # Compute squared distances
  spatial_sq <- (Hx^2 + Hy^2) / (2 * length_space^2)
  temporal_sq <- (Ht^2) / (2 * length_time^2)
  
  # Apply RBF kernel
  kernel_array <- sigma^2 * exp(- (spatial_sq + temporal_sq))
  
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
  
  # Coordinate vectors (cell centers)
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

# ------------------------------------------------------------------------------------------------
# SIMULATE DATA

#set.seed(3)

# parameters
x_range <- c(-180, 180)
y_range <- c(-90, 90)
t_range <- c(0, 100)
nx <- 120
ny <- 80
nt <- 100
x_buffer <- 100
y_buffer <- 100
t_buffer <- 100
length_space <- 30
length_time <- 30
sigma <- 2
mu <- -3

# simulate z values in space-time grid
sim_grid <- sim_spacetime_spectral_circulant(x_range = x_range,
                                             y_range = y_range,
                                             t_range = t_range,
                                             nx = nx,
                                             ny = ny,
                                             nt = nt,
                                             x_buffer = x_buffer,
                                             y_buffer = y_buffer,
                                             t_buffer = t_buffer,
                                             length_space = length_space,
                                             length_time = length_time,
                                             mu = mu,
                                             sigma = sigma)

# get into long form data.frame, and compute prevalence as logistic(z)
df_sim_grid <- expand_grid(t = sim_grid$t,
                           y = sim_grid$y,
                           x = sim_grid$x) |>
  mutate(z = as.vector(sim_grid$grid),
         p = plogis(z))

mean(df_sim_grid$z)
sd(df_sim_grid$z)

# plot at a few select times
t_plot <- c(1, 30, 60, 100)

df_sim_grid |>
  filter(ceiling(t) %in% t_plot) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = p)) +
  facet_wrap(~t) +
  #scale_fill_viridis_c(option = "magma")
  scale_fill_viridis_c(option = "magma", limits = c(0, 1))

# simulate data from grid
n_data <- 2e2
df_dat <- array_draw_data(sim_array = sim_grid,
                          n_data = n_data,
                          N = 100)

mean(df_dat$z)
sd(df_dat$z)

# ------------------------------------------------
# exact likelihood

if (FALSE) {
  
  df_mu_sigma <- expand_grid(mu = seq(-4, 4, l = 101),
                            sigma = seq(0.1, 6, l = 101),
                            ll = NA)
  
  Sigma <- get_unscaled_rbf_kernel(x1 = df_dat$x, y1 = df_dat$y, t1 = df_dat$t,
                                   x2 = df_dat$x, y2 = df_dat$y, t2 = df_dat$t, 
                                   length_space = length_space, 
                                   length_time = length_time)
  
  for (i in 1:nrow(df_mu_sigma)) {
    df_mu_sigma$ll[i] <- mvtnorm::dmvnorm(x = df_dat$z,
                                          mean = rep(df_mu_sigma$mu[i], n_data), 
                                          sigma = df_mu_sigma$sigma[i]^2 * Sigma, 
                                          log = TRUE)
  }
  
  df_mu_sigma |>
    mutate(like = exp(ll - max(ll))) |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = mu, y = sigma, fill = like)) +
    geom_point(x = mu, y = sigma, col = "red")
  
}

# ------------------------------------------------

# draw features
n_features <- 100
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
omega_time <- rbind(rnorm(n_features))

# make feature maps
X <- df_dat |>
  dplyr::select(x, y) |>
  as.matrix()

feat <- cbind(cos(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time)) * sqrt(1 / n_features)

# prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  feat = feat,
  k = df_dat$k,
  n_trials = df_dat$N,
  mu_mean = 0,
  mu_sd = 5,
  sigma_shape = 1,
  sigma_rate = 0.1
)

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/fixedlength_v4.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

# extract MCMC samples
draws <- extract(fit)
n_draws <- nrow(draws$beta)

beta_samples <- draws$beta
mu_samples <- draws$mu
sigma_samples <- draws$sigma

hist(mu_samples, breaks = 100)
abline(v = mu, col = 2, lwd = 3)

hist(sigma_samples, breaks = 100)
abline(v = sigma, col = 2, lwd = 3)

# ------------------------------------------------

# make feature maps for prediction
X_pred <- expand.grid(x = sim_grid$x,
                      y = sim_grid$y) |>
  as.matrix()

l <- list()
for (i in seq_along(t_plot)) {
  message(sprintf("%s of %s", i, length(t_plot)))
  
  # fix prediction at this time
  t_pred <- rep(t_plot[i], nrow(X_pred))
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    sweep(MARGIN = 2, STATS = mu_samples, FUN = "+") |>
    plogis()
  
  # summarise posterior predictions
  pred_quant <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
    t() |>
    as.data.frame()
  
  # get combined data.frame from inferred distribution and truth
  df_pred <- data.frame(x = X_pred[,1],
                        y = X_pred[,2],
                        t = t_pred) |>
    bind_cols(pred_quant) |>
    mutate(MOE = `97.5%` - `2.5%`,
           p_pred = `50%`,
           #mask = (MOE > 0.3), # mask based on total MOE
           #mask = (`97.5%`/p_pred > 2) | (`2.5%`/p_pred < 0.5), # mask based on multiplicative prevalence
           mask = (`97.5%` - p_pred) > 0.1 | (p_pred - `2.5%`) > 0.1, # mask based on absolute prevalence
           p_true = df_sim_grid |>
             filter(ceiling(t) == t_plot[i]) |>
             pull(p))
  
  l[[i]] <- df_pred
}
df_combined <- bind_rows(l)

# produce plots
df_masked <- df_combined |>
  filter(mask == TRUE) |>
  mutate(Type = factor("p_pred", levels = c("p_true", "p_pred")))
dx <- diff(sim_grid$x)[1]
dy <- diff(sim_grid$y)[1]

t_eps <- 5
df_dat_plot <- df_dat |>
  expand_grid(t_plot = t_plot) |>
  filter(t > (t_plot - t_eps) & t < (t_plot + t_eps)) |>
  mutate(t = t_plot)

df_combined |>
  dplyr::select(x, y, t, p_pred, p_true) |>
  pivot_longer(cols = c(p_pred, p_true), names_to = "Type") |>
  mutate(Type = factor(Type, levels = c("p_true", "p_pred"))) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = value)) +
  #geom_tile(aes(x = x, y = y), width = dx, height = dy, fill = grey(0.5), alpha = 0.5, data = df_masked) +
  geom_point(aes(x = x, y = y, fill = p_est), pch = 21, colour = grey(0.5), size = 2, data = df_dat_plot) +
  facet_grid(Type ~ t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  ggtitle("Inferred vs. True Prevalence")


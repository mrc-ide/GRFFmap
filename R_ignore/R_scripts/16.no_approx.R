# 16.no_approx.R
#
# Author: Bob Verity
# Date: 2025-07-30
#
#
# Purpose:
# Experiment with using exact solution, no RFF approximation. Uses full GRF likelihood in STAN code.
# Conclusion: unlike RFF, this method works well for hyper-parameter estimation (e.g. sill). However,
# it is very slow. Only really viable for 100-200 sites, going to 1000 would be way too many.
#
# ------------------------------------------------------------------

library(tidyverse)
library(plotly)
library(rstan)

get_unscaled_rbf_kernel <- function(x1, y1, t1, x2, y2, t2, length_space, length_time) {
  D2_space <- (outer(x1, x2, "-"))^2 + (outer(y1, y2, "-"))^2
  D2_time <- (outer(t1, t2, "-"))^2
  K <- exp(-D2_space / (2 * length_space^2)) *
    exp(-D2_time / (2 * length_time^2))
  return(K)
}

# ------------------------------------------------------------------------------------------------
# SIMULATE DATA

set.seed(2)

# parameters
x_range <- c(-180, 180)
y_range <- c(-90, 90)
t_range <- c(1, 100)
nx <- 120
ny <- 80
length_space <- 30
length_time <- 30
sigma <- 2
mu <- -4

# simulate z values in space-time grid
sim_grid <- sim_spacetime_spectral(x_range = x_range,
                                   y_range = y_range,
                                   t_range = t_range,
                                   nx = nx,
                                   ny = ny,
                                   length_space = length_space,
                                   length_time = length_time,
                                   x_buffer = 50,
                                   y_buffer = 50,
                                   t_buffer = 50,
                                   sill = sigma^2,
                                   mu = mu)

# get into long form data.frame, and compute prevalence as logistic(z)
df_sim_grid <- expand_grid(t = sim_grid$t,
                           y = sim_grid$y,
                           x = sim_grid$x) |>
  mutate(z = as.vector(sim_grid$grid),
         p = plogis(z))

# plot at a few select times
t_plot <- c(1, 30, 60, 100)

df_sim_grid |>
  filter(t %in% t_plot) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = p)) +
  facet_wrap(~t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1))

# simulate data from grid
n_data <- 1e2
df_dat <- array_draw_data(sim_array = sim_grid,
                          n_data = n_data,
                          N = 1e2)

# ------------------------------------------------

K_oo <- sill * get_unscaled_rbf_kernel(x1 = df_dat$x, y1 = df_dat$y, t1 = df_dat$t,
                                       x2 = df_dat$x, y2 = df_dat$y, t2 = df_dat$t, 
                                       length_space = length_space, length_time = length_time)

# prepare data list for Stan
data_list <- list(
  N = n_data,
  k = df_dat$k,
  n = df_dat$N,
  K = K_oo
)

# run stan model
t0 <- Sys.time()
fit <- stan(
  file = "R_ignore/R_scripts/stan/no_approx.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)
Sys.time() - t0

# extract MCMC samples
draws <- extract(fit)
n_draws <- nrow(draws$beta)

z_samples <- draws$z
mu_samples <- draws$mu
sigma_samples <- draws$sigma

hist(mu_samples, breaks = 100)
abline(v = mu, lwd = 2, col = 2)

hist(sigma_samples, breaks = 100)
abline(v = sigma, lwd = 2, col = 2)

# ------------------------------------------------

# Compute covariance matrix between observations and predictions
K_po <- sill * get_unscaled_rbf_kernel(x1 = df_dat$x, y1 = df_dat$y, t1 = df_dat$t,
                                       x2 = df_sim_grid$x, y2 = df_sim_grid$y, t2 = df_sim_grid$t, 
                                       length_space = length_space, length_time = length_time)

K_oo_inv <- solve(K_oo)
tmp <- t(K_po) %*% K_oo_inv

#i <- 1

mu_pred <- tmp %*% colMeans(z_samples)
mu_pred <- mean(mu_samples) + mean(sigma_samples)*mu_pred

t_eps <- 5
df_dat_plot <- df_dat |>
  expand_grid(t_plot = t_plot) |>
  filter(t > (t_plot - t_eps) & t < (t_plot + t_eps)) |>
  mutate(t = t_plot)

df_sim_grid |>
  mutate(p_pred = plogis(mu_pred)) |>
  filter(t %in% t_plot) |>
  pivot_longer(cols = c(p, p_pred), names_to = "Type") |>
  mutate(Type = factor(Type, levels = c("p", "p_pred"))) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = value)) +
  geom_point(aes(x = x, y = y, fill = p_est), pch = 21, colour = grey(0.5), size = 2, data = df_dat_plot) +
  facet_grid(Type ~ t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1))

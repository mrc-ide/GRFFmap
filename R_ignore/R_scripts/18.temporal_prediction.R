# 18.temporal_prediction.R
#
# Author: Bob Verity
# Date: 2025-08-13
#
# Purpose:
# First attempt at moving beyond static time-slice predictions to true temporal
# prediction. After fitting the model to spatially sampled data, posterior draws
# are projected forward (and backward) in time to recover the changing regional
# average prevalence. Results now show full temporal trajectories with credible
# intervals, rather than just isolated snapshots.
#
# ------------------------------------------------------------------

library(tidyverse)
library(plotly)
library(rstan)
library(epitools)

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
length_space <- 50
length_time <- 30
sigma <- 2
mu <- -2

# simulate z values in space-time grid
sim_grid <- sim_spacetime_spectral(x_range = x_range,
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

# plot at a few select times
t_plot <- c(1, 30, 60, 100)

df_sim_grid |>
  filter(ceiling(t) %in% t_plot) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = p)) +
  facet_wrap(~t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1))

# imagine our region of interest is a circle of radius 50 centred at [0,0].
# Plot just this region
df_sim_grid |>
  filter((x^2 + y^2) <= 50^2) |>
  filter(ceiling(t) %in% t_plot) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = p)) +
  facet_wrap(~t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  xlim(c(-100, 100)) + ylim(c(-100, 100))

# simulate data by drawing from grid
n_data <- 1e2
df_dat <- array_draw_data(sim_array = sim_grid,
                          n_data = n_data,
                          n_samp = 100)

# plot true average prevalence in the region over time. Overlay raw data CIs for
# sites within this region
df_sim_time <- df_sim_grid |>
  filter((x^2 + y^2) <= 50^2) |>
  group_by(t) |>
  summarise(p = mean(p))

df_dat_time <- df_dat |>
  filter((x^2 + y^2) <= 50^2) |>
  mutate(lower = epitools::binom.exact(x = n_pos, n = n_samp)$lower,
         upper = epitools::binom.exact(x = n_pos, n = n_samp)$upper)
  
df_sim_time |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p)) +
  geom_pointrange(aes(x = t, y = p_est, ymin = lower, ymax = upper), data = df_dat_time) +
  ylim(c(0, 1)) +
  ggtitle("Regional average prevalence (line)")


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
  k = df_dat$n_pos,
  n_trials = df_dat$n_samp,
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

# ------------------------------------------------

# make feature maps for prediction
X_pred <- expand.grid(x = sim_grid$x,
                      y = sim_grid$y) |>
  as.matrix()

t_pred_vec <- df_sim_time$t

l <- list()
for (i in seq_along(t_pred_vec)) {
  message(sprintf("%s of %s", i, length(t_pred_vec)))
  
  # fix prediction at this time
  t_pred <- rep(t_pred_vec[i], nrow(X_pred))
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    sweep(MARGIN = 2, STATS = mu_samples, FUN = "+") |>
    plogis()
  
  # subset to grid cells in the region, and take the mean. Do this *for each posterior draw*
  w <- which((X_pred[,1]^2 + X_pred[,2]^2) < 50^2)
  l[[i]] <- data.frame(t = t_pred_vec[i],
                       post_draw = 1:n_draws,
                       p = colMeans(post_all[w,]))
  
}
df_combined <- bind_rows(l)

# now we can calculate quantiles over posterior draws
df_quantiles <- df_combined |>
  group_by(t) |>
  summarise(Q2.5 = quantile(p, probs = 0.025),
            Q50 = quantile(p, probs = 0.5),
            Q97.5 = quantile(p, probs = 0.975))

# plot quantiles
df_quantiles |>
  ggplot() + theme_bw() +
  geom_ribbon(aes(x = t, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.5) +
  geom_line(aes(x = t, y = Q50), col = "red") +
  geom_line(aes(x = t, y = p), data = df_sim_time) +
  geom_pointrange(aes(x = t, y = p_est, ymin = lower, ymax = upper), data = df_dat_time) +
  ylim(c(0, 1))

# spaghetti plot
df_combined |>
  filter(post_draw < 100) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p, group = post_draw), col = "red", alpha = 0.2) +
  geom_line(aes(x = t, y = p), data = df_sim_time) +
  geom_pointrange(aes(x = t, y = p_est, ymin = lower, ymax = upper), data = df_dat_time) +
  ylim(c(0, 1))

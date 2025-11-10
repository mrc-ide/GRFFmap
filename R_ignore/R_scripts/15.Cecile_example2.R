
library(tidyverse)
library(plotly)
library(rstan)

# ------------------------------------------------------------------------------------------------
# SIMULATE DATA

#set.seed(3)

# parameters
x_range <- c(-180, 180)
y_range <- c(-90, 90)
t_range <- c(1, 100)
nx <- 120
ny <- 80
length_space <- 50
length_time <- 100
sill <- 2
mu <- -2

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
                                   sill = sill,
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
n_data <- 200
df_dat <- array_draw_data(sim_array = sim_grid,
                          n_data = n_data,
                          N = 50)

# get distance between observations in space, time, and z-value
data_dist <- get_data_dist(x = df_dat$x,
                           y = df_dat$y,
                           t = df_dat$t,
                           z = df_dat$p_emplogit)

# ------------------------------------------------

# draw features
n_features <- 50
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
omega_time <- rbind(rnorm(n_features))

# make feature maps
X <- df_dat |>
  select(x, y) |>
  as.matrix()

feat <- cbind(cos(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time)) / sqrt(n_features)

# prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  feat = feat,
  k = df_dat$k,
  n_trials = df_dat$N,
  sill_shape = 1,
  sill_rate = 0.1,
  mu_mean = 0,
  mu_sd = 3
)

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/fixedlength_v3.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

# extract MCMC samples
draws <- extract(fit)
n_draws <- nrow(draws$beta)

beta_samples <- draws$beta
sill_samples <- draws$sill
mu_samples <- draws$mu

hist(sill_samples, breaks = 100)
abline(v = sill, lwd = 2, col = 2)

hist(mu_samples, breaks = 100)
abline(v = mu, lwd = 2, col = 2)

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
             filter(t == t_plot[i]) |>
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
  select(x, y, t, p_pred, p_true) |>
  pivot_longer(cols = c(p_pred, p_true), names_to = "Type") |>
  mutate(Type = factor(Type, levels = c("p_true", "p_pred"))) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = value)) +
  geom_tile(aes(x = x, y = y), width = dx, height = dy, fill = grey(0.5), alpha = 0.5, data = df_masked) +
  geom_point(aes(x = x, y = y, fill = p_est), pch = 21, colour = grey(0.5), size = 2, data = df_dat_plot) +
  facet_grid(Type ~ t) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  ggtitle("Inferred vs. True Prevalence")



# 09.sphere_estsigma.R
#
# Author: Bob Verity
# Date: 2024-07-26
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# Fits GRF via RFF on surface of a sphere using geodesic distances rather than raw coordinates.
# Allows for rbf or exponential kernel through mag object.
# Problems when estimating the length_scale - sometimes the posterior appears to
# be multimodal, therefore stan struggles. Drjacoby struggles for other reasons,
# because of the level of correlelation.
# Proposed solution - do not estimate length_scale like this. Instead, fix and
# the MCMC should become much more efficient. Could even consider running the
# MCMC separately over a range of length_scale values and doing some sort of
# model comparison.
#
# ------------------------------------------------------------------

library(geosphere)
library(tidyverse)
library(rstanarm)
library(rstan)
library(drjacoby)

# 95% interval
quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# calculate geodesic distances. Assumes columns lat/lon
dist_geodesic <- function(A, B, R = 6378137) {
  A <- A[,2:1]
  B <- B[,2:1]
  ret <- matrix(0, nrow = nrow(A), ncol = nrow(B))
  for (i in 1:nrow(A)) {
    for (j in 1:nrow(B)) {
      ret[i, j] <- distHaversine(A[i,], B[j,], r = R)
    }
  }
  return(ret)
}

# convert cartesian coords to lat/lon
cartesian_to_latlon <- function(cartesian_coords) {
  x <- cartesian_coords[, 1]
  y <- cartesian_coords[, 2]
  z <- cartesian_coords[, 3]
  
  r <- sqrt(x^2 + y^2 + z^2)
  lat <- asin(z / r) * 180 / pi
  lon <- atan2(y, x) * 180 / pi
  
  cbind(lat, lon)
}

# sample points uniformly on a sphere
sample_sphere <- function(n, R) {
  phi <- runif(n, 0, 2 * pi)
  costheta <- runif(n, -1, 1)
  theta <- acos(costheta)
  
  x <- R * sin(theta) * cos(phi)
  y <- R * sin(theta) * sin(phi)
  z <- R * cos(theta)
  
  return(cbind(x, y, z))
}

# ------------------------------------------------

# Set parameters
R <- 1
n_features <- 30
n_data <- 100
length_space <- 0.3

# Generate random points on the sphere (for demonstration purposes)
points <- sample_sphere(n_data, R)
points_latlon <- cartesian_to_latlon(points)

# Generate random Fourier features
rff <- sample_sphere(n_features, R)
rff_latlon <- cartesian_to_latlon(rff)

# get geodesic distances between all observed points and features
dist_mat <- dist_geodesic(A = points_latlon,
                          B = rff_latlon,
                          R = R)

# construct feature map
#mag <- rcauchy(n_features)
mag <- rep(1, n_features)
dist_scaled <- sweep(dist_mat, 2, mag, FUN = "*")

feature_maps <- cbind(cos(dist_scaled / length_space), sin(dist_scaled / length_space)) / sqrt(n_features)

# Generate some synthetic data for demonstration purposes
true_beta <- rnorm(n_features * 2)
y <- sweep(feature_maps, 2, true_beta, FUN = "*") |>
  rowSums()

points_latlon |>
  as.data.frame() |>
  mutate(p = plogis(y)) |>
  ggplot() + theme_bw() +
  geom_point(aes(x = lon, y = lat, col = p), size = 3) +
  scale_color_viridis_c(option = "magma", limits = c(0, 1))

# plot_ly(x = ~x, y = ~y, z = ~z, color = ~p, type = 'scatter3d', mode = 'markers',
#         marker = list(size = 5),
#         data = data.frame(x = points[,"x"],
#                           y = points[,"y"],
#                           z = points[,"z"],
#                           p = plogis(y)))

# ------------------------------------------------

dat <- data.frame(N = rep(100, n_data)) |>
  mutate(prev_true = plogis(y),
         k = rbinom(n(), size = N, prob = prev_true))

# Prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  dist_mat = dist_scaled,
  k = dat$k,
  n_trials = dat$N,
  length_space_meanlog = 0,
  length_space_sdlog = 1,
  beta_mean = 0,
  beta_sd = 1,
  tmp = cbind(cos(dist_scaled), sin(dist_scaled)) / sqrt(n_features)
)

# Compile and fit the model
fit <- stan(
  file = "R_ignore/R_scripts/stan/spatial_dist.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

print(fit, pars = c("beta", "length_space"))

if (FALSE) {
  
  df_params <- data.frame(name = sprintf("beta_%s", 1:(2*n_features)),
                          min = -Inf,
                          max = Inf)
  df_params <- define_params(name = "length_space", min = 0, max = 10) |>
    bind_rows(df_params)
  
  r_loglike <- function(params, data, misc) {
    length_space <- params[1]
    feat <- cbind(cos(misc$dist_mat / length_space), sin(misc$dist_mat / length_space)) / sqrt(misc$n_features)
    y_mod <- colSums(t(feat) * params[-1])
    #y_mod <- (misc$tmp %*% params[-1])[,1]
    sum(dbinom(data$k, size = data$N, prob = plogis(y_mod), log = TRUE))
  }
  r_logprior <- function(params, misc) {
    sum(dnorm(params, sd = 10, log = TRUE))
  }
  
  mcmc <- run_mcmc(data = dat,
                   df_params = df_params,
                   misc = list(dist_mat = data_list$dist_mat,
                               n_features = n_features,
                               tmp = data_list$tmp),
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin = 1e3,
                   samples = 1e3,
                   chains = 2,
                   beta_manual = beta_tuned)

#beta_tuned <- seq(0, 1, l = 40)^2
#beta_tuned <- tune_rungs(mcmc, target_acceptance = 0.4)

#plot_mc_acceptance(mcmc)
#plot_trace(mcmc)
#mcmc$diagnostics
}

# ------------------------------------------------

# Extract samples for beta and std
draws <- extract(fit)
beta_samples <- draws$beta
length_space_samples <- draws$length_space
n_draws <- nrow(beta_samples)

# n_draws <- 100
# draws <- sample_chains(mcmc, sample_n = n_draws)
# beta_samples <- draws[,sprintf("beta_%s", 1:(2*n_features))] |>
#   as.matrix()
# length_space_samples <- rep(1, n_draws)

# Plot histogram
#length_space_samples <- rep(1, n_draws)
hist(length_space_samples)

# prediction points
grid_x <- 30
grid_y <- 30
dy <- 180 / grid_y
dx <- 360 / grid_x
points_pred_latlon <- expand_grid(lat = seq(-90 + dy/2, 90 - dy/2, by = dy),
                                  lon = seq(-180 + dx/2, 180 - dx/2, by = dx))

# get geodesic distances between all prediction points and features
dist_mat_pred <- dist_geodesic(A = points_pred_latlon,
                               B = rff_latlon,
                               R = R)

# choose a posterior draw
draw_i <- sample(n_draws, 1)

# construct feature map
length_space_i <- length_space_samples[draw_i]
dist_mat_pred_scaled <- sweep(dist_mat_pred, 2, mag, FUN = "*")
feature_maps_pred <- cbind(cos(dist_mat_pred_scaled / length_space_i), sin(dist_mat_pred_scaled / length_space_i)) / sqrt(n_features)

# predict at each point
y_pred <- sweep(feature_maps_pred, 2, beta_samples[draw_i,], FUN = "*") |>
  rowSums()

tmp <- points_latlon |>
  as.data.frame() |>
  mutate(p = dat$k / dat$N,
         p2 = plogis(y))

points_pred_latlon |>
  mutate(p = plogis(y_pred)) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = lon, y = lat, fill = p)) +
  geom_point(aes(x = lon, y = lat, fill = p), col = grey(0.5), pch = 21, size = 3, data = tmp) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1))


library(tidyverse)
library(rstanarm)
library(rstan)

# 95% interval
quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# ------------------------------------------------

set.seed(2)

n_data <- 100
dat <- data.frame(x = runif(n_data),
                  y = runif(n_data),
                  t = runif(n_data)) |>
  mutate(prev_true = y*(1 - t)*x^3*(1 - x)^3 / 0.5^6 + (1 - y)*t*(x^3 + (1 - x)^3),
         N = ifelse(x < 0.5, 1000, 10),
         k = rbinom(n(), size = N, prob = prev_true))

X <- dat |>
  dplyr::select(x, y) |>
  as.matrix()

n_features <- 20
omega_space <- rbind(rcauchy(n_features),
                     rcauchy(n_features))
omega_time <- rbind(rcauchy(n_features))

# ------------------------------------------------

# Prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  space_mat = X,
  time_vec = cbind(dat$t),
  k = dat$k,
  n_trials = dat$N,
  omega_space = omega_space,
  omega_time = omega_time
)

# Compile and fit the model
fit <- stan(
  file = "R_ignore/R_scripts/stan/spacetime_sigma.stan",
  data = data_list,
  iter = 100,
  chains = 1
)

# Extract samples for beta and std
draws <- extract(fit)
beta_samples <- draws$beta
sigma_space_samples <- draws$sigma_space
sigma_time_samples <- draws$sigma_time

hist(sigma_space_samples)
hist(sigma_time_samples)

# get the posterior mean sigmas
sigma_space_mean <- mean(sigma_space_samples)
sigma_time_mean <- mean(sigma_time_samples)

# Compute the features using the mean of the posterior sigmas
feat <- cbind(cos(sigma_space_mean * X %*% omega_space + sigma_time_mean * dat$t %*% omega_time),
              sin(sigma_space_mean * X %*% omega_space + sigma_time_mean * dat$t %*% omega_time)) / sqrt(n_features)

# Compute the predicted log-odds (z) for all posterior samples
z_samples <- feat %*% t(beta_samples)

# Compute the predicted probabilities for all posterior samples
pred_probs_matrix <- plogis(z_samples)

# Compute the posterior predictive mean probabilities
pred_probs_mean <- rowMeans(pred_probs_matrix)

# Compare the posterior predictive mean probabilities with the observed proportions
plot(pred_probs_mean, dat$k / dat$N, ylim = c(0, 1), xlim = c(0, 1))
abline(a = 0, b = 1, col = "red")

# ------------------------------------------------

n_grid <- 30
X_pred <- expand_grid(x = seq(0, 1, l = n_grid),
                      y = seq(0, 1, l = n_grid)) |>
  as.matrix()

t_pred <- 1
t_vec <- rep(t_pred, nrow(X_pred))

n_draws <- length(sigma_space_samples)

i <- sample(n_draws, 1)

feat_pred <- cbind(cos(sigma_space_samples[i] * X_pred %*% omega_space + sigma_time_samples[i] * t_vec %*% omega_time),
                   sin(sigma_space_samples[i] * X_pred %*% omega_space + sigma_time_samples[i] * t_vec %*% omega_time)) / sqrt(n_features)
#dim(feat_pred)

z <- feat_pred %*% beta_samples[i,]

X_pred |>
  as.data.frame() |>
  mutate(pred = z[,1]) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = plogis(pred))) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  geom_point(aes(x = x, y = y, fill = k / N), pch = 21, col = grey(0.5), size = 4, data = dat)


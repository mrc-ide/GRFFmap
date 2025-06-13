
library(tidyverse)
library(rstanarm)
library(rstan)

# 95% interval
quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# ------------------------------------------------

# Set seed for reproducibility
set.seed(2)

# Generate data
n_data <- 100
dat <- data.frame(x = runif(n_data),
                  y = runif(n_data)) %>%
  mutate(prev_true = plogis(1 * sin(10 * x) + 1 * cos(10 * y)),
         N = ifelse(x < 0.5, 100, 100),
         k = rbinom(n(), size = N, prob = prev_true))

X <- dat %>%
  select(x, y) %>%
  as.matrix()

n_features <- 30
omega <- rbind(rcauchy(n_features),
               rcauchy(n_features))

# Prepare data list for Stan
data_list <- list(
  N = nrow(dat),
  n_features = n_features,
  X = X,
  k = dat$k,
  n_trials = dat$N,
  omega = omega
)

# Compile and fit the model
fit <- stan(
  file = "R_ignore/R_scripts/stan/spatial_sigmaest.stan",
  data = data_list,
  iter = 200,
  chains = 1
)

# Extract samples for beta and std
draws <- extract(fit)
beta_samples <- draws$beta
std_samples <- draws$std

hist(std_samples)

# Compute the mean of std_samples
std_mean <- mean(std_samples)

# Scale omega by the mean of std_samples
omega_mean <- std_mean * omega

# Compute the features using the mean of the posterior samples of std
feat <- cbind(cos(X %*% omega_mean),
              sin(X %*% omega_mean)) / sqrt(n_features)

# Compute the predicted log-odds (z) for all posterior samples
z_samples <- feat %*% t(beta_samples)

# Compute the predicted probabilities for all posterior samples
pred_probs_matrix <- plogis(z_samples)

# Compute the posterior predictive mean probabilities
pred_probs_mean <- rowMeans(pred_probs_matrix)

# Compare the posterior predictive mean probabilities with the observed proportions
plot(pred_probs_mean, dat$k / dat$N, ylim = c(0, 1), xlab = "Predicted Probabilities", ylab = "Observed Proportions")
abline(a = 0, b = 1, col = "red")

# ------------------------------------------------

n_grid <- 30
X_pred <- expand_grid(x = seq(0, 1, l = n_grid),
                      y = seq(0, 1, l = n_grid)) |>
  as.matrix()

n_draws <- length(std_samples)

i <- sample(n_draws, 1)

feat_pred <- cbind(cos(X_pred %*% omega * std_samples[i]),
                   sin(X_pred %*% omega * std_samples[i])) / sqrt(n_features)
#dim(feat_pred)

z <- feat_pred %*% beta_samples[i,]

X_pred |>
  as.data.frame() |>
  mutate(pred = z[,1]) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = plogis(pred))) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  #geom_point(aes(x = x, y = y, fill = prev_true), pch = 21, col = grey(0.5), size = 4, data = dat)
  geom_point(aes(x = x, y = y, fill = k / N), pch = 21, col = grey(0.5), size = 4, data = dat)


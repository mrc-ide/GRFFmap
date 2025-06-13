
library(tidyverse)
library(rstanarm)

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
  mutate(prev_true = plogis(t**sin(10*x) + (1 - t)*cos(10*y)),
         N = ifelse(x < 0.5, 100, 100),
         k = rbinom(n(), size = N, prob = prev_true))

X <- dat |>
  select(x, y) |>
  as.matrix()

n_features <- 30
sigma_x <- 5
omega <- cbind(rnorm(n_features, 0, sigma_x),
               rnorm(n_features, 0, sigma_x))

sigma_t <- 5
#beta <- rcauchy(n_features, location = 0, scale = sigma_t)
beta <- rnorm(n_features, mean = 0, sd = sigma_t)

feat <- cbind(cos(X %*% t(omega) + dat$t %*% t(beta)),
              sin(X %*% t(omega) + dat$t %*% t(beta))) / sqrt(n_features)

# ------------------------------------------------

prior <- normal(location = 0, scale = 10)

#formula_string <- paste0('k / N ~ -1 + ', paste0('x.',1:(2*n_features), collapse = " + "))
formula_string <- paste0('cbind(k, N - k) ~ -1 + ', paste0('x.', 1:(2 * n_features), collapse = " + "))
post1 <- stan_glm(formula = as.formula(formula_string),
                  data = data.frame(k = dat$k, N = dat$N, x = feat),
                  family = binomial(link = "logit"),
                  iter = 1e3,
                  chains = 1,
                  prior = prior,
                  thin = 1)

summary(post1)

z <- sweep(feat, MARGIN = 2, STATS = coef(post1), FUN = "*") |>
  rowSums()

plot(plogis(z), dat$k / dat$N)
abline(a = 0, b = 1)

# ------------------------------------------------

draws <- post1 |>
  as.matrix()

z <- feat %*% t(draws[,1:ncol(feat)])
CrI <- apply(plogis(z), 1, quantile_95)

plot(dat$k / dat$N, ylim = c(0,1))
segments(x0 = 1:n_data, x1 = 1:n_data, y0 = CrI[1,], y1 = CrI[3,], col = "red")

# ------------------------------------------------

n_grid <- 30
X_pred <- expand_grid(x = seq(0, 1, l = n_grid),
                      y = seq(0, 1, l = n_grid)) |>
  as.matrix()

t_pred <- 1
t_vec <- rep(t_pred, nrow(X_pred))

feat_pred <- cbind(cos(X_pred %*%t(omega) + t_vec %*% t(beta)),
                   sin(X_pred %*%t(omega) + t_vec %*% t(beta))) / sqrt(n_features)
dim(feat_pred)

z <- feat_pred %*% t(draws[,1:ncol(feat)])
dim(z)

i <- sample(ncol(z), 1)
X_pred |>
  as.data.frame() |>
  mutate(pred = z[,i]) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = plogis(pred))) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  #geom_point(aes(x = x, y = y, fill = prev_true), pch = 21, col = grey(0.5), size = 4, data = dat)
  geom_point(aes(x = x, y = y, fill = k / N), pch = 21, col = grey(0.5), size = 4, data = dat)

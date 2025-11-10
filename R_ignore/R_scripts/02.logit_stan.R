
library(tidyverse)
library(rstanarm)

# empirical logit two numbers
emplogit <- function(Y, N) {
  log(Y + 0.5) - log(N - Y + 0.5)
}

# 95% interval
quantile_95 <- function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}

# ------------------------------------------------

set.seed(2)

n_data <- 50
dat <- data.frame(x = runif(n_data),
                  y = runif(n_data)) |>
  mutate(prev_true = plogis(1*sin(10*x) + 1*cos(10*y)),
         N = 100,
         k = rbinom(n(), size = N, prob = prev_true),
         ye = emplogit(k, N))

X <- dat |>
  select(x, y) |>
  as.matrix()

n_features <- 100
std <- 5
omega <- cbind(rnorm(n_features, 0, std),
               rnorm(n_features, 0, std))

feat <- cbind(cos(X %*% t(omega)),
              sin(X %*% t(omega))) / sqrt(n_features)

# ------------------------------------------------

prior <- normal(location = 0, scale = 10)

formula_string <- paste0('y ~ -1 + ', paste0('x.',1:(2*n_features), collapse = " + "))
post1 <- stan_glm(as.formula(formula_string),
                  data = data.frame(y = dat$ye, x = feat),
                  family = gaussian(link = "identity"), iter = 1e3, chains = 1, prior = prior, thin = 1)

summary(post1)

z <- sweep(feat, MARGIN = 2, STATS = coef(post1), FUN = "*") |>
  rowSums()

plot(z, dat$ye)
abline(a = 0, b = 1)

# ------------------------------------------------

draws <- post1 |>
  as.matrix()

z <- feat %*% t(draws[,1:ncol(feat)])
CrI <- apply(z, 1, quantile_95)

plot(dat$ye)
segments(x0 = 1:n_data, x1 = 1:n_data, y0 = CrI[1,], y1 = CrI[3,], col = "red")

# ------------------------------------------------

n_grid <- 30
X_pred <- expand_grid(x = seq(0, 1, l = n_grid),
                      y = seq(0, 1, l = n_grid)) |>
  as.matrix()

feat_pred <- cbind(cos(X_pred %*%t(omega)),
                   sin(X_pred %*%t(omega))) / sqrt(n_features)
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
  geom_point(aes(x = x, y = y, fill = prev_true), pch = 21, col = grey(0.5), size = 4, data = dat)

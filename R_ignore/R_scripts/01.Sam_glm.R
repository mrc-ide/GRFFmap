library(tidyverse)
library(glmnet)

# empirical logit two numbers
emplogit <- function(Y, N) {
  log(Y + 0.5) - log(N - Y + 0.5)
}

#set.seed(1)

n_data <- 50
dat <- data.frame(x = runif(n_data),
                  y = runif(n_data)) |>
  mutate(prev_true = y,
         N = 100,
         k = rbinom(n(), size = N, prob = prev_true),
         ye = emplogit(k, N))

X <- dat |>
  select(x, y) |>
  as.matrix()

n_features <- 10
std <- 1
omega <- cbind(rnorm(n_features, 0, std),
               rnorm(n_features, 0, std))

feat <- cbind(cos(X %*% t(omega)),
              sin(X %*% t(omega))) / sqrt(n_features)
dim(feat)

opt_lambda <- 1e-1
fit <- glmnet(feat, dat$ye, alpha = 0, lambda = opt_lambda, intercept = TRUE)

n_grid <- 30
X_pred <- expand_grid(x = seq(0, 1, l = n_grid),
                      y = seq(0, 1, l = n_grid)) |>
  as.matrix()

feat_pred <- cbind(cos(X_pred %*%t(omega)),
                   sin(X_pred %*%t(omega))) / sqrt(n_features)
dim(feat_pred)

f3 <- predict(fit, s = opt_lambda, newx = feat_pred)

X_pred |>
  as.data.frame() |>
  mutate(f3 = f3) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = x, y = y, fill = plogis(f3))) +
  scale_fill_viridis_c(option = "magma") +
  geom_point(aes(x = x, y = y, fill = prev_true), pch = 21, col = grey(0.5), size = 4, data = dat)


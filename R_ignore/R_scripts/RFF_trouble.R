# RFF_trouble.R
#
# Author: Bob Verity
# Date: 2025-07-27
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# I'm having some trouble wih random fourier features. I'm trying to estimate
# parameters of the GP (mainly the variance of the process, or "sill") but am
# finding that the RFF approximation is so bad that estimates are horribly
# biased. I can't seem to do anything to improve the approximation, which makes
# me think perhaps this is a bug or mistake.
#
# This script gives a minimal example of the problem
#
# ------------------------------------------------------------------

library(tidyverse)
library(mvtnorm)
library(MASS)
library(randtoolbox)

# ----------------------------------------------------------------
# Simulate some data

#set.seed(1)

# function to simulate from a 1D GP with RBF kernel
simulate_gp_rbf <- function(x, length_scale = 1, sill = 1, mu = 0) {
  n <- length(x)
  dists_sq <- as.matrix(dist(x))^2
  K <- sill * exp(-dists_sq / (2 * length_scale^2))
  y <- as.vector(MASS::mvrnorm(1, mu = rep(mu, n), Sigma = K))
  data.frame(x = x, y = y)
}

# set some parameters
length_scale <- 0.5
sill <- 1 # (variance of GP)
nugget <- 1e-3

# simulate data
df_dat <- simulate_gp_rbf(x = scale(seq(0, 10, l = 11)),
                          length_scale = length_scale,
                          sill = sill)
n <- nrow(df_dat)

# plot the data
df_dat |>
  ggplot() + theme_bw() +
  geom_point(aes(x = x, y = y)) +
  ylim(c(-3, 3)*sqrt(sill))

# ----------------------------------------------------------------
# Estimate the sill using the true likelihood

# calculate the kernel matrix between points
dists_sq <- as.matrix(dist(df_dat$x))^2
K <- exp(-dists_sq / (2 * length_scale^2))

# loop through range of values of sill. For each one calculate loglikelihood
df_sill <- data.frame(sill = seq(0.1, 5, 0.01), loglike = NA)
for (i in 1:nrow(df_sill)) {
  Sigma <- df_sill$sill[i]*K + diag(nugget, n)
  df_sill$loglike[i] <- mvtnorm::dmvnorm(x = df_dat$y, sigma = Sigma, log = TRUE)
}

# plot likelihood function
df_sill |>
  mutate(like = exp(loglike - max(loglike))) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = sill, y = like)) +
  geom_vline(xintercept = sill, linetype = "dashed")

# so far so good...the likelihood peaks around the true value of the sill

# ----------------------------------------------------------------
# Random fourier features approach

# draw features and create approximate version of the kernel matrix
n_features <- 1e2
#omega <- rnorm(n_features, sd = 1/length_scale)
omega <- qnorm(halton(n_features), sd = 1/length_scale)
proj <- df_dat$x %*% t(omega)
phi <- cbind(cos(proj), sin(proj)) / sqrt(n_features)
K_approx <- phi %*% t(phi)

# perform same likelihood calculation as before, but now using the approximated kernel matrix
df_sill$loglike2 <- NA
for (i in 1:nrow(df_sill)) {
  Sigma <- df_sill$sill[i]*K_approx + diag(nugget, n)
  df_sill$loglike2[i] <- mvtnorm::dmvnorm(x = df_dat$y, sigma = Sigma, log = TRUE)
}

# compare the likelihoods from the true vs. RFF approximation
df_sill |>
  mutate(exact = exp(loglike - max(loglike)),
         RFF = exp(loglike2 - max(loglike2))) |>
  dplyr::select(sill, exact, RFF) |>
  pivot_longer(cols = -1) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = sill, y = value, col = name))

# The likelihood based on the RFF approximation of the kernel matrix is miles off

# ----------------------------------------------------------------
# Explore issues related to the matrix

# for any row in the kernel matrix, compare the true and approximated values
i <- 10
plot(K[i,], type = 'l')
lines(K_approx[i,], col = "red")

# this approximation seems strangely bad to me, and doesn't seem to improve with more features. Bug?

# ----------------------------------------------------------------
# Is the fit OK, even if sill estimation is poor?

# estimate sill via maximum likelihood
opt <- optim(1, fn = function(v) {
  Sigma <- v*K_approx + diag(nugget, n)
  -mvtnorm::dmvnorm(x = df_dat$y, sigma = Sigma, log = TRUE)
}, method = "Brent", lower = 0.01, upper = 1e2)
sill_ML <- opt$par


sill_ML

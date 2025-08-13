data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2*n_features] feat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real mu_mean;
  real mu_sd;
  real<lower=0> sigma_shape;
  real<lower=0> sigma_rate;
}

parameters {
  vector[2 * n_features] beta;
  real mu;
  real<lower=0> sigma;
}

model {
  // Priors
  beta ~ normal(0, sigma);
  mu ~ normal(mu_mean, mu_sd);
  sigma ~ gamma(sigma_shape, sigma_rate);
  
  // Likelihood
  k ~ binomial_logit(n_trials, mu + feat * beta);
}


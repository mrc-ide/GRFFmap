data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2*n_features] feat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real<lower=0> sill_shape;
  real<lower=0> sill_rate;
  real<lower=0> nugget_shape;
  real<lower=0> nugget_rate;
  real<lower=0> mu_mean;
  real<lower=0> mu_sd;
}

parameters {
  vector[N] z;
  vector[2 * n_features] beta;
  real<lower=0> sill;
  real<lower=0> nugget;
  real mu;
}

model {
  // Priors
  sill ~ gamma(sill_shape, sill_rate);
  nugget ~ gamma(nugget_shape, nugget_rate);
  beta ~ normal(0, sqrt(sill));
  mu ~ normal(mu_mean, mu_sd);
  
  // Likelihood
  z ~ normal(mu + feat * beta, sqrt(nugget));
  k ~ binomial_logit(n_trials, z);
}


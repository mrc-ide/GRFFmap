data {
  int<lower=1> n_sites;
  int<lower=1> n_features;
  matrix[n_sites, 2*n_features] feat;
  int<lower=0> n_pos[n_sites];
  int<lower=0> n_samp[n_sites];
  real mu_mean;
  real<lower=0> mu_sd;
  real<lower=0> sigma_shape;
  real<lower=0> sigma_rate;
  real<lower=0> pi_nugget_shape1;
  real<lower=0> pi_nugget_shape2;
}

parameters {
  vector[2 * n_features] beta;
  real mu;
  real<lower=0> sigma;
  real<lower=0, upper=1> pi_nugget;
  vector[n_sites] epsilon;
}

model {
  // Priors
  beta ~ normal(0, sigma);
  mu ~ normal(mu_mean, mu_sd);
  sigma ~ gamma(sigma_shape, sigma_rate);
  pi_nugget ~ beta(pi_nugget_shape1, pi_nugget_shape2);
  epsilon ~ normal(0, sqrt(pi_nugget)*sigma);
  
  // Likelihood
  n_pos ~ binomial_logit(n_samp, mu + feat * beta + epsilon);
}


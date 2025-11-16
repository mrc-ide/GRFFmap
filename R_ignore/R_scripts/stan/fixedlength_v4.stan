data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2*n_features] feat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real<lower=0> sill_shape;
  real<lower=0> sill_rate;
  real<lower=0> mu_mean;
  real<lower=0> mu_sd;
}

parameters {
  vector[2 * n_features] beta;
  real<lower=0> sill;
  real mu;
}

model {
  // Priors
  beta ~ normal(0, sqrt(sill));
  sill ~ gamma(sill_shape, sill_rate);
  mu ~ normal(mu_mean, mu_sd);
  
  // Likelihood
  k ~ binomial_logit(n_trials, mu + feat * beta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(k[n] | n_trials[n], mu + dot_product(feat[n], beta));
  }
}

data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2*n_features] feat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real<lower=0> beta_mean;
  real<lower=0> beta_sd;
}

parameters {
  vector[2 * n_features] beta;
}

model {
  // Priors
  beta ~ normal(beta_mean, beta_sd);

  // Likelihood
  k ~ binomial_logit(n_trials, feat * beta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(k[n] | n_trials[n], dot_product(feat[n], beta));
  }
}

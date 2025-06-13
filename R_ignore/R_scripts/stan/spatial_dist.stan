data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, n_features] dist_mat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real<lower=0> length_space_meanlog;
  real<lower=0> length_space_sdlog;
  real<lower=0> beta_mean;
  real<lower=0> beta_sd;
  matrix[N, 2*n_features] tmp;
}

parameters {
  real<lower=0> length_space;
  vector[2 * n_features] beta;     // Regression coefficients
}

transformed parameters {
  matrix[N, 2 * n_features] feat;  // Fourier features matrix
  
  matrix[N, n_features] cos_feat = cos(dist_mat / length_space);
  matrix[N, n_features] sin_feat = sin(dist_mat / length_space);
  
  feat = append_col(cos_feat, sin_feat) / sqrt(n_features);
}

model {
  // Priors
  //length_space ~ lognormal(length_space_meanlog, length_space_sdlog);
  length_space ~ uniform(0, 10);
  beta ~ normal(beta_mean, beta_sd);

  // Likelihood
  k ~ binomial_logit(n_trials, feat * beta);
  //k ~ binomial_logit(n_trials, tmp * beta);
}


data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2] X;                  // Input features (x, y)
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  matrix[2, n_features] omega; // Fixed raw omega
}

parameters {
  real<lower=0> std;               // Standard deviation for omega
  vector[2 * n_features] beta;     // Regression coefficients
}

transformed parameters {
  matrix[N, 2 * n_features] feat;
  
  // Compute Fourier features using matrix multiplication
  matrix[N, n_features] cos_feat = cos(X * std * omega);
  matrix[N, n_features] sin_feat = sin(X * std * omega);
  feat = append_col(cos_feat, sin_feat) / sqrt(n_features);
}

model {
  // Priors
  std ~ normal(0, 10);
  beta ~ normal(0, 10);

  // Likelihood
  k ~ binomial_logit(n_trials, feat * beta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(k[n] | n_trials[n], dot_product(feat[n], beta));
  }
}

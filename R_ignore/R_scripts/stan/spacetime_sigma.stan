data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2] space_mat;
  matrix[N, 1] time_vec;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  matrix[2, n_features] omega_space;
  matrix[1, n_features] omega_time;
}

parameters {
  real<lower=0> sigma_space;
  real<lower=0> sigma_time;
  vector[2 * n_features] beta;
}

transformed parameters {
  matrix[N, 2 * n_features] feat;
  
  // Compute Fourier features using matrix multiplication
  matrix[N, n_features] cos_feat = cos(space_mat * sigma_space * omega_space + time_vec * sigma_time * omega_time);
  matrix[N, n_features] sin_feat = sin(space_mat * sigma_space * omega_space + time_vec * sigma_time * omega_time);
  feat = append_col(cos_feat, sin_feat) / sqrt(n_features);
}

model {
  // Priors
  sigma_space ~ normal(0, 10);
  sigma_time ~ normal(0, 10);
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

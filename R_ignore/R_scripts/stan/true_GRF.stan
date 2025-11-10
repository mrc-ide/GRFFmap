data {
  int<lower=0> N;           // Number of observations
  int<lower=0> x[N];        // Numerator (successes)
  int<lower=0> n[N];        // Denominator (trials)
  matrix[N, N] dist_space;  // Precomputed spatial distances
  matrix[N, N] dist_time;   // Precomputed temporal distances
}

parameters {
  real alpha;               // Intercept
  real<lower=0> sigma;      // Gaussian Process variance
  real<lower=0> rho_space;  // Spatial correlation length
  real<lower=0> rho_time;   // Temporal correlation length
  vector[N] eta;            // Latent field
}

transformed parameters {
  vector[N] eta_std;
  matrix[N, N] K;
  
  // Covariance matrix for the Gaussian Process
  K = square(sigma) * exp(-dist_space / rho_space - dist_time / rho_time);
  //for (i in 1:N) {
  //  K[i, i] = K[i, i] + 1e-6;
  //}
  
  // Standardize eta
  eta_std = cholesky_decompose(K) \ eta;
}

model {
  // Priors
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 1);
  rho_space ~ normal(0, 5);
  rho_time ~ normal(0, 5);
  eta ~ normal(0, 1);
  
  // Likelihood
  for (i in 1:N) {
    x[i] ~ binomial_logit(n[i], alpha + eta_std[i]);
  }
}

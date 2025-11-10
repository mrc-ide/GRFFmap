data {
  int<lower=1> N;                  // Number of data points
  int<lower=1> n_features;         // Number of Fourier features
  matrix[N, 2*n_features] feat;
  int<lower=0> k[N];               // Observed successes
  int<lower=0> n_trials[N];        // Number of trials
  real<lower=0, upper=1> pi_nugget;
}

parameters {
  vector[2 * n_features] beta;
  real mu;
  real<lower=0> sigma;
  vector[N] epsilon;
}

model {
  // Priors
  beta ~ normal(0, sigma);
  mu ~ normal(0, 5);
  sigma ~ gamma(1, 0.1);
  epsilon ~ normal(0, sqrt(pi_nugget)*sigma);
  
  // Likelihood
  k ~ binomial_logit(n_trials, mu + feat * beta + epsilon);
}


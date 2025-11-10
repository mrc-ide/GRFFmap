data {
  int<lower=1> N;
  int<lower=1> n_features;
  matrix[N, 2*n_features] feat;
  vector[N] z;
}

parameters {
  vector[2*n_features] beta;
  real<lower=0> sigma;
}

model {
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  z ~ normal(feat * beta, sigma);
}

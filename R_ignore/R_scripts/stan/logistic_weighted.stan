data {
  int<lower=1> n_data;
  int<lower=1> n_features;
  matrix[n_data, 2*n_features] feat;
  vector[n_data] n;
  vector[n_data] z;
}

parameters {
  vector[2*n_features] beta;
  real<lower=0> sigma;
  real mu;
}

model {
  beta ~ normal(0, 1);
  sigma ~ normal(0, 100);
  mu ~ normal(0, 100);
  
  z ~ normal(feat * beta + mu, sigma / n);
}

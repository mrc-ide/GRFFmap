data {
  int<lower=1> N;               // Number of data points
  int<lower=0> k[N];            // Observed counts
  int<lower=1> n[N];            // Binomial totals
  matrix[N, N] K;               // Precomputed RBF kernel matrix
}

parameters {
  vector[N] z;                  // Latent field values
  real mu;
  real<lower=0> sigma;
}

model {
  // Prior over latent field
  z ~ multi_normal(rep_vector(0, N), K);
  mu ~ normal(0, 5);
  sigma ~ gamma(1, 0.1);

  // Binomial likelihood with logit link
  k ~ binomial_logit(n, mu + sigma*z);
}

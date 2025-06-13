data {
  int<lower=1> N;          // number of data points
  matrix[N, N] cor_mat;
  vector[N] z;             // response variable
}
parameters {
  real<lower=0> psill;
  real mu;
}
model {
  z ~ multi_normal(rep_vector(mu, N), psill * cor_mat);
}

data {
  int<lower=0> d;
  real<lower=0> nu;
  matrix[d, d] S;
}

transformed data {
  vector[d] mu = rep_vector(0.0, d);
}

parameters {
  vector[d] x1;
  real<lower=0> x2;
}

transformed parameters {
  vector[d] x = x1 * sqrt(x2);;
}

model {
  x1 ~ multi_normal(mu, S);
  x2 ~ inv_gamma(nu * 0.5, nu * 0.5);
}

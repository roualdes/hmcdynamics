data {
  int<lower=0> d;
  real<lower=0> nu;
  matrix[d, d] S;
}

transformed data {
  vector[d] mu = rep_vector(0.0, d);
  real<lower=0.> nnu = 0.5 * nu;
}

parameters {
  vector[d] mvnormal;
  real<lower=0> invgamma;
}

transformed parameters {
  vector[d] x = mvnormal * sqrt(invgamma);;
}

model {
  mvnormal ~ multi_normal(mu, S);
  invgamma ~ inv_gamma(nnu, nnu);
}

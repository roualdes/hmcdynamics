data {
  real<lower=0> nu;
  real<lower=0> s;
}
parameters {
  real x1;
  real<lower=0> x2;
}

transformed parameters {
  real x = x1 * sqrt(x2);
}

model {
  x1 ~ normal(0, 1);
  x2 ~ inv_gamma(nu * 0.5, nu * s * s * 0.5);
}

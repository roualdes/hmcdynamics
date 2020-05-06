data {
  int<lower=0> d;
  matrix[d, d] S;
}

transformed data {
  vector[d] mu = rep_vector(0.0, d);
}

parameters {
  vector[d] x;
}

model {
  x ~ multi_normal(mu, S);
}

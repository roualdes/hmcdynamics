data {
  int<lower=0> d;
  vector[d] mu;
  matrix[d, d] S;
}
parameters {
  vector[d] x;
}
model {
  x ~ multi_normal(mu, S);
}

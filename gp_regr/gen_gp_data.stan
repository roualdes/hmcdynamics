transformed data {
  real rho = 5.5;
  real alpha = 3;
  real sigma = 2;
}

parameters {}
model {}

generated quantities {
  int<lower=1> N = 11;
  real x[11] = {-10, -8, -6, -4, -2, 0.0, 2, 4, 6, 8, 10};
  real f[11];
  vector[11] y;
  int k[11];

    matrix[11, 11] cov =   cov_exp_quad(x, alpha, rho)
                       + diag_matrix(rep_vector(1e-10, 11));
    matrix[11, 11] L_cov = cholesky_decompose(cov);
    f = multi_normal_cholesky_rng(rep_vector(0, 11), L_cov);

    for (n in 1:11) {
      y[n] = normal_rng(f[n], sigma);
      k[n] = poisson_rng(exp(f[n]));
    }

}

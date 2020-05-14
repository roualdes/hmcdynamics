// Adapted from
// https://github.com/google-research/google-research/blob/master/neutra/logistic_reg_pystan.py
//
// Copyright 2020 The Google Research Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

data {
  int<lower=0> n;               // number of observations
  int<lower=0> d;               // number of predictors
  int<lower=0,upper=1> y[n];    // outputs
  matrix[n,d] x;                // inputs
}
parameters {
  vector[d] z;
  vector<lower=0>[d] local_scale;
  real<lower=0> global_scale;
}
model {
  z ~ normal(0, 1);
  local_scale ~ gamma(0.5, 0.5);
  global_scale ~ gamma(0.5, 0.5);

  y ~ bernoulli_logit(x * (z .* local_scale * global_scale));
}

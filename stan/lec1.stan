data {
  int<lower=1> N;                       // Total observations (300)
  int<lower=1> Npop;                    // Number of populations (15)
  int<lower=1, upper=Npop> Index[N];    // Population index for each observation (length N)
  vector[N] Lengths;
  vector[N] Weights;
}

parameters {
  vector[Npop] log_a;       // log(a_i) for each population
  vector[Npop] log_b;       // log(b_i) for each population
  real mu_log_a;            // Hypermean for log(a_i)
  real mu_log_b;            // Hypermean for log(b_i)
  real<lower=0> sigma_a;    // Hyper-sd for log(a_i)
  real<lower=0> sigma_b;    // Hyper-sd for log(b_i)
  real<lower=0> sigma;      // Residual error
}

model {
  // Hyperpriors (from your prior image)
  mu_log_a ~ normal(0, 1000);
  mu_log_b ~ normal(0, 1000);
  sigma_a ~ cauchy(0, 2.5);
  sigma_b ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 2.5);

  // Population-level priors
  log_a ~ normal(mu_log_a, sigma_a);
  log_b ~ normal(mu_log_b, sigma_b);

  // Likelihood: log(W_ij) = log(a_i) + b_i * log(L_ij) + error
  for (j in 1:N) {
    target += normal_lpdf(
      log(Weights[j]) | log_a[Index[j]] + log_b[Index[j]] * log(Lengths[j]), sigma);
  }
}

generated quantities {
  vector[N] log_weight_pred;
  for (j in 1:N) {
    log_weight_pred[j] = normal_rng(log_a[Index[j]] + log_b[Index[j]] * log(Lengths[j]), sigma);
  }
}
data {
  int<lower=0> N;              // Number of obs
  int<lower=0> S;              // Number of sites
  int<lower=0> K;              // Number of statuses
  int<lower=0> G;              // Number of prey groups
  vector[N] TL;                // Fish length
  int<lower=1, upper=S> site[N];    // Site index (1:S)
  int<lower=1, upper=K> status[N];  // Status index (1:K)
  int<lower=0, upper=1> prey_mat[N, G]; // N x G prey presence/absence
  real<lower=0> meanTLs_site[S];     // mean TL per site (for prediction)
  real<lower=0> meanTLs_status[K];   // mean TL per status (for prediction)
}
parameters {
  vector[G] alpha;             // Intercept for each prey group
  matrix[G, S] beta_site;      // Site effect for each prey group
  matrix[G, K] beta_status;    // Status effect for each prey group
  vector[G] beta_TL;           // TL effect for each prey group
}
model {
  alpha ~ normal(0, 1);
  to_vector(beta_site) ~ normal(0, 1);
  to_vector(beta_status) ~ normal(0, 1);
  beta_TL ~ normal(0, 1);

  for (i in 1:N) {
    for (g in 1:G) {
      prey_mat[i, g] ~ bernoulli_logit(
        alpha[g]
        + beta_site[g, site[i]]
        + beta_status[g, status[i]]
        + beta_TL[g] * TL[i]
      );
    }
  }
}
generated quantities {
  // Site-wise posterior predictions for each prey group
  matrix[G, S] p_site;
  for (g in 1:G) {
    for (s in 1:S) {
      p_site[g, s] = inv_logit(alpha[g]
                               + beta_site[g, s]
                               + beta_TL[g] * meanTLs_site[s]);
    }
  }
  // Status-wise posterior predictions for each prey group
  matrix[G, K] p_status;
  for (g in 1:G) {
    for (k in 1:K) {
      p_status[g, k] = inv_logit(alpha[g]
                                 + beta_status[g, k]
                                 + beta_TL[g] * meanTLs_status[k]);
    }
  }
}

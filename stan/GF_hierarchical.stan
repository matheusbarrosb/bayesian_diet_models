data {
  int<lower=1> N;            // number of observations
  int<lower=1> S;            // number of species
  int<lower=1> I;            // number of sites
  int<lower=1> TT;            // number of restoration statuses
  int<lower=1,upper=S> species[N];
  int<lower=1,upper=I> site[N];
  int<lower=1,upper=TT> status[N];
  vector[N] TL;
  vector[N] FW;
  vector[N] GW;
}

parameters {
  // Intercepts for each species
  vector[S] alpha;

  // Regression coefficients
  vector[S] beta_TL;
  vector[S] beta_FW;

  // Site main effects
  matrix[S, I] beta_site;

  // Restoration status main effects
  matrix[S, TT] beta_RS;

  // Site × Restoration status interaction
  array[S] matrix[I, TT] beta_site_RS;

  // Error SD per species
  vector<lower=0>[S] sigma;

  // Hyperpriors for the group-level SDs (each coefficient group has its own scale)
  real<lower=0> tau_alpha;
  real<lower=0> tau_TL;
  real<lower=0> tau_FW;
  vector<lower=0>[S] tau_site;
  vector<lower=0>[S] tau_RS;
  array[S] vector<lower=0>[I] tau_site_RS;
}

model {
  // Hyperpriors (Cauchy for SDs, as is common for scale)
  tau_alpha ~ cauchy(0, 2.5);
  tau_TL    ~ cauchy(0, 2.5);
  tau_FW    ~ cauchy(0, 2.5);
  tau_site  ~ cauchy(0, 2.5);
  tau_RS    ~ cauchy(0, 2.5);
  for (s in 1:S) tau_site_RS[s] ~ cauchy(0, 2.5);

  // Priors
  alpha   ~ normal(0, tau_alpha);
  beta_TL ~ normal(0, tau_TL);
  beta_FW ~ normal(0, tau_FW);

  for (s in 1:S) {
    beta_site[s] ~ normal(0, tau_site[s]);
    beta_RS[s] ~ normal(0, tau_RS[s]);
    for (i in 1:I)
      beta_site_RS[s, i] ~ normal(0, tau_site_RS[s][i]);
    sigma[s] ~ cauchy(0, 2.5);
  }

  // Likelihood
  for (n in 1:N) {
    GW[n] ~ lognormal(
      alpha[species[n]]
      + beta_site[species[n], site[n]]
      + beta_RS[species[n], status[n]]
      + beta_site_RS[species[n], site[n], status[n]]
      + beta_TL[species[n]] * TL[n]
      + beta_FW[species[n]] * FW[n],
      sigma[species[n]]
    );
  }
}

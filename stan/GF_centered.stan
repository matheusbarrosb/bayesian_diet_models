data {
  int<lower=1> N;            // number of observations
  int<lower=1> S;            // number of species
  int<lower=1> I;            // number of sites
  int<lower=1> TT;           // number of restoration statuses
  int<lower=1,upper=S> species[N];
  int<lower=1,upper=I> site[N];
  int<lower=1,upper=TT> status[N];
  vector[N] TL;
  vector[N] FW;
  vector[N] GW;
}

parameters {
  vector[S] alpha;                     // Baseline for each species

  // Deviations for site (except baseline)
  matrix[S, I-1] beta_site;            // S species, I-1 site deviations

  // Deviations for restoration status (except baseline)
  matrix[S, TT-1] beta_RS;             // S species, TT-1 status deviations

  // Deviations for site × status interaction (excluding baseline site and status)
  array[S] matrix[I-1, TT-1] beta_site_RS; // S species, I-1 x TT-1 deviations

  // Continuous predictors
  vector[S] beta_TL;
  vector[S] beta_FW;

  // Error SD per species
  vector<lower=0>[S] sigma;

  // Hyperpriors for group-level SDs
  real<lower=0> tau_alpha;
  real<lower=0> tau_TL;
  real<lower=0> tau_FW;
  vector<lower=0>[S] tau_site;
  vector<lower=0>[S] tau_RS;
  array[S] vector<lower=0>[I-1] tau_site_RS;
}

model {
  // Hyperpriors
  tau_alpha ~ cauchy(0, 2.5);
  tau_TL    ~ cauchy(0, 2.5);
  tau_FW    ~ cauchy(0, 2.5);
  tau_site  ~ cauchy(0, 2.5);
  tau_RS    ~ cauchy(0, 2.5);
  for (s in 1:S) tau_site_RS[s] ~ cauchy(0, 2.5);

  // Priors for main effects
  alpha   ~ normal(0, tau_alpha);
  beta_TL ~ normal(0, tau_TL);
  beta_FW ~ normal(0, tau_FW);

  for (s in 1:S) {
    beta_site[s] ~ normal(0, tau_site[s]);
    beta_RS[s]   ~ normal(0, tau_RS[s]);
    for (i in 1:(I-1))
      beta_site_RS[s, i] ~ normal(0, tau_site_RS[s][i]);
    sigma[s] ~ cauchy(0, 2.5);
  }

  // Likelihood
  for (n in 1:N) {
    // Site effect: 0 if baseline, else use deviation
    real site_eff = site[n] == 1 ? 0 : beta_site[species[n], site[n] - 1];
    // Status effect: 0 if baseline, else use deviation
    real status_eff = status[n] == 1 ? 0 : beta_RS[species[n], status[n] - 1];
    // Interaction: only for non-baseline site and status, else 0
    real inter_eff = (site[n] == 1 || status[n] == 1) ? 0
      : beta_site_RS[species[n], site[n] - 1, status[n] - 1];

    GW[n] ~ lognormal(
      alpha[species[n]]
      + site_eff
      + status_eff
      + inter_eff
      + beta_TL[species[n]] * TL[n]
      + beta_FW[species[n]] * log(FW[n]),
      sigma[species[n]]
    );
  }
}
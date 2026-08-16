data {
  int<lower=1> N;         // number of data points
  int<lower=1> NS;        // number of sites
  int<lower=1> S;         // number of species
  int<lower=1> TT;        // number of restoration statuses
  int<lower=1, upper=S> spp[N];
  int<lower=1, upper=NS> site[N];
  int<lower=1, upper=TT> RS[N];
  real<lower=0> GW[N];    // gut weight
  real<lower=0> FW[N];    // wet weight
  real<lower=0> TL[N];    // total length
}
parameters {
  real alpha[S];             // intercept per species
  real beta_site[S, NS];     // site effect per species
  real beta_RS[S, TT];       // restoration status effect per species
  real beta_TL[S];           // TL effect per species
  real beta_FW[S];           // FW effect per species
  real<lower=0> sigma[S];    // residual SD per species
}
model {
  to_array_1d(alpha) ~ normal(0, 5);
  to_array_1d(beta_site) ~ normal(0, 2);
  to_array_1d(beta_RS) ~ normal(0, 2);
  to_array_1d(beta_TL) ~ normal(0, 2);
  to_array_1d(beta_FW) ~ normal(0, 2);
  to_array_1d(sigma) ~ cauchy(0, 0.1);

  for (n in 1:N) {
    int s = spp[n];
    int i = site[n];
    int t = RS[n];
    GW[n] ~ lognormal(
      alpha[s]
      + beta_site[s, i]
      + beta_RS[s, t]
      + beta_TL[s] * TL[n]
      + beta_FW[s] * FW[n],
      sigma[s]
    );
  }
}
 
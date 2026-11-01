data {
  
  int<lower=1> N;                       
  int<lower=1> S;                       
  int<lower=1> I;                       
  int<lower=1> TT;                      
  
  int<lower=1,upper=S> species[N];      
  int<lower=1,upper=I> site[N];         
  int<lower=1,upper=TT> status[N];      
  
  vector<lower=0,upper=1>[N] GF;        
  vector[N] TL;         
  vector[S] global_mean_TL;
  
}
parameters {
  
  vector[S]           alpha_c;
  matrix[S, I-1]      z_site;
  matrix[S, TT-1]     z_RS;
  array[S] matrix[I-1, TT-1] z_int;
  vector[S]           beta_TL;          
  
  vector<lower=0, upper=pi()/2>[S] u_phi;
  vector<lower=0, upper=pi()/2>[S] u_site;
  vector<lower=0, upper=pi()/2>[S] u_RS;
  vector<lower=0, upper=pi()/2>[S] u_int;
  
}
transformed parameters {
  
  vector[S] phi        = 10 * tan(u_phi);
  vector[S] sigma_site = 1.0 * tan(u_site);
  vector[S] sigma_RS   = 1.0 * tan(u_RS);
  vector[S] sigma_int  = 1.0 * tan(u_int);
  
  vector[S] alpha = alpha_c - beta_TL .* global_mean_TL;
  matrix[S, I-1]  beta_site = diag_pre_multiply(sigma_site, z_site);   // DEFLECTIONS FOR NON-BASELINE SITES
  matrix[S, TT-1] beta_RS   = diag_pre_multiply(sigma_RS, z_RS);       // DEFLECTIONS FOR NON-BASELINE STATUS
  array[S] matrix[I-1, TT-1] beta_int;                                 // NON-BASELINE INTERACTIONS
  
  vector[N] A;                          
  vector[N] B;                          
  
  for (s in 1:S) {
    beta_int[s] = sigma_int[s] * z_int[s];
  }
  
  for (n in 1:N) {
    int sp = species[n];
    
    real site_eff   = (site[n] == 1)   ? 0.0 : beta_site[sp, site[n] - 1];
    real status_eff = (status[n] == 1) ? 0.0 : beta_RS[sp, status[n] - 1];
    real int_eff    = (site[n] == 1 || status[n] == 1) ? 0.0 : beta_int[sp, site[n] - 1, status[n] - 1];
    
    real mu = inv_logit(alpha_c[sp] + site_eff + status_eff + int_eff + beta_TL[sp] * (TL[n] - global_mean_TL[sp]));
    
    A[n] = mu * phi[sp];
    B[n] = (1.0 - mu) * phi[sp];
  }
  
}
model {
  
  for (n in 1:N) {
    target += beta_lpdf(GF[n] | A[n] , B[n]); 
  }

  alpha      ~ normal(0, 2);                
  beta_TL    ~ normal(0, 2);
  
  to_vector(z_site) ~ std_normal();
  to_vector(z_RS)   ~ std_normal();
  for (s in 1:S) {
    to_vector(z_int[s]) ~ std_normal();
  }
  
}
generated quantities {
  
  array[S] matrix[I, TT] y_hat_int;      
  vector[S]              mu_GF;          
  
  { 
    for (s in 1:S) {
      for (i in 1:I) {
        for (t in 1:TT) {
          real site_eff   = (i == 1) ? 0.0 : beta_site[s, i - 1];
          real status_eff = (t == 1) ? 0.0 : beta_RS[s, t - 1];
          real int_eff    = (i == 1 || t == 1) ? 0.0 : beta_int[s, i - 1, t - 1];
          
          real mu_hat = inv_logit(alpha[s] + site_eff + status_eff + int_eff + beta_TL[s] * global_mean_TL[s]);
          y_hat_int[s, i, t] = beta_rng(mu_hat * phi[s], (1.0 - mu_hat) * phi[s]);
        }
      }
      mu_GF[s] = mean(to_vector(y_hat_int[s])); 
    }
  }
  
}
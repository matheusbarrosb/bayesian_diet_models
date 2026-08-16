data {
  
  int<lower=1> N;                       
  int<lower=1> S;                       
  int<lower=1> I;                       
  int<lower=1> TT;                      
  
  int<lower=1,upper=S> species[N];      
  int<lower=1,upper=I> site[N];         
  int<lower=1,upper=TT> status[N];      
  
  vector<lower=0,upper=1>[N] GF;        
  int<lower=0,upper=1> Z[N];              
  vector[N] TL;         
  vector[S] global_mean_TL;
  
}
parameters {
  
  vector[S]           alpha;            
  matrix[S, I-1]      beta_site;        // DEFLECTIONS FOR NON-BASELINE SITES
  matrix[S, TT-1]     beta_RS;          // DEFLECTIONS FOR NON-BASELINE STATUS
  array[S] matrix[I-1, TT-1] beta_int;  // NON-BASELINE INTERACTIONS
  vector[S]           beta_TL;          
  vector<lower=0>[S]  phi;              
  
  vector[S]           omega;            
  vector[S]           omega_TL;         
  
  vector<lower=0>[S]  sigma_site;       
  vector<lower=0>[S]  sigma_RS;         
  vector<lower=0>[S]  sigma_int;        
  
}
transformed parameters {
  
  vector[N] A;                          
  vector[N] B;                          
  vector[N] LP_omega;                   
  
  for (n in 1:N) {
    int sp = species[n];
    
    real site_eff   = (site[n] == 1)   ? 0.0 : beta_site[sp, site[n] - 1];
    real status_eff = (status[n] == 1) ? 0.0 : beta_RS[sp, status[n] - 1];
    real int_eff    = (site[n] == 1 || status[n] == 1) ? 0.0 : beta_int[sp, site[n] - 1, status[n] - 1];
    
    real mu = inv_logit(alpha[sp] + site_eff + status_eff + int_eff + beta_TL[sp] * TL[n]);
    
    A[n] = mu * phi[sp];
    B[n] = (1.0 - mu) * phi[sp];
    
    LP_omega[n] = omega[sp] + omega_TL[sp] * TL[n];
  }
  
}
model {
  
  for (n in 1:N) {
    if (Z[n] == 0) {
      target += beta_lpdf(GF[n] | A[n] , B[n]); 
    }
  }
  
  Z ~ bernoulli_logit(LP_omega); 

  alpha      ~ normal(0, 2);                
  beta_TL    ~ normal(0, 2);
  omega      ~ normal(0, 2);
  omega_TL   ~ normal(0, 2);
  phi        ~ cauchy(0, 10);
  
  sigma_site ~ cauchy(0, 1);
  sigma_RS   ~ cauchy(0, 1);
  sigma_int  ~ cauchy(0, 1);
  
  for (s in 1:S) {
    beta_site[s] ~ normal(0, sigma_site[s]);
    beta_RS[s]   ~ normal(0, sigma_RS[s]);
    for (i in 1:(I-1)) {
      beta_int[s, i] ~ normal(0, sigma_int[s]);
    }
  }
  
}
generated quantities {
  
  array[S] matrix[I, TT] y_hat_int;      
  vector[S]              mu_GF;          
  
  { 
    vector[S] mean_theta = rep_vector(0.0, S);
    
    for (s in 1:S) {
      mean_theta[s] = inv_logit(omega[s] + omega_TL[s] * global_mean_TL[s]);
    }
    
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
      mu_GF[s] = mean(to_vector(y_hat_int[s])) * (1.0 - mean_theta[s]); 
    }
  }
  
}
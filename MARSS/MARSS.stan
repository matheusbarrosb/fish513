data {
  
  int<lower=0>              N;              // number of years
  int<lower=0>              M;              // number of time series
  array [M]int<lower=0>     states;
  int<lower=0>              S;              // number of states
  array [M]int<lower=0>     obsVariances;   // observation variance map
  int<lower=0>              n_obsvar;       // number of observation variances
  array [S+1]int<lower=0>   proVariances;   // process variance map
  int<lower=0>              n_provar;       // number of process variances
  array [S+1]int<lower=0>   trends;         // trend map
  int<lower=0>              n_trends;       // number of trends
  int<lower=0>              n_pos;          // number of non-NA values
  array [n_pos]int<lower=0> col_indx_pos;
  array [n_pos]int<lower=0> row_indx_pos;
  int<lower=0>              est_trend;
  int<lower=0>              n_A;
  int<lower=0>              est_nu;
  array [n_A+2]int<lower=0> est_A;
  vector[n_pos]             y;              // data
  int                       family;         // 1 = normal, 2 = binomial, 3 = poisson, 4 = gamma, 5 = lognormal

}
parameters {
  
  vector<lower=0>[S]                     x0;              // initial states
  vector<lower=2>[est_nu]                nu;              // nu, constrainted to be > 2
  matrix <lower=-3,upper=3>[N-1,S]       pro_dev;         // process deviations
  vector[n_trends]                       U;
  vector[n_A]                            A; 
  array [n_provar]real<lower=0,upper=1>  sigma_process;   // process variation (SD)
  array [n_obsvar]real<lower=0>          sigma_obs;       // observation variation
  real<lower=0,upper=1>                  phi;

}
transformed parameters {
  
  real<lower=0>          sigma_process_real = sigma_obs[1]*sigma_process[1];
  matrix[N, M]           pred;
  matrix[N,S]            x; 
  vector[S]              Uvec;
  vector[M]              Avec;
  real<lower=-1,upper=1> phi_real = 2*phi - 1;
  
  for(i in 1:M) Avec[i] = 0;
  for(i in 1:n_A) Avec[est_A[i]] = A[i];
  
  for(i in 1:S) {
    
    if(est_trend) {
      
      Uvec[i] = U[trends[i]]; // map shared trends
      
    } else {
      
     Uvec[i] = 0;
     
    }
  }

  
  for(s in 1:S) {
    
    x[1,s] = x0[s]; // assign initial states
    
  }  


for(t in 2:N) {
    for(s in 1:S) {
      // process equation
      x[t,s] = x[t-1,s] + pro_dev[t-1,s] * sigma_process[proVariances[s]] * sigma_obs[1];
      
      if(est_trend == 1) {
      
        x[t,s] = x[t-1,s] + Uvec[s];
     
      
      }
    }
  }

  // map predicted states to time series
  for(m in 1:M) {
    for(t in 1:N) {
      
      pred[t,m] = x[t,states[m]] + Avec[m];
      
    }
  }
}
model {

// PRIORS ---------------------------------------------------------------------
  phi ~ beta(1.5,1.5) T[1e-4,1];
  
  for(i in 1:n_obsvar) {
    
    sigma_obs[i] ~ gamma(1,0.01); 
    
  }
  
  for(s in 1:n_provar) {
    
    sigma_process[s] ~ gamma(1,0.01);
  }
  
  for(i in 1:n_trends) {
      U[i] ~ normal(0, 0.1); 
  }
  
  
  if(est_nu ==1) {
    
    nu[1] ~ gamma(2, 0.1);
    
  for(t in 1:(N-1)) {
      for(s in 1:S) {
        
        pro_dev[t,s] ~ student_t(nu[1], 0, 0.1);
        
      }
    }
  } else {

    for(s in 1:S) {
      
      pro_dev[1,s] ~ normal(0, 2);
      
      }
      
    for(t in 2:(N-1)) {
      for(s in 1:S) {
        
       // pro_dev[t,s] ~ normal((2*phi-1)*pro_dev[t-1,s],10);
       pro_dev[t,s] ~ normal(0, 2);
        
      }
    }
  }


  // LIKELIHOODS ---------------------------------------------------------------
  if(family == 1) {
    
    for(i in 1:n_pos) y[i] ~ normal(pred[col_indx_pos[i], row_indx_pos[i]], sigma_obs[obsVariances[row_indx_pos[i]]]);
    
  }
  // if(family == 2) {
  //   for(i in 1:n_pos) y_int[i] ~ bernoulli_logit(pred[col_indx_pos[i], row_indx_pos[i]]);
  // }
  //  if(family == 3) {
  //    for(i in 1:n_pos) y_int[i] ~ poisson_log(pred[col_indx_pos[i], row_indx_pos[i]]);
  // }
  if(family == 4) {
    for(i in 1:n_pos) y[i] ~ gamma(sigma_obs[obsVariances[row_indx_pos[i]]], sigma_obs[obsVariances[row_indx_pos[i]]] ./ pred[col_indx_pos[i], row_indx_pos[i]]);
  }
  if(family == 5) {
    for(i in 1:n_pos) y[i] ~ lognormal(pred[col_indx_pos[i], row_indx_pos[i]], sigma_obs[obsVariances[row_indx_pos[i]]]);
  }
}
generated quantities {
  
  vector[n_pos] log_lik;
  matrix[N, S]   lambda;

  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);
  // if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(pred[col_indx_pos[n], row_indx_pos[n]]));
  // if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(pred[col_indx_pos[n], row_indx_pos[n]]));
  if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | sigma_obs[obsVariances[row_indx_pos[n]]], sigma_obs[obsVariances[row_indx_pos[n]]] ./ exp(pred[col_indx_pos[n], row_indx_pos[n]]));
  if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);
  
  for (s in 1:S){ 
    for (t in 2:N) lambda[t,s] = (pred[t,s] - pred[t-1,s])/(t - (t-1));
    }
  
}
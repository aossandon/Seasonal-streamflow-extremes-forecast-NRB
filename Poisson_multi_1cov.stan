functions {

}

data {
  int<lower=0> N;
  int<lower=0> S;
  int<lower=0> y1[N];
  int<lower=0> y2[N];
  int<lower=0> y3[N];
  int<lower=0> y4[N];
  int<lower=0> y5[N];
  vector [N] covar1;
}

transformed data {
  
}
parameters {
  vector[S] alpha_loc0;
  vector[S] alpha_loc1;
}
model {
  vector[N] tmp_lamb;
  matrix [N, S] lp;
  vector[N] lp_cop;
  vector[S] zer;
  
  zer=rep_vector(0,S);
  alpha_loc0 ~ normal(0,1000);
  alpha_loc1 ~ normal(0,1000);
  // alpha_loc2 ~ normal(0,100);
  for (r in 1:S){
    tmp_lamb=exp(alpha_loc0[r] +alpha_loc1[r]*covar1);
    for (i in 1:N){
      if (r == 1){
      lp[i,1]=poisson_lpmf(y1[i] | tmp_lamb[i]);
      }else if (r == 2){
      lp[i,2]=poisson_lpmf(y2[i] | tmp_lamb[i]);
      }else if (r == 3){
      lp[i,3]=poisson_lpmf(y3[i] | tmp_lamb[i]);
      }else if (r == 4){
      lp[i,4]=poisson_lpmf(y4[i] | tmp_lamb[i]);
      }else{
      lp[i,5]=poisson_lpmf(y5[i] | tmp_lamb[i]);
      }
    
    }
  }
  
  
target += sum(col(lp,1)+col(lp,2)+col(lp,3)+col(lp,4)+col(lp,5)); 
}



generated quantities{
  vector[N] tmp_lamb;
  matrix<lower=0>[N, S] y_rep;
  matrix [N, S]  lp;
  vector[N] log_lik;
  
  for (r in 1:S){
    tmp_lamb=exp(alpha_loc0[r] +alpha_loc1[r]*covar1);
    for (i in 1:N){
      if (r == 1){
      lp[i,1]=poisson_lpmf(y1[i] | tmp_lamb[i]);
      y_rep[i,1]=poisson_rng(tmp_lamb[i]);
      }else if (r == 2){
      lp[i,2]=poisson_lpmf(y2[i] | tmp_lamb[i]);
      y_rep[i,2]=poisson_rng(tmp_lamb[i]);
      }else if (r == 3){
      lp[i,3]=poisson_lpmf(y3[i] | tmp_lamb[i]);
      y_rep[i,3]=poisson_rng(tmp_lamb[i]);
      }else if (r == 4){
      lp[i,4]=poisson_lpmf(y4[i] | tmp_lamb[i]);
      y_rep[i,4]=poisson_rng(tmp_lamb[i]);
      }else{
      lp[i,5]=poisson_lpmf(y5[i] | tmp_lamb[i]);
      y_rep[i,5]=poisson_rng(tmp_lamb[i]);
      }
    
    }
  }
  
log_lik=col(lp,1)+col(lp,2)+col(lp,3)+col(lp,4)+col(lp,5); 
}

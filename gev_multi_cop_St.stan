functions {


   vector gev_v_log_v_2cov (vector y, vector covar1, vector covar2, real alpha_loc0, real alpha_loc1, real alpha_loc2, real log_scale0, real log_scale1 , real log_scale2, real shape){
    real inv_shape;
    real inv_shape_p1;
    real neg_inv_shape;
    vector[rows(y)] loc;
    vector[rows(y)] scale;
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;

    N = rows(y);
    loc= alpha_loc0 + alpha_loc1*covar1 + alpha_loc2*covar2;
    scale=exp(log_scale0 +log_scale1*covar1 +log_scale2*covar2);
    if(shape == 0){
      z = (y - loc)./scale;
      for(n in 1:N){
        lp[n] = log(scale[n]) + z[n] + exp(-z[n]);
      }
    }else{
      inv_shape = 1.0/shape;
      inv_shape_p1 = inv_shape + 1;
      neg_inv_shape = -inv_shape;

      z = 1 + (y - loc) * shape ./ scale;
      for(n in 1:N){
        lp[n] = log(scale[n]) + inv_shape_p1*log(z[n]) + pow(z[n],neg_inv_shape);
      }
    }
    // print("loc ",loc);
    // print("scale ",scale);
    // print("shape ",shape);
    return -lp;  
  }
   
   vector gev_v_log_v_1cov (vector y, vector covar1, real alpha_loc0, real alpha_loc1, real log_scale0, real log_scale1, real shape){
    real inv_shape;
    real inv_shape_p1;
    real neg_inv_shape;
    vector[rows(y)] loc;
    vector[rows(y)] scale;
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;

    N = rows(y);
    loc= alpha_loc0 + alpha_loc1*covar1;
    scale=exp(log_scale0 +log_scale1*covar1);
    if(shape == 0){
      z = (y - loc)./scale;
      for(n in 1:N){
        lp[n] = log(scale[n]) + z[n] + exp(-z[n]);
      }
    }else{
      inv_shape = 1.0/shape;
      inv_shape_p1 = inv_shape + 1;
      neg_inv_shape = -inv_shape;

      z = 1 + (y - loc) * shape ./ scale;
      for(n in 1:N){
        lp[n] = log(scale[n]) + inv_shape_p1*log(z[n]) + pow(z[n],neg_inv_shape);
      }
    }
    // print("loc ",loc);
    // print("scale ",scale);
    // print("shape ",shape);
    return -lp;  
  }
  
  
  vector gev_v_log_v_St (vector y, real alpha_loc0, real log_scale0, real shape){
    real inv_shape;
    real inv_shape_p1;
    real neg_inv_shape;
    real loc;
    real scale;
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;

    N = rows(y);
    loc= alpha_loc0;
    scale=exp(log_scale0);
    if(shape == 0){
      z = (y - loc)./scale;
      for(n in 1:N){
        lp[n] = log(scale) + z[n] + exp(-z[n]);
      }
    }else{
      inv_shape = 1.0/shape;
      inv_shape_p1 = inv_shape + 1;
      neg_inv_shape = -inv_shape;

      z = 1 + (y - loc) * shape ./ scale;
      for(n in 1:N){
        lp[n] = log(scale) + inv_shape_p1*log(z[n]) + pow(z[n],neg_inv_shape);
      }
    }
    // print("loc ",loc);
    // print("scale ",scale);
    // print("shape ",shape);
    return -lp;  
  }
   
   vector gev_v_log_v_2cov_SgST (vector y, vector covar1, vector covar2, real alpha_loc0, real alpha_loc1, real alpha_loc2, real log_scale, real shape){
    real inv_shape;
    real scale;
    real inv_shape_p1;
    real neg_inv_shape;
    vector[rows(y)] loc;
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;

    N = rows(y);
    loc= alpha_loc0 + alpha_loc1*covar1 + alpha_loc2*covar2;
    scale=exp(log_scale);
    if(shape == 0){
      z = (y - loc)/scale;
      lp = log(scale) + z + exp(-z);
    }else{
      inv_shape = 1.0/shape;
      inv_shape_p1 = inv_shape + 1;
      neg_inv_shape = -inv_shape;

      z = 1 + (y - loc) * shape / scale;
      for(n in 1:N){
        lp[n] = log(scale) + inv_shape_p1*log(z[n]) + pow(z[n],neg_inv_shape);
      }
    }
    // print("loc ",loc);
    // print("scale ",scale);
    // print("shape ",shape);
    return -lp;  
  }
  
  
  vector gev_v_log_v_1cov_SgST (vector y, vector covar1, real alpha_loc0, real alpha_loc1, real log_scale, real shape){
    real inv_shape;
    real scale;
    real inv_shape_p1;
    real neg_inv_shape;
    vector[rows(y)] loc;
    vector[rows(y)] z;
    vector[rows(y)] lp;
    int N;

    N = rows(y);
    loc= alpha_loc0 + alpha_loc1*covar1;
    scale=exp(log_scale);
    if(shape == 0){
      z = (y - loc)/scale;
      lp = log(scale) + z + exp(-z);
    }else{
      inv_shape = 1.0/shape;
      inv_shape_p1 = inv_shape + 1;
      neg_inv_shape = -inv_shape;

      z = 1 + (y - loc) * shape / scale;
      for(n in 1:N){
        lp[n] = log(scale) + inv_shape_p1*log(z[n]) + pow(z[n],neg_inv_shape);
      }
    }
    // print("loc ",loc);
    // print("scale ",scale);
    // print("shape ",shape);
    return -lp;  
  }
  
  vector gev_v_cdf_v_2cov (vector q, vector covar1, vector covar2, real alpha_loc0, real alpha_loc1, real alpha_loc2, real log_scale0, real log_scale1 , real log_scale2, real shape) {
    
    vector[rows(q)] p;
    vector[rows(q)] r; 
    vector[rows(q)] loc;
    vector[rows(q)] scale;
    real neg_inv_shape;
    vector[rows(q)] s;
    real es;
    int N;
    
    N = rows(q);
    loc= alpha_loc0 + alpha_loc1*covar1 + alpha_loc2*covar2;
    scale=exp(log_scale0 +log_scale1*covar1 +log_scale2*covar2);
    neg_inv_shape = -1.0/shape;
    s = (q - loc) ./ scale;

    if(shape == 0.0){
        p = exp(-exp(-(q-loc)./scale));
    }else{
      for (n in 1:N){
        // s= 1 +  scale / shape * (q[n] - loc);
        // s[n] = (q[n] - loc[n]) / scale;
        es = s[n]*shape;
        if(es<= -1){
          if(shape>0){
            p[n] = 0;
          }else{
            p[n] = 1;
          }
        }else{
          p[n] = exp(-pow(1 + es,neg_inv_shape));
        }
        // p[n] = exp(-pow(fmax(s,0),neg_inv_shape));
      }
      
      
        
    }
    return p;
  }
  
  vector gev_v_cdf_v_1cov (vector q, vector covar1, real alpha_loc0, real alpha_loc1, real log_scale0, real log_scale1, real shape) {
    
    vector[rows(q)] p;
    vector[rows(q)] r; 
    vector[rows(q)] loc;
    vector[rows(q)] scale;
    real neg_inv_shape;
    vector[rows(q)] s;
    real es;
    int N;
    
    N = rows(q);
    loc= alpha_loc0 + alpha_loc1*covar1;
    scale=exp(log_scale0 +log_scale1*covar1);
    neg_inv_shape = -1.0/shape;
    s = (q - loc) ./ scale;

    if(shape == 0.0){
        p = exp(-exp(-(q-loc)./scale));
    }else{
      for (n in 1:N){
        // s= 1 +  scale / shape * (q[n] - loc);
        // s[n] = (q[n] - loc[n]) / scale;
        es = s[n]*shape;
        if(es<= -1){
          if(shape>0){
            p[n] = 0;
          }else{
            p[n] = 1;
          }
        }else{
          p[n] = exp(-pow(1 + es,neg_inv_shape));
        }
        // p[n] = exp(-pow(fmax(s,0),neg_inv_shape));
      }
      
      
        
    }
    return p;
  }
  
  vector gev_v_cdf_v_St (vector q, real alpha_loc0, real log_scale0, real shape) {
    
    vector[rows(q)] p;
    vector[rows(q)] r; 
    real loc;
    real scale;
    real neg_inv_shape;
    vector[rows(q)] s;
    real es;
    int N;
    
    N = rows(q);
    loc= alpha_loc0;
    scale=exp(log_scale0);
    neg_inv_shape = -1.0/shape;
    s = (q - loc) ./ scale;

    if(shape == 0.0){
        p = exp(-exp(-(q-loc)./scale));
    }else{
      for (n in 1:N){
        // s= 1 +  scale / shape * (q[n] - loc);
        // s[n] = (q[n] - loc[n]) / scale;
        es = s[n]*shape;
        if(es<= -1){
          if(shape>0){
            p[n] = 0;
          }else{
            p[n] = 1;
          }
        }else{
          p[n] = exp(-pow(1 + es,neg_inv_shape));
        }
        // p[n] = exp(-pow(fmax(s,0),neg_inv_shape));
      }
      
      
        
    }
    return p;
  }
  
  vector gev_v_cdf_v_2cov_SgST (vector q, vector covar1, vector covar2, real alpha_loc0, real alpha_loc1, real alpha_loc2, real log_scale, real shape) {
    
    vector[rows(q)] p;
    vector[rows(q)] r; 
    vector[rows(q)] loc;
    real scale;
    real neg_inv_shape;
    vector[rows(q)] s;
    real es;
    int N;
    
    N = rows(q);
    loc= alpha_loc0 + alpha_loc1*covar1 + alpha_loc2*covar2;
    scale=exp(log_scale);
    neg_inv_shape = -1.0/shape;
    s = (q - loc) / scale;

    if(shape == 0.0){
        p = exp(-exp(-(q-loc)/scale));
    }else{
      for (n in 1:N){
        // s= 1 +  scale / shape * (q[n] - loc);
        // s[n] = (q[n] - loc[n]) / scale;
        es = s[n]*shape;
        if(es<= -1){
          if(shape>0){
            p[n] = 0;
          }else{
            p[n] = 1;
          }
        }else{
          p[n] = exp(-pow(1 + es,neg_inv_shape));
        }
        // p[n] = exp(-pow(fmax(s,0),neg_inv_shape));
      }
      
      
        
    }
    return p;
  }
  
  
  vector gev_v_cdf_v_1cov_SgST (vector q, vector covar1, real alpha_loc0, real alpha_loc1, real log_scale, real shape) {
    
    vector[rows(q)] p;
    vector[rows(q)] r; 
    vector[rows(q)] loc;
    real scale;
    real neg_inv_shape;
    vector[rows(q)] s;
    real es;
    int N;
    
    N = rows(q);
    loc= alpha_loc0 + alpha_loc1*covar1;
    scale=exp(log_scale);
    neg_inv_shape = -1.0/shape;
    s = (q - loc) / scale;

    if(shape == 0.0){
        p = exp(-exp(-(q-loc)/scale));
    }else{
      for (n in 1:N){
        // s= 1 +  scale / shape * (q[n] - loc);
        // s[n] = (q[n] - loc[n]) / scale;
        es = s[n]*shape;
        if(es<= -1){
          if(shape>0){
            p[n] = 0;
          }else{
            p[n] = 1;
          }
        }else{
          p[n] = exp(-pow(1 + es,neg_inv_shape));
        }
        // p[n] = exp(-pow(fmax(s,0),neg_inv_shape));
      }
      
      
        
    }
    return p;
  }
  
  vector eliptical_copula_log_v (matrix logf, matrix F, matrix L, vector zeros){
    vector[dims(F)[1]] lp;
    vector[dims(F)[2]] u;
    vector[dims(F)[2]] log_phi_u;
    int N;
    int T;

    // logf and F are N x T matricies
    N = dims(F)[2];
    T = dims(F)[1];

    for(t in 1:T){
      for(n in 1:N){
        u[n] = inv_Phi(F[t,n]);
        log_phi_u[n] = normal_lpdf(u[n] | 0, 1);
      }

      // print("sum(logf[t]) ",logf[t]);
      // print("sum(log_phi_u)",log_phi_u);
      // print("multi normal ", multi_normal_cholesky_lpdf(u | zeros, L));

      lp[t] = sum(logf[t]) - sum(log_phi_u) + multi_normal_cholesky_lpdf(u | zeros, L);
    }
    return lp;  
  }

}

data {
  int<lower=0> N;
  int<lower=0> S;
  int<lower=0> Nc;
  matrix [N, S] y;
}

transformed data {
  
}
parameters {
  vector<lower=-0.5,upper=0.5>[S] shape;
  vector[S] Log_scale0;
  // vector[S] Log_scale1;
  // vector[S] Log_scale2;
  
  vector[S] alpha_loc0;
  // vector[S] alpha_loc1;
  // vector[S] alpha_loc2;
  corr_matrix[S] cop_mat;
  // cov_matrix[S] covmat_scale;
  // cov_matrix[S] covmat_loc0;
  // cov_matrix[S] covmat_shape;
}
model {
  vector[N] tmp2;
  vector[N] tmp3;
  matrix [N, S] lp;
  matrix [N, S] Fs;
  vector[S] zer;
  matrix[S,S] L_copula;
  vector[N] eclp;
  
  zer=rep_vector(0,S);
  // priors
  cop_mat ~ lkj_corr(2);
  Log_scale0 ~ normal(0,100);
  // Log_scale1 ~ normal(0,100);
  // Log_scale2 ~ normal(0,100);
  shape ~ normal(0,1);
  alpha_loc0 ~ normal(0,100);
  // alpha_loc1 ~ normal(0,100);
  // alpha_loc2 ~ normal(0,100);
  for (r in 1:S){
    // print("station ",r);
    tmp2 = gev_v_log_v_St(col(y,r), alpha_loc0[r], Log_scale0[r], shape[r]);
    tmp3 = gev_v_cdf_v_St(col(y,r), alpha_loc0[r], Log_scale0[r], shape[r]);
    for (i in 1:N){
      lp[i,r] = tmp2[i];
      Fs[i,r] = tmp3[i]; 
    }
    
  }
  
  L_copula = cholesky_decompose(cop_mat);
  eclp = eliptical_copula_log_v(lp, Fs, L_copula, zer);
  

  target += sum(eclp); 
}



generated quantities{
  vector[N] tmp2;
  vector[N] tmp3;
  matrix [N, S] lp;
  matrix [N, S] Fs;
  vector[N] log_lik;
  vector[S] zer;
  matrix[S,S] L_copula;

  zer=rep_vector(0,S);
  for (r in 1:S){
    tmp2 = gev_v_log_v_St(col(y,r), alpha_loc0[r], Log_scale0[r], shape[r]);
    tmp3 = gev_v_cdf_v_St(col(y,r), alpha_loc0[r], Log_scale0[r], shape[r]);
    for (i in 1:N){
      lp[i,r] = tmp2[i];
      Fs[i,r] = tmp3[i]; 
    }
    
  }
  
  L_copula = cholesky_decompose(cop_mat);
  log_lik = eliptical_copula_log_v(lp, Fs, L_copula, zer);

}

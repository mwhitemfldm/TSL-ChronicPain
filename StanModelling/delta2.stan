data {
  int<lower=1> N; // number of subjects
  int<lower=1> Tn; // number of trials
  matrix[Tn, N] seq; // sequence outcomes
  matrix[Tn, N] pred; // participant estimation
  vector<lower=0>[N] m_init; //  initial estimation mu, first prediction
}

parameters {
  //Individual level
  vector<lower=0,upper=1>[N] alphaCoeff; //learning rate
  vector<lower=0>[N] precCoeff; //participant prediction error
  
  // Group-level parameters
  vector<lower=0>[2]  muGroup; // group mean
  vector<lower=0>[2]  varGroup; // group variance
}

transformed parameters {
  vector<lower=0,upper=1>[N] alphaSubj;
  vector<lower=0>[N] precSubj;
  
  matrix[Tn, N] estMu;
  matrix[Tn, N] predErr;
  
  for (s in 1:N) {
    
    alphaSubj[s] = Phi_approx(muGroup[1] + varGroup[1] * alphaCoeff[s]);
    precSubj[s] = muGroup[2] + varGroup[2] * precCoeff[s];
    
    estMu[1, s] = m_init[s];
    predErr[1,s] = 0;
    
    for (t in 2:Tn){
      predErr[t,s] = seq[(t-1),s] - estMu[(t-1),s];
      estMu[t,s] = estMu[(t-1),s] + alphaSubj[s] * predErr[t,s];
    }
  }
  
}

model {
  
  for (s in 1:N){
    
    for (t in 1:Tn){
      
      target += normal_lpdf(pred[t,s] | estMu[t,s], precSubj[s]);
    }
    
   /// subject level priors
   target += normal_lpdf(alphaCoeff[s]|0, 1.0);
   target += normal_lpdf(precCoeff[s]|0, 1.0);
  }
  
  //////// hyperpriors
  target += normal_lpdf(muGroup|0,1);
  //target += normal_lpdf(varGroup|1,0.2);
  //target += cauchy_lpdf(varGroup|0,5);
  target += student_t_lpdf(varGroup|3,0,1);
}


generated quantities {
  matrix[Tn, N] y_sim;
  matrix[Tn, N] log_lik;
  for (s in 1:N){
    for (trial in 1:Tn) {
      y_sim[trial,s] = normal_rng(estMu[trial,s], precSubj[s]);
      log_lik[trial,s] = normal_lpdf(pred[trial,s] | estMu[trial,s], precSubj[s]);
    }
  }
  
}

data {
  int<lower=1> N; // number of subjects
  int<lower=1> Tn; // number of trials
  matrix[Tn, N] seq; // sequence
  matrix[Tn, N] pred; // predictions
  vector<lower=0>[N] m_init; // initial estimate
  int<lower=1> xVar1; // initial prior variance
}

parameters {
  //Subject level parameters
  vector<lower=0>[N] driftCoeff; // random walk drift variance
  vector<lower=0>[N] noiseCoeff; // outcome noise variance
  vector<lower=0>[N] precCoeff; // precision error variance

  // Group-level parameters
  vector<lower=0>[3]  muGroup; // group drift variance
  vector<lower=0>[3]  varGroup; // group noise variance
}

transformed parameters {
  
  vector<lower=0>[N] driftSubj;
  vector<lower=0>[N] noiseSubj;
  vector<lower=0>[N] precSubj;
  
  matrix[Tn, N] k_gain;
  matrix[Tn, N] x_mu; // mean estimate
  matrix[Tn, N] x_var; // estimate variance
  
 
  for (s in 1:N) {
    
    // Non centred paramaterisation
    //driftSubj[s] = Phi_approx(muGroup[1] + varGroup[1] * driftCoeff[s])*30;
    driftSubj[s] = muGroup[1] + varGroup[1] * driftCoeff[s];
    noiseSubj[s] = muGroup[2] + varGroup[2] * noiseCoeff[s];
    precSubj[s] = muGroup[3] + varGroup[3] * precCoeff[s];
    
    // Initialise
    x_mu[1, s] = m_init[s]; // initial estimate is first prediction
    x_var[1,s] = xVar1; // initial variance of prior
    k_gain[1,s] = 1; // K gain initially 1
    
    for (t in 2:(Tn)){
      x_mu[t,s] = (noiseSubj[s]*x_mu[(t-1),s] + x_var[(t-1),s]*seq[(t-1),s])/(noiseSubj[s] + x_var[(t-1),s]); //update mean estimate
      k_gain[t,s] = (x_var[(t-1),s])/(x_var[(t-1),s]+noiseSubj[s]);
      x_var[t,s] = noiseSubj[s]*x_var[(t-1),s]/(noiseSubj[s] + x_var[(t-1),s]) + driftSubj[s]; //update variance estimate
    }
  }
}

model {  
  for (s in 1:N){
  
    for (t in 1:Tn){
      target += normal_lpdf(pred[t,s] | x_mu[t,s], precSubj[s]); // prediction drawn from gaussian, mean - estmate, var predError
    }
    
    /// subject level priors
   target += normal_lpdf(driftCoeff[s]|0, 1.0);
   target += normal_lpdf(noiseCoeff[s]|0, 1.0);
   target += normal_lpdf(precCoeff[s]|0, 1.0);
    
  }
  
  
  //////// hyperpriors
  target += normal_lpdf(muGroup|0,1);
 // target += normal_lpdf(varGroup|1,0.2);
  //target += cauchy_lpdf(varGroup|0, 5);
  target += student_t_lpdf(varGroup|3,0,1);
}

generated quantities {
  matrix[Tn, N] log_lik;
  matrix[Tn, N] y_sim;
  for (s in 1:N){
    for (t in 1:Tn) {
      y_sim[t,s] = normal_rng(x_mu[t,s], precSubj[s]);
      log_lik[t,s] = normal_lpdf(pred[t,s] | x_mu[t,s], precSubj[s]);
    }
  }
}













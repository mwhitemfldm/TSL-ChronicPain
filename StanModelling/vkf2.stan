////////////////////////////////////
  //
  // Volatile Kalman Filter adapted from:
  // https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007963#pcbi.1007963.e001
//
  ////////////////////////////////////
  
data {
  int<lower=1> N; //No. participants
  int<lower=1> Tn; //No. trials
  matrix[Tn, N] seq; // sequence outcomes
  matrix[Tn, N] pred; // participant estimation
  vector<lower=0>[N] m_init; //initial m
}

parameters {
  // subject level parameters
  
  vector<lower=0>[N] lambdaCoeff; //volatility update parameter
  vector<lower=0>[N] v1Coeff; // initial process noise variance
  vector<lower=0>[N] w1Coeff; // initial variance estimate
  vector<lower=0>[N] sigmaCoeff; // outcome noise variance
  //vector<lower=0.01>[N] precCoeff; // precision error variance
  vector<lower=0>[N] precCoeff; // precision error variance

  // Group-level parameters
  vector<lower=0>[5]  muGroup; // group drift variance
  vector<lower=0>[5]  varGroup; // group noise variance
}

transformed parameters {
  
  // subject level transformed parameters
  vector<lower=0>[N] v1Subj; // initial process noise variance
  vector<lower=0>[N] w1Subj; // initial variance estimate
  vector<lower=0>[N] sigmaSubj; // outcome noise variance
  vector<lower=0, upper=1>[N] lambdaSubj; //volatility update parameter
  vector<lower=0>[N] precSubj; // precision error variance
  
  
  matrix[Tn, N] k;//kalman gain
  matrix[Tn, N] v; // volatility
  matrix[Tn, N] x_mu; // mean estimate
  matrix[Tn, N] w; //  estimated variance
  
  for (s in 1:N) {
    
    // Non centred paramaterisation
    w1Subj[s] = muGroup[1] + varGroup[1] * w1Coeff[s];
    v1Subj[s] = muGroup[2] + varGroup[2] * v1Coeff[s];
    sigmaSubj[s] = muGroup[3] + varGroup[3] * sigmaCoeff[s];
    lambdaSubj[s] = Phi_approx(muGroup[4] + varGroup[4] * lambdaCoeff[s]);
    precSubj[s] = muGroup[5] + varGroup[5] * precCoeff[s];
    
    x_mu[1, s] = m_init[s]; // initial estimate is first prediction
    w[1,s] = w1Subj[s]; // initial variance of prior
    v[1,s] = v1Subj[s]; // initial volatility
    k[1,s] = 1; // K gain initially 1

    
    for (t in 2:(Tn)){
      
      k[t,s] = (w[(t-1),s] + v[(t-1),s])/(w[(t-1),s] + v[(t-1),s] + sigmaSubj[s]);// kalman gain
      x_mu[t,s] = x_mu[(t-1),s] + k[t,s]*(seq[(t-1),s] - x_mu[(t-1),s]); //update mean estimate
      w[t,s] = (1 - k[t,s])*(w[(t-1),s] + v[(t-1),s]); //update variance estimate
      v[t,s] = v[(t-1),s] + lambdaSubj[s]*(square(x_mu[t,s] - x_mu[(t-1),s]) + w[(t-1),s] + w[(t),s] - 2*((1-k[t,s])* w[(t-1),s]) - v[(t-1),s]);
    }
  }
}


model {
  //////// prior: hyperparameters
  
  for (s in 1:N){
  
    for (t in 1:Tn){
      target += normal_lpdf(pred[t,s] | x_mu[t,s], precSubj[s]); // prediction drawn from gaussian, mean - estmate, var predError
    }
    
   /// subject level priors
   target += normal_lpdf(w1Coeff[s]|0, 2.5);
   target += normal_lpdf(v1Coeff[s]|0, 2.5);
   target += normal_lpdf(sigmaCoeff[s]|0, 2.5);
   target += normal_lpdf(lambdaCoeff[s]|0, 2.5);
   target += normal_lpdf(precCoeff[s]|0, 2.5);
  }
  
  //////// hyperpriors
  target += normal_lpdf(muGroup|0,2.5);
  //target += normal_lpdf(varGroup|1,0.2);
  //target += cauchy_lpdf(varGroup|0,10);
  //target += normal_lpdf(errorVar|1,0.2);
  target += student_t_lpdf(varGroup|3,0,1);
}
  
generated quantities {
  matrix[Tn, N] y_sim;
  matrix[Tn, N] log_lik;
  for (s in 1:N){
    for (t in 1:Tn) {
      y_sim[t,s] = normal_rng( x_mu[t,s], precSubj[s]);
      log_lik[t,s] = normal_lpdf(pred[t,s] | x_mu[t,s], precSubj[s]);
    }
  }
}
  



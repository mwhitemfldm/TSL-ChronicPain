//PearceHall-RL model

data {
  int<lower=1> N; // number of subjects
  int<lower=1> Tn; // number of trials
  matrix[Tn, N] seq; // sequence outcomes
  matrix[Tn, N] pred; // participant estimation
  vector<lower=0>[N] m_init; //  initial estimation mu, first prediction
}

parameters {
  // Subject level parameters
  vector<lower=0>[N] kappaCoeff; // fixed component of learning rate
  vector<lower=0>[N] etaCoeff; // decay parameter
  vector<lower=0>[N] assoc1Coeff; // initial associability term
  vector<lower=0>[N] precCoeff; //participant prediction error
  
  // Group-level parameters
  vector<lower=0>[4] muGroup; //mean of alpha
  vector<lower=0>[4] varGroup; //variance of alpha
}

transformed parameters {

  vector<lower=0>[N] assoc1Subj; // subj initiial associability
  vector<lower=0>[N] etaSubj;
  vector<lower=0>[N] kappaSubj;
  vector<lower=0>[N] precSubj; // subj precCoeffision error
  
  // initialise trialwise parameters
  matrix[Tn, N] estMu; // participant estimate of mean
  matrix[Tn, N] predErr; // trialwise predicition error
  matrix[Tn, N] assoc; // time dependent associatibilty

  for (s in 1:N) {

    // Non centred paramaterisation
    assoc1Subj[s] = Phi_approx(muGroup[1] + varGroup[1] * assoc1Coeff[s]);
    etaSubj[s] = Phi_approx(muGroup[2] + varGroup[2] * etaCoeff[s]);
    kappaSubj[s] = Phi_approx(muGroup[3] + varGroup[3] * kappaCoeff[s]);
    precSubj[s] = muGroup[4] + varGroup[4] * precCoeff[s];
    
    //initialise
    estMu[1, s] = m_init[s];
    assoc[1,s] = assoc1Subj[s];
    
    for (t in 2:Tn){
      // update estimate
      estMu[t,s] = estMu[(t-1),s] + kappaSubj[s]*assoc[(t-1), s]*(seq[(t-1),s] - estMu[(t-1),s]); 
      // update associability - paper has a factor of 0.01 in first term??
      assoc[t,s] = 0.01*etaSubj[s]*fabs(seq[(t-1),s] - estMu[(t-1),s]) + (1 - etaSubj[s])*assoc[(t-1), s];
      //assoc[t,s] = etaSubjSubj[s]*fabs(seq[(t-1),s] - estMu[(t-1),s]) + (1 - etaSubj[s])*assoc[(t-1), s];
    }
  }
}

model {
  
  for (s in 1:N){
    for (t in 1:Tn){
      target += normal_lpdf(pred[t,s] | estMu[t,s], precSubj[s]);
    }
    
   /// subject level priors
   target += normal_lpdf(assoc1Coeff[s]|0, 1.0);
   target += normal_lpdf(etaCoeff[s]|0, 1.0);
   target += normal_lpdf(kappaCoeff[s]|0, 1.0);
   target += normal_lpdf(precCoeff[s]|0, 1.0);
  }
  
  //////// hyperpriors
  target += normal_lpdf(muGroup|0,1);
  //target += cauchy_lpdf(varGroup|0,10);
  //target += normal_lpdf(varGroup|1,0.2);
  target += student_t_lpdf(varGroup|3,0,1);
}

generated quantities {
  matrix[Tn, N] log_lik;
  matrix[Tn, N] y_sim;
  for (s in 1:N){
    for (t in 1:Tn) {
      y_sim[t,s] = normal_rng(estMu[t,s], precSubj[s]);
      log_lik[t,s] = normal_lpdf(pred[t,s] | estMu[t,s], precSubj[s]);
    }
  }
}


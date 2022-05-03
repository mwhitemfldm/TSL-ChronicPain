library("rstan") # observe startup messages
library("tidyverse")
library("rstanarm")
library("bayesplot")
library("loo")
library("cowplot")
library("readr")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read Data --------------------

Hseq <- select(read.csv('data/healthyPop_seqs.csv'), -X)
Hpred <- select(read.csv('data/healthyPop_preds.csv'), -X)

## Prepare data for fitting---------------------

## Healthy

#hN <- ncol(Hseq) # no. participants
hN <- 8
hTn <- 200 # no. trials
hseq=array(0,dim=c(hTn,1))
hpred=array(0,dim=c(hTn,1))
hm_init=array(0,dim=c(1))

hm_init[1]=Hpred[1,hN]

for (t in 1:hTn){
    hseq[t,1] = Hseq[t,hN]
    hpred[t,1]= Hpred[t,hN]
  }


Hdata <- list(
  N = 1,
  Tn = hTn,
  seq = hseq,
  pred = hpred,
  m_init = hm_init)

## Model fitting arguments  -----------------

# stan fit arguments
chains_value = 4          
warmup_value = 2000        
iter_value = 4000                      
cores_value = 4              

adapt_delta_value = 0.95
max_treedepth_value = 10

## Delta fit -------------------

## select pars
delta_pars = c("alphaCoeff", "alphaSubj", 
               "precCoeff", "precSubj", 
               "muGroup", "varGroup", "log_lik", "y_sim")

## fit

fit_deltaH <- stan(
  
  file = "delta2.stan",  # Stan program
  data = Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  #refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = delta_pars
  
)

## PH fit -------------

## select pars

ph_pars = c("etaCoeff", "etaSubj", 
            "kappaCoeff", "kappaSubj", 
            "assoc1Coeff", "assoc1Subj", 
            "precCoeff", "precSubj", 
            "muGroup", "varGroup", "log_lik", "y_sim")

fit_phH <- stan(
  
  file = "pearcehall2.stan",  # Stan program
  data = Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = ph_pars
  
)

## KF fit --------------

xVar1 <- 1

kalman_Hdata <- list(
  N = 1,
  Tn = hTn,
  seq = hseq,
  pred = hpred,
  m_init = hm_init,
  xVar1 = xVar1)


## select pars

kalman_pars = c("driftCoeff", "driftSubj",
                "noiseCoeff", "noiseSubj", 
                "precCoeff", "precSubj", 
                "muGroup", "varGroup", "log_lik", "y_sim")

fit_kalmanH <- stan(
  
  file = "kalman3.stan",  # Stan program
  data = kalman_Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  #refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = kalman_pars
  
)

## VKF fit ------------

## select pars
vkf_pars = c("v1Coeff", "v1Subj", 
             "w1Coeff", "w1Subj",
             "sigmaCoeff", "sigmaSubj", 
             "lambdaCoeff", "lambdaSubj",
             "precCoeff", "precSubj", 
             "muGroup", "varGroup", "log_lik", "y_sim")

fit_vkfH <- stan(
  
  file = "vkf2.stan",  # Stan program
  data = Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  #refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = vkf_pars
  
)

## Extract simulations
delta_sim <- as.matrix(summary(fit_deltaH, pars = c("y_sim"))$summary)[,1]
ph_sim <- as.matrix(summary(fit_phH, pars = c("y_sim"))$summary)[,1]
kf_sim <- as.matrix(summary(fit_kalmanH, pars = c("y_sim"))$summary)[,1]
vkf_sim <- as.matrix(summary(fit_vkfH, pars = c("y_sim"))$summary)[,1]

trials = seq(1,200)

hsim = data.frame(trials, delta_sim, ph_sim, kf_sim, vkf_sim, hseq)
# Save data
write.csv(hsim, 'simdata.csv')

library("rstan") # observe startup messages
library("tidyverse")
library("rstanarm")
library("bayesplot")
library("loo")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Fit Options -------------------

chains_value = 4          
warmup_value = 2000       
iter_value = 4000                      
cores_value = 4              

adapt_delta_value = 0.95
max_treedepth_value = 10

## Read Data --------------------

Hseq <- select(read.csv('data/healthyPop_seqs.csv'), -X)
Hpred <- select(read.csv('data/healthyPop_preds.csv'), -X)

Pseq <- select(read.csv('data/painPop_seqs.csv'), -X)
Ppred <- select(read.csv('data/painPop_preds.csv'), -X)


#hN <- ncol(Hseq) # no. participants
hN <- 1
hTn <- 200 # no. trials
hseq=array(0,dim=c(hTn,hN))
hpred=array(0,dim=c(hTn,hN))
hm_init=array(0,dim=c(hN))

for (s in 1:hN){
  hm_init[s]=Hpred[1,s]
  for (t in 1:hTn){
    hseq[t,s] = Hseq[t,s]
    hpred[t,s]= Hpred[t,s]
  }
}


## Pain

pN <- ncol(Pseq) # no. participants
#pN <- 5
pTn <- 200 # no. trials
pseq=array(0,dim=c(pTn,pN))
ppred=array(0,dim=c(pTn,pN))
pm_init=array(0,dim=c(pN))

for (s in 1:pN){
  pm_init[s]=Ppred[1,s]
  for (t in 1:pTn){
    pseq[t,s] = Pseq[t,s]
    ppred[t,s]= Ppred[t,s]
  }
}

## Prepare data ---------------------

xVar1 <- 1

kalman_Hdata <- list(
  N = hN,
  Tn = hTn,
  seq = hseq,
  pred = hpred,
  m_init = hm_init,
  xVar1 = xVar1)

kalman_Pdata <- list(
  N = pN,
  Tn = pTn,
  seq = pseq,
  pred = ppred,
  m_init = pm_init,
  xVar1 = xVar1)

## select pars

kalman_pars = c( "y_sim","muGroup", "varGroup", "log_lik")
## fitting ------------

print("fitting KALMAN")


fit_kalmanH <- stan(
  
  file = "kalman2.stan",  # Stan program
  data = kalman_Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = kalman_pars
  
)

fit_kalmanP <- stan(
  
  file = "kalman2.stan",  # Stan program
  data = kalman_Pdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = kalman_pars
  
)

## Diagnostics ------------
source("checkfit.R")

check_all_diagnostics(fit_kalmanH)
check_all_diagnostics(fit_kalmanP)

## Score model fit --------

loo_kalmanH <- loo(fit_kalmanH, save_psis = TRUE)
print(loo_kalmanH )
loo_kalmanP <- loo(fit_kalmanP, save_psis = TRUE)
print(loo_kalmanP )

## save files ----------

saveRDS(fit_kalmanP, file = "fit_kalmanP.rds")
saveRDS(loo_kalmanP, file = "loo_kalmanP.rds")
saveRDS(fit_kalmanH, file = "fit_kalmanH.rds")
saveRDS(loo_kalmanH, file = "loo_kalmanH.rds")





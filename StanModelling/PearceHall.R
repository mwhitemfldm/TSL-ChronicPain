library("rstan") # observe startup messages
library("tidyverse")
library("rstanarm")
library("bayesplot")
library("loo")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Read Data --------------------

Hseq <- select(read.csv('data/healthyPop_seqs.csv'), -X)
Hpred <- select(read.csv('data/healthyPop_preds.csv'), -X)

Pseq <- select(read.csv('data/painPop_seqs.csv'), -X)
Ppred <- select(read.csv('data/painPop_preds.csv'), -X)

## Prepare data for fitting---------------------

## Healthy

hN <- ncol(Hseq) # no. participants
#hN <- 5
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

Hdata <- list(
  N = hN,
  Tn = hTn,
  seq = hseq,
  pred = hpred,
  m_init = hm_init)

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

Pdata <- list(
  N = pN,
  Tn = pTn,
  seq = pseq,
  pred = ppred,
  m_init = pm_init)



## Model fitting arguments  -----------------

# stan fit arguments
chains_value = 4          
warmup_value = 2000        
iter_value = 4000                      
cores_value = 4              

adapt_delta_value = 0.95
max_treedepth_value = 10

## PH fit -------------

print("fitting  ph")

## select pars

ph_pars = c("etaCoeff", "etaSubj", 
            "kappaCoeff", "kappaSubj", 
            "assoc1Coeff", "assoc1Subj", 
            "precCoeff", "precSubj", 
            "muGroup", "varGroup", "log_lik")

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


fit_phP <- stan(
  
  file = "pearcehall2.stan",  # Stan program
  data = Pdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = ph_pars
  
)


## Diagnostics ------------
source("checkfit.R")

check_all_diagnostics(fit_phH)
check_all_diagnostics(fit_phP)

## Score model fit --------

loo_phH <- loo(fit_phH, save_psis = TRUE)
print(loo_phH )
loo_phP <- loo(fit_phP, save_psis = TRUE)
print(loo_phP )


## save files ----------

saveRDS(fit_phP, file = "fit_phP.rds")
saveRDS(loo_phP, file = "loo_phP.rds")
saveRDS(fit_phH, file = "fit_phH.rds")
saveRDS(loo_phH, file = "loo_phH.rds")



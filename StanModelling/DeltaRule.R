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

## Prepare data ---------------------

## Healthy data

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

delta_Hdata <- list(
  N = hN,
  Tn = hTn,
  seq = hseq,
  pred = hpred,
  m_init = hm_init)

## Pain data

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

delta_Pdata <- list(
  N = pN,
  Tn = pTn,
  seq = pseq,
  pred = ppred,
  m_init = pm_init)


## select pars
delta_pars = c("alphaCoeff", "alphaSubj", 
         "precCoeff", "precSubj", 
         "muGroup", "varGroup", "log_lik","y_sim")

## Fit models -----------------

print("fitting Delta")



fit_deltaH <- stan(
  
  file = "delta2.stan",  # Stan program
  data = delta_Hdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = delta_pars
  
  )

fit_deltaP <- stan(
  
  file = "delta2.stan",  # Stan program
  data = delta_Pdata,    # named list of data
  chains = chains_value,             # number of Markov chains
  warmup = warmup_value,        # number of warmup iterations per chain
  iter = iter_value,            # total number of iterations per chain
  cores = cores_value,              # number of cores (could use one per chain)
  refresh = 0,
  control = list(adapt_delta=adapt_delta_value,
                 max_treedepth = max_treedepth_value),
  pars = delta_pars
  
)

## Diagnostics ------------
source("checkfit.R")
print("Healthy diagnostics")
check_all_diagnostics(fit_deltaH)


print("Healthy diagnostics")
check_all_diagnostics(fit_deltaP)



## Score model fit --------
loo_deltaH <- loo(fit_deltaH, save_psis = TRUE)
print(loo_deltaH )

loo_deltaP <- loo(fit_deltaP, save_psis = TRUE)
print(loo_deltaP )


## save files ----------

saveRDS(fit_deltaP, file = "fit_deltaP.rds")
saveRDS(loo_deltaP, file = "loo_deltaP.rds")

saveRDS(fit_deltaH, file = "fit_deltaH.rds")
saveRDS(loo_deltaH, file = "loo_deltaH.rds")


## get simulated data ---------------

fit_summary <- summary(fit_deltaH, pars = c("y_sim"))$summary

simdata = fit_summary[,1]
trials = seq(1,200)
hsim = data.frame(trials, hseq, simdata)


ggplot(hsim, aes(x=trials))+
  geom_line(aes(y=hseq), color = "darkred") +
  geom_line(aes(y=simdata), color="blue")
length(hsim)

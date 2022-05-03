library("tidyverse")
library("ggpubr")
library("ggplot2")
library("cowplot")
library("readr")

## initial stuff --------------

Hconf = read.csv('data/Hconfpe.csv')
Pconf = read.csv('data/Pconfpe.csv')
p1 = ggscatter(Pconf, x = "painconf", y = "painseq", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "confidence", ylab = "sequence outcome",
               title="Confidence and sequence outcome",
               size = 0.3)
p2 = ggscatter(Pconf, x = "painconf", y = "painPE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "confidence", ylab = "prediction error",
          title="Confidence and prediction error",
          size = 0.3)

plot_grid(p1, p2,  ncol = 1, nrow =2)


h1 = ggscatter(Hconf, x = "healthconf", y = "healthseq", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "confidence", ylab = "sequence outcome",
               title="Confidence and sequence outcome",
               size = 0.3)
h2 = ggscatter(Hconf, x = "healthconf", y = "healthPE", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "confidence", ylab = "prediction error",
               title="Confidence and prediction error",
               size = 0.3)

plot_grid(h1, h2,  ncol = 1, nrow =2)



## get Healthy and Pain ids corresponding to order --------------

## Healthy ids
Hcols = colnames(read.csv('data/healthyPop_seqs.csv'))[-1]
hN = length(Hcols)
Hids = array(0,dim=c(hN))

for (s in 1:hN){
  Hids[s] = substring(Hcols[s],2)
  }
## Pain ids
Pcols = colnames(read.csv('data/painPop_seqs.csv'))[-1]
pN = length(Pcols)
Pids = array(0,dim=c(pN))

for (s in 1:pN){
  Pids[s] = substring(Pcols[s],2)
}

PROLIFIC_PID = c(Hids, Pids)

id_data = data.frame(PROLIFIC_PID)
  
## Get questionnaire Data

Hqnrs = read.csv('qnrAnalysis/Hqnrdata.csv')
Pqnrs = read.csv('qnrAnalysis/Pqnrdata.csv')
ratings = read.csv('qnrAnalysis/ratings2.csv')
qnrdata = rbind(Hqnrs, Pqnrs)

Hmcon = read.csv('Pmeanconf.csv') %>%
        select(c("prolific_id", "mean_conf")) %>%
        rename("PROLIFIC_PID" = "prolific_id")

Pmcon = read.csv('Hmeanconf.csv') %>%
        select(c("prolific_id", "mean_conf")) %>%
        rename("PROLIFIC_PID" = "prolific_id")

mcon = rbind(Hmcon,Pmcon)

totdata = inner_join(id_data, qnrdata, by="PROLIFIC_PID")
totdata = left_join(totdata, ratings, by="PROLIFIC_PID")
totdata = left_join(totdata, mcon, by="PROLIFIC_PID")
dim(totdata)
dim(qnrdata)

noise = read.csv('data/noise_data.csv')
drift = read.csv('data/drift_data.csv')

qnr_tot = data.frame(totdata, noise, drift)
write.csv(qnr_tot,'data/qnr_tot.csv')
colnames(qnr_tot)
## noise ---------------------
n1 = ggscatter(qnr_tot, x = "STAI", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "STAI score", ylab = "observation noise",
          title="Anxiety")

n2 = ggscatter(qnr_tot, x = "PHQ9", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PHQ9 score", ylab = "observation noise",
          title="Depression")


n3 = ggscatter(qnr_tot, x = "PASS20", y = "mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PASS20 score", ylab = "observation noise",
          title="Pain Anxiety")

n4 = ggscatter(qnr_tot, x = "PCS", y = "mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PCS score", ylab = "observation noise",
          title="Pain Catastrophizing")

n5 = ggscatter(qnr_tot, x = "START", y = "mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "START score", ylab = "observation noise",
          title="Back Pain")


n6 = ggscatter(qnr_tot, x = "pain", y = "mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pain rating", ylab = "observation noise",
          title="Pain rating")

n7 = ggscatter(qnr_tot, x = "fatigue", y = "mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "fatigue rating", ylab = "observation noise",
          title = "Fatigue rating")

n8 = ggscatter(qnr_tot, x = "Age", y = "mean",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Age/ yr", ylab = "observation noise",
               title = "Age")


plot_grid(n1, n2, n3, n4,  ncol = 2, nrow = 2)
plot_grid(n5, n8, n6, n7,  ncol = 2, nrow = 2)

c1 = ggscatter(qnr_tot, x = "mean_conf", y = "sd.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mean confidence rating", ylab = "drift variance posterior SD",
          title = "Drift Variance Posterior Uncertainty and Confidence")
c2 = ggscatter(qnr_tot, x = "mean_conf", y = "sd",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "mean confidence rating", ylab = "observation noise posterior SD",
               title = "Observation Noise Posterior Uncertainty and Confidence")
plot_grid(c1, c2, ncol = 1, nrow = 2)

## Drift ---------------

d1 = ggscatter(qnr_tot, x = "STAI", y = "mean.1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "STAI score", ylab = "drift variance",
          title="Anxiety")

  
d2 = ggscatter(qnr_tot, x = "PHQ9", y = "mean.1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PHQ9 score", ylab = "drift variance",
          title="Depression")

d3 = ggscatter(qnr_tot, x = "PASS20", y = "mean.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PASS20 score", ylab = "drift variance",
          title="Pain Anxiety")

d4 = ggscatter(qnr_tot, x = "PCS", y = "mean.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PCS score", ylab = "drift variance",
          title="Pain Catastophizing")

d5 = ggscatter(qnr_tot, x = "START" , y = "mean.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "START score", ylab = "drift variance",
          title="Back Pain")

d6 = ggscatter(qnr_tot, x = "pain", y = "mean.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pain rating", ylab = "drift variance", 
          title="Pain Rating")

d7 = ggscatter(qnr_tot, x = "fatigue", y = "mean.1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "fatigue rating", ylab = "drift variance",
          title="Fatigue Rating")

d8 = ggscatter(qnr_tot, x = "Age", y = "mean.1",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Age/ yr", ylab = "drift variance",
               title="Age")

plot_grid(d1, d2, d3, d4,  ncol = 2, nrow = 2)
plot_grid(d5, d8, d6, d7,  ncol = 2, nrow = 2)

d8 

## Linear regression --------------

reg2<-lm(mean.1~START+STAI+PHQ9+PASS20+PCS+Age, data = qnr_tot)
summary(reg2)

reg3<-lm(mean.1~Age, data = qnr_tot)
summary(reg3)

ggscatter(qnr_tot, x = "mean_conf", y = "sd.1",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Age/ yr", ylab = "drift variance",
               title="Age")


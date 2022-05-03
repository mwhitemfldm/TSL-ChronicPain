library("rstan") # observe startup messages
library("tidyverse")
library("rstanarm")
library("bayesplot")
library("loo")
library(hBayesDM)
library("cowplot")
library("gridExtra")

## Kalman group comparison ----------------------

fit_kalmanH = readRDS("fit_kalmanH.rds")
fit_kalmanP = readRDS("fit_kalmanP.rds")

## Mu 1 = Drift Variance ------------

pain_group =rstan::extract(fit_kalmanP, pars="muGroup[1]")[[1]]
painCI = HDIofMCMC(pain_group, credMass = 0.95)

control_group =rstan::extract(fit_kalmanH, pars="muGroup[1]")[[1]]
controlCI = HDIofMCMC(control_group, credMass = 0.95)

groupdiff = control_group - pain_group  
diffCI = HDIofMCMC(groupdiff, credMass = 0.95)
diffCI 

hdiH = plotHDI(
  sample = control_group,
  credMass = 0.95,
  xLab = "Drift variance",
  yLab = "Density",
  Title = "Control Group mean",
  binSize = 30,
) 

hdiP = plotHDI(
  sample = pain_group,
  credMass = 0.95,
  xLab = "Drift variance",
  yLab = "Density",
  Title = "Chronic pain group mean",
  binSize = 30,
) 

hdiD = plotHDI(
  sample = groupdiff ,
  credMass = 0.95,
  xLab = "Drift variance",
  yLab = "Density",
  Title = "Group mean difference",
  binSize = 30,

)


plot_grid(hdiP, hdiH, hdiD,  ncol = 1, nrow =3)


## plot both together 

group = c(rep('Pain group',length(pain_group)),
          rep('Control group', length(control_group)))

xdata = c(pain_group, control_group)

group2 = data.frame(group, xdata)





ggplot(group2, aes(x=xdata, color=group, fill=group)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, bins=60, position = "identity") +
  labs(title="Group-Level Mean Drift Variance Posteriors",x="Drift Variance", y = "Count") +
  scale_color_manual(values = c("green", "orange")) +
  scale_fill_manual(values = c("green", "orange")) +
  xlim(-0.005,0.17) + theme_classic() 
  #theme(legend.position="bottom")



## Mu 2 ----------
pain_group =rstan::extract(fit_kalmanH, pars="muGroup[2]")[[1]]
painCI = HDIofMCMC(pain_group, credMass = 0.95)

control_group =rstan::extract(fit_kalmanP, pars="muGroup[2]")[[1]]
controlCI = HDIofMCMC(control_group, credMass = 0.95)

groupdiff = control_group - pain_group 
diffCI = HDIofMCMC(groupdiff, credMass = 0.95)
diffCI 


hdiH = plotHDI(
  sample = control_group,
  credMass = 0.95,
  xLab = "Observation noise",
  yLab = "Density",
  Title = "Control group mean",
  binSize = 30,
) 

hdiP = plotHDI(
  sample = pain_group,
  credMass = 0.95,
  xLab = "Observation noise",
  yLab = "Density",
  Title = "Pain group mean",
  binSize = 30,
) 

hdiD = plotHDI(
  sample = groupdiff ,
  credMass = 0.95,
  xLab = "Observation noise",
  yLab = "Density",
  Title = "Group mean difference",
  binSize = 30,
  
)

plot_grid(hdiP, hdiH, hdiD,  ncol = 1, nrow = 3)


group = c(rep('Pain group',length(pain_group)),
          rep('Control group', length(control_group)))

xdata = c(pain_group, control_group)

group2 = data.frame(group, xdata)

ggplot(group2, aes(x=xdata, color=group, fill=group)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, bins=60, position = "identity") +
  labs(title="Group-Level Mean Observation Noise Posteriors",x="Observation noise", y = "Count") +
  scale_color_manual(values = c("green", "orange")) +
  scale_fill_manual(values = c("green", "orange")) +
  xlim(-0.005,0.21) + theme_classic() 



## group 3 -------------

pain_group =rstan::extract(fit_kalmanH, pars="muGroup[3]")[[1]]
painCI = HDIofMCMC(pain_group, credMass = 0.95)

control_group =rstan::extract(fit_kalmanP, pars="muGroup[3]")[[1]]
controlCI = HDIofMCMC(control_group, credMass = 0.95)
group = c(rep('Pain group',length(pain_group)),
          rep('Control group', length(control_group)))

xdata = c(pain_group, control_group)

group2 = data.frame(group, xdata)

ggplot(group2, aes(x=xdata, color=group, fill=group)) +
  geom_histogram( alpha=0.4, binwidth = 0.15, position = "identity") +
  labs(title="Weight histogram plot",x="muGroup[2]", y = "Count") +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_segment(aes(x = painCI[[1]], y = 330, xend = painCI[[2]], yend = 330), colour='red', size=1) +
  geom_segment(aes(x = controlCI[[1]], y = 320, xend = controlCI[[2]], yend = 320),  colour='blue', size=1) +
  xlim(-0.15, 5) + theme_classic()


plotHDI(
  sample = control_group,
  credMass = 0.95,
  xLab = "Value",
  yLab = "Density",
  Title = "Control group mean observation noise distribution",
  binSize = 30,
) 

plotHDI(
  sample = pain_group,
  credMass = 0.95,
  xLab = "Value",
  yLab = "Density",
  Title = "Pain group mean observation noise distribution",
  binSize = 30,
) 


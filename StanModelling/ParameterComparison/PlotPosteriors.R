library("rstan") # observe startup messages
library("tidyverse")
library("rstanarm")
library("bayesplot")
library("loo")
library("ggplot2")
library("cowplot")
library("readr")


source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

## Read data and extract parameters
#fit_vkfH <- readRDS("fit_delta.rds")

## Plot parameter distribution ---------------------------

## summary: get mean of params

Hlambda <- summary(fit_vkfH, pars="lambdaSubj")$summary
Plambda <- summary(fit_vkfP, pars="lambdaSubj")$summary

NsubjH <- nrow(Hlambda)
NsubjP <- nrow(Plambda)

Hlambdmat <- as.matrix(Hlambda)[,c(1,2,3,6)]
Plambdmat <- as.matrix(Plambda)[,c(1,2,3,6)]
group <- c(rep('H', NsubjH), rep('P', NsubjP))

lambda_data <- data.frame(rbind(Hlambdmat, Plambdmat)) %>%
               cbind(group)


lambda_data = read.csv("data/lambda_data.csv")
# Rainclouds for repeated measures, continued 
lambda_plot <- ggplot(lambda_data, aes(x = group, y = mean, fill = group)) +
  geom_flat_violin(aes(fill = group),position = position_nudge(x = .1, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  geom_point(aes(x = group, y = mean, colour = group),
                   position = position_jitter(width = .05), size = 1, shape = 20)  +
  geom_boxplot(aes(x = group, y = mean, fill = group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Lambda parameters")

lambda_plot

#ggsave('10repanvplot.png', width = w, height = h)
#coord_flip()+

## Plot sigma -------------------------


Hsigma <- summary(fit_vkfH, pars="sigmaSubj")$summary
Psigma <- summary(fit_vkfP, pars="sigmaSubj")$summary

NsubjH <- nrow(Hsigma)
NsubjP <- nrow(Psigma)

Hlambdmat <- as.matrix(Hsigma)[,c(1,2,3,6)]
Plambdmat <- as.matrix(Psigma)[,c(1,2,3,6)]
group <- c(rep('H', NsubjH), rep('P', NsubjP))

sigma_data <- data.frame(rbind(Hlambdmat, Plambdmat)) %>%
  cbind(group)


# Rainclouds for repeated measures, continued 
sigma_plot <- ggplot(sigma_data, aes(x = group, y = mean, fill = group)) +
  geom_flat_violin(aes(fill = group),position = position_nudge(x = .1, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  geom_point(aes(x = group, y = mean, colour = group),
             position = position_jitter(width = .05), size = 1, shape = 20)  +
  geom_boxplot(aes(x = group, y = mean, fill = group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("sigma parameters")

sigma_plot

## kalman plots --------------------------
Hdrift <- summary(fit_kalmanH, pars="driftSubj")$summary
Pdrift <- summary(fit_kalmanP, pars="driftSubj")$summary

NsubjH <- nrow(Hdrift)
NsubjP <- nrow(Pdrift)

Hdriftmat <- as.matrix(Hdrift)[,c(1,2,3,6)]
Pdriftmat <- as.matrix(Pdrift)[,c(1,2,3,6)]
group <- c(rep('Control', NsubjH), rep('Pain', NsubjP))

drift_data <- data.frame(rbind(Hdriftmat, Pdriftmat)) %>%
  cbind(group)

# Rainclouds for repeated measures, continued 
drift_plot <- ggplot(drift_data, aes(x = group, y = mean, fill = group)) +
  geom_flat_violin(aes(fill = group),position = position_nudge(x = .1, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  geom_point(aes(x = group, y = mean, colour = group),
             position = position_jitter(width = .05), size = 1, shape = 20)  +
  geom_boxplot(aes(x = group, y = mean, fill = group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Subject level drift variance")+
  xlab("Group") + ylab("subject drift variance") + theme_classic()

drift_plot


noise_data = read.csv("data/noise_data.csv")

# Rainclouds for repeated measures, continued 
noise_plot <- ggplot(noise_data, aes(x = group, y = mean, fill = group)) +
  geom_flat_violin(aes(fill = group),position = position_nudge(x = .1, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA) +
  geom_point(aes(x = group, y = mean, colour = group),
             position = position_jitter(width = .05), size = 1, shape = 20)  +
  geom_boxplot(aes(x = group, y = mean, fill = group),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Kalman filter fitted observation noise parameter")

noise_plot

## Plot posteriors -----------------

plot(fit_kalmanH, plotfun = "hist",pars="muGroup", bins=50)
plot(fit_kalmanP, plotfun = "hist",pars="muGroup", bins=50)

plot(fit_kalmanH, plotfun = "hist",pars="varGroup", bins=50)
plot(fit_kalmanP, plotfun = "hist",pars="varGroup", bins=50)

Hdrift <- summary(fit_kalmanH, pars="driftSubj")$summary[,1]
Pdrift <- summary(fit_kalmanP, pars="driftSubj")$summary[,1]

drifts <- rbind(Hdrift,Pdrift)
hist(Hdrift)

hist(Pdrift)



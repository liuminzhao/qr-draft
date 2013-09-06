#!/bin/Rscript
##' Time-stamp: <liuminzhao 09/01/2013 15:48:44>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/08/31 using bqrpt package

rm(list=ls())
library(bqrpt)
library(quantreg)
source('../R/BQRboot.R')

TOURS <- read.csv('~/Documents/qr-draft/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1/100
weight2 <- TOURS$wtkg2/100
weight3 <- TOURS$wtkg3/100
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-50)/5
race3 <- as.numeric(TOURS$RACE == 3)

y <- (weight2 - weight1)*10
n <- length(y)
X <- matrix(0, n, 4)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3
X[,4] <- weight1

dat <- data.frame(weight1, weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)

## PLOT
library(ggplot2)
library(gridExtra)
ggplot(data = dat, aes(x = age_center, y = y)) + geom_point()

ggplot(data = dat, aes(x = weight1, y = y)) + geom_point()

ggplot(data = dat, aes(x = race, y = y)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) at 6 months') + xlab('Race')

ggplot(data = dat, aes(x = age_center, y = y, color = race)) + geom_point()

box3 <- ggplot(data = dat, aes(x = race, y = y)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) at 6 months') + xlab('Race')

scat <- ggplot(data = dat, aes(x = age_center, y = y, color = race)) + geom_point()

pdf('weight-age-race-base.pdf')
sds <- grid.arrange(box1, box2, box3, box4, box5, box6, nrow = 2, ncol = 3)
dev.off()

###############
## ANALYSIS
###############

mcmc <- list(nburn=30000, nskip=2, nsave=10000, ndisp=10000, arate=0.2)
## mcmc <- list(nburn=30000, nskip=1, nsave=30000, ndisp=10000, arate=0.25)
quan <- c(0.1, 0.3, 0.5, 0.7, 0.9)

modss <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'ss')

mod <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'normal')

summary(mod)
summary(modss)

## RQ

modrq <- rq(y ~ age_center + race3, tau = c(0.1,0.3,0.5,0.7,0.9))

## BQR
modbqr1 <- BayesQReg(y, X, 0.1)
modbqr3 <- BayesQReg(y, X, 0.3)
modbqr5 <- BayesQReg(y, X, 0.5)
modbqr7 <- BayesQReg(y, X, 0.7)
modbqr9 <- BayesQReg(y, X, 0.9)

BQRCoef(modbqr1)
BQRCI(modbqr1)

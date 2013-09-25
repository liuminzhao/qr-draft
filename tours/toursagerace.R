#!/bin/Rscript
##' Time-stamp: <liuminzhao 09/25/2013 16:28:11>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/08/31 using bqrpt package

rm(list=ls())
library(bqrpt)
library(quantreg)
source('../R/BQRboot.R')
set.seed(1)

TOURS <- read.csv('~/Documents/qr-draft/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-50)/5
race3 <- as.numeric(TOURS$RACE == 3)

y <- weight1 - weight2
n <- length(y)
X <- matrix(0, n, 3)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3

dat <- data.frame(loss = y, weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)


###############
## PLOT
###############
library(ggplot2)
library(gridExtra)

box1 <- ggplot(data = dat, aes(x = age, y = loss)) + geom_point() + ylab('Weight Loss (Kg)') + xlab('Age') + ggtitle("Weight Loss vs Age")

box2 <- ggplot(data = dat, aes(x = race, y = loss)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight Loss (Kg)') + xlab('Race') + ggtitle("Weight Loss vs Race")

box3 <- ggplot(data = dat, aes(x = age, y = loss, color = race)) + geom_point() + ylab('Weight Loss (Kg)') + xlab('Age') + ggtitle("Weight Loss vs Age") + ggtitle('Weight Loss vs Age and Race')

pdf('../image/weight-age-race.pdf', width = 20, height = 7)
sds <- grid.arrange(box1, box2, box3, nrow = 1, ncol = 3)
dev.off()

###############
## ANALYSIS
###############
p <- dim(X)[2]
tuneinit <- c(rep(0.3,p),1, c(0.3, p-1), 0.04, 0.1)
mcmc <- list(nburn=30000, nskip=5, nsave=30000, ndisp=10000, arate=0.2, tuneinit = tun\
einit)
quan <- c(0.1, 0.3, 0.5, 0.7, 0.9)

modss <- HeterPTlmMH(y, X, mcmc = mcmc, quan = quan, den = TRUE, method = 'ss')

mod <- HeterPTlmMH(y, X, mcmc = mcmc, quan = quan, den = TRUE, method = 'normal')

summary(mod)
summary(modss)

coefmodss <- c(t(coef(modss)$betatauMean))
coefmod <- c(t(coef(mod)$betatauMean))

cimodsslbd <- c(t(coef(modss)$betatauCIlbd))
cimodssubd <- c(t(coef(modss)$betatauCIubd))

cimodlbd <- c(t(coef(mod)$betatauCIlbd))
cimodubd <- c(t(coef(mod)$betatauCIubd))

cimodss <- cbind(coefmodss, cimodsslbd, cimodssubd)
cimod <- cbind(coefmod, cimodlbd, cimodubd)

## RQ

modrq <- rq(y ~ age_center + race3, tau = c(0.1,0.3,0.5,0.7,0.9))

coefmodrq <- t(coef(modrq))

cimodrq <- read.table('rq.txt', header = F)

## BQR
modbqr1 <- BayesQReg(y, X, 0.1)
modbqr3 <- BayesQReg(y, X, 0.3)
modbqr5 <- BayesQReg(y, X, 0.5)
modbqr7 <- BayesQReg(y, X, 0.7)
modbqr9 <- BayesQReg(y, X, 0.9)

coefmodbqr <- c(BQRCoef(modbqr1), BQRCoef(modbqr3), BQRCoef(modbqr5),
                BQRCoef(modbqr7), BQRCoef(modbqr9))

cimodbqrlbd <- c(BQRCIlbd(modbqr1),
                 BQRCIlbd(modbqr3),
                 BQRCIlbd(modbqr5),
                 BQRCIlbd(modbqr7),
                 BQRCIlbd(modbqr9))

cimodbqrubd <- c(BQRCIubd(modbqr1),
                 BQRCIubd(modbqr3),
                 BQRCIubd(modbqr5),
                 BQRCIubd(modbqr7),
                 BQRCIubd(modbqr9))

cimodbqr <- cbind(coefmodbqr, cimodbqrlbd, cimodbqrubd)

## put together

totalci <- cbind(cimod, cimodss, cimodrq, cimodbqr)

## posterior prob of coefficient are zero for each quantile
zero <- rbind(coef(modss)$deltabetaprop, coef(modss)$deltagammaprop)

library(xtable)
sink('tours-result-0925.txt')
print(xtable(totalci))
print(xtable(zero))
sink()

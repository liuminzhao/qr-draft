#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/31/2013 14:58:01>
##' manipulate data TOURS
##' 2013/06/05 focus on AGE and RACE
##' 2013/06/22 add baseline y0 as a covariate
##' 2013/08/31 using bqrpt package

rm(list=ls())
library(bqrpt)

TOURS <- read.csv('~/Documents/qr-draft/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-50)/5
race3 <- as.numeric(TOURS$RACE == 3)

y <- weight2 - weight1
n <- length(y)
X <- matrix(0, n, 4)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3
X[,4] <- weight1

dat <- data.frame(weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)

## PLOT
## library(ggplot2)
## library(gridExtra)
## ggplot(data = dat, aes(x = age_center, y = weight2)) + geom_point()
## ggplot(data = dat, aes(x = age_center, y = weight3)) + geom_point()

## box1 <- ggplot(data = dat, aes(x = age, y = weight2*100)) + geom_point()+ ylab('Weight (Kg) at 6 months') + xlab('Age')
## box2 <- ggplot(data = dat, aes(x = age, y = weight3*100)) + geom_point()+ ylab('Weight (Kg) at 18 months') + xlab('Age')
## box3 <- ggplot(data = dat, aes(x = race, y = weight2*100)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight (Kg) at 6 months') + xlab('Race')
## box4 <- ggplot(data = dat, aes(race, y = weight3*100)) + geom_boxplot() + ylab('Weight (Kg) at 18 months') + xlab('Race') + scale_x_discrete(labels=c("Black", "White"))
## box5 <- ggplot(data = dat, aes(x = weight1*100, y = weight2*100)) + geom_point() + ylab('Weight (Kg) at 6 months') + xlab('Weight (Kg) at baseline')
## box6 <- ggplot(data = dat, aes(x = weight1*100, y = weight3*100)) + geom_point() + ylab('Weight (Kg) at 18 months') + xlab('Weight (Kg) at baseline')

## scat1 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat2 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()
## scat3 <- ggplot(data = dat, aes(x = age_center, y = weight2, color = race)) + geom_point()
## scat4 <- ggplot(data = dat, aes(x = age_center, y = weight3, color = race)) + geom_point()

## pdf('weight-age-race-base.pdf')
## sds <- grid.arrange(box1, box2, box3, box4, box5, box6, nrow = 2, ncol = 3)
## dev.off()

###############
## ANALYSIS
###############

mcmc <- list(nburn=10000, nskip=1, nsave=10000, ndisp=10000, arate=0.25)
quan <- c(0.1, 0.3, 0.5, 0.7, 0.9)

modss <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'ss')

mod <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'normal')

summary(mod)
summary(modss)
## age race

X <- X[, 1:3]

mod3ss <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'ss')

mod3 <- HeterPTlm(y, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), den = TRUE, method = 'normal')

summary(mod3)
summary(mod3ss)

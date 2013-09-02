#!/bin/Rscript
##' Time-stamp: <liuminzhao 09/02/2013 14:08:41>
##' 2013/08/31 simulation M3

sink('sim-m3-0831.txt')
rm(list = ls())
library(bqrpt)
library(quantreg)
library(xtable)
source('sendEmail.R')
source('RQboot.R')
source('BQRboot.R')
library(doMC)
registerDoMC()
options(cores=10)
set.seed(1)

rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}

###############
## PARAMETERS
###############
n <- 200
mcmc <- list(nburn=20000, nskip=1, nsave=20000, ndisp=20000, arate=0.25)
b1 <- b2 <- 1
quan <- c(0.5, 0.9)
###############
## SIMULATION
###############

boot <- 100
start <- proc.time()[3]

result <- foreach(icount(boot), .combine=rbind) %dopar% {

  x1 <- runif(n)
  x2 <- runif(n)
  e1 <- rMN3(n)

  y1 <- 1 + x1*b1 + x2*b2 + e1

  X <- cbind(1,x1,x2)

  ## rq
  modrq5 <- rq(y1 ~ x1 + x2, 0.5)
  modrq9 <- rq(y1 ~ x1 + x2, 0.9)

  foo1.9 <- summary(rq(y1~x1+x2, 0.9))

  ## bqr
  modbqr5 <- BayesQReg(y1, X, 0.5)
  modbqr9 <- BayesQReg(y1, X, 0.9)

  ## pt
  modpt <- HeterPTlm(y1, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6))

  ## pt ss
  modptss <- HeterPTlm(y1, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), method = 'ss')

  ## coef
  coefrq5 <- coef(modrq5)
  coefrq9 <- coef(modrq9)

  coefbqr5 <- BQRCoef(modbqr5)
  coefbqr9 <- BQRCoef(modbqr9)

  coefpt5 <- coef(modpt)$betatauMean[1, ]
  coefpt9 <- coef(modpt)$betatauMean[2, ]

  coefptss5 <- coef(modptss)$betatauMean[1, ]
  coefptss9 <- coef(modptss)$betatauMean[2, ]

  ans <- c(coefrq5, coefrq9,
           coefbqr5, coefbqr9,
           coefpt5, coefpt9,
           coefptss5, coefptss9)
}

write.table(result, file="sim-m3-result-0831.txt", row.names = F, col.names = F)
sendEmail(subject = "simulation-m3", text = "done", address = "liuminzhao@gmail.com")



###############
## TRUE VALUE
###############
result <- read.table('sim-m3-result-0831.txt')
truebetatau5 <- c(1,1,1)
truebetatau9 <- c(1+2 + qnorm(0.8),1,1)
truebetatau <- rep(c(truebetatau5, truebetatau9), 4)

library(xtable)

## MSE
mse <- rep(0, 24)
for (i in 1:24){
  mse[i] <- mean((result[,i] - truebetatau[i])^2)*100
}
mse <- matrix(mse, 6, 4)
colnames(mse) <- c('RQ', 'BQR', 'PT', 'PTSS')
print(xtable(mse))
print(mse)

## BIAS
bias <- rep(0, 24)
for (i in 1:24){
  bias[i] <- mean((result[,i] - truebetatau[i]))*100
}
bias <- matrix(bias, 6, 4)
colnames(bias) <- c('RQ', 'BQR', 'PT', 'PTSS')
print(xtable(bias))
print(bias)

## time
cat("Time: ", proc.time()[3] - start, '\n')
sink()

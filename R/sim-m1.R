#!/bin/Rscript
##' Time-stamp: <liuminzhao 08/31/2013 00:56:13>
##' 2013/08/31 simulation M1

sink('sim-m1-0831.txt')
rm(list = ls())
library(bqrpt)
library(quantreg)
source('sendEmail.R')
source('RQboot.R')
source('BQRboot.R')
library(doMC)
registerDoMC()
options(cores=10)
set.seed(1)


###############
## PARAMETERS
###############
n <- 200
mcmc <- list(nburn=10000, nskip=1, nsave=10000, ndisp=10000, arate=0.25)
b1 <- b2 <- 1
quan <- c(0.5, 0.9)
###############
## SIMULATION
###############

boot <- 100

result <- foreach(icount(boot), .combine=cbind) %dopar% {

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  e1 <- rnorm(n)

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

  ans <- c(coefrq5, coefbqr5, coefpt5, coefptss5,
           coefrq9, coefbqr9, coefpt9, coefptss9)
}

write.table(result, file="sim-m1-result-0831.txt", row.names = F, col.names = F)
sendEmail(subject = "simulation-m1", text = "done", address = "liuminzhao@gmail.com")

truebetatau5 <- c(1,1,1)
truebetatau9 <- c(1+qnorm(0.9),1,1)
truebetatau1 <- rep(c(truebetatau5, truebetatau9), each = 4)

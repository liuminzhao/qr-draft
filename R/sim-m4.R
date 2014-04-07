#!/bin/Rscript
##' Time-stamp: <liuminzhao 04/06/2014 22:10:01>
##' 2013/08/31 simulation M4

sink('sim-m4-0406.txt')
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


rMN <-function(n){
  posneg<-rbinom(n,1,0.8)
  posneg*rnorm(n) + (1-posneg)*rnorm(n, 3, sqrt(3))
}


###############
## PARAMETERS
###############
n <- 200
tuneinit <- c(0.3, 0.3, 1, 0.3, 0.04, 0.1)
mcmc <- list(nburn=0, nskip=5, nsave=30000, ndisp=30000, arate=0.2, tuneinit = tuneinit)

b1 <- 1
quan <- c(0.5, 0.9)
###############
## SIMULATION
###############

boot <- 100
start <- proc.time()[3]

result <- foreach(icount(boot), .combine=rbind) %dopar% {

  x1 <- runif(n, min = -1, max = 1)
  e1 <- rMN(n)

  X <- cbind(1,x1)

  y1 <- 1 + x1*b1 + e1

  ## rq
  modrq5 <- rq(y1 ~ x1, 0.5)
  modrq9 <- rq(y1 ~ x1, 0.9)

  ## bqr
  modbqr5 <- BayesQReg(y1, X, 0.5)
  modbqr9 <- BayesQReg(y1, X, 0.9)

  ## pt
  modpt <- HeterPTlmMH(y1, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6))

  ## pt ss
  modptss <- HeterPTlmMH(y1, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6), method = 'ss')

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

write.table(result, file="sim-m4-result-0406.txt", row.names = F, col.names = F)
sendEmail(subject = "simulation-m4", text = "done", address = "liuminzhao@gmail.com")



###############
## TRUE VALUE
###############
result <- read.table('sim-m4-result-0406.txt')
pMN.5 <- function(x){
  0.8*pnorm(x)+0.2*pnorm(x,3,sqrt(3))-0.5
}

pMN.9 <- function(x){
  0.8*pnorm(x)+0.2*pnorm(x,3,sqrt(3))-0.9
}

q5 <- uniroot(pMN.5, c(-5,5))$root
q9 <- uniroot(pMN.9, c(-5,5))$root

truebetatau5 <- c(1,1) + c(1, 0)*q5
truebetatau9 <- c(1,1) + c(1, 0)*q9
truebetatau <- rep(c(truebetatau5, truebetatau9), 4)

library(xtable)

## MSE
mse <- rep(0, 16)
for (i in 1:16){
  mse[i] <- mean((result[,i] - truebetatau[i])^2)*100
}
mse <- matrix(mse, 4, 4)
colnames(mse) <- c('RQ', 'BQR', 'PT', 'PTSS')
print(xtable(mse))
print(mse)
## MCSE
mcse <- rep(0, 16)
for (i in 1:16){
  mcse[i] <- sd((result[,i] - truebetatau[i])^2)/sqrt(boot) * 100
}
mcse <- matrix(mcse, 4, 4)
colnames(mcse) <- c('RQ', 'BQR', 'PT', 'PTSS')
print(xtable(mcse))
print(mcse)

## combine mse and MCSE
msemcse <- matrix(0, 4, 8)
for (i in 1:4){
  msemcse[, i*2 - 1] <- mse[, i]
  msemcse[, i*2] <- mcse[, i]
}
print(msemcse)
print(xtable(msemcse))

## BIAS
bias <- rep(0, 16)
for (i in 1:16){
  bias[i] <- mean((result[,i] - truebetatau[i]))*100
}
bias <- matrix(bias, 4, 4)
colnames(bias) <- c('RQ', 'BQR', 'PT', 'PTSS')
print(xtable(bias))
print(bias)


## time
cat("Time: ", proc.time()[3] - start, '\n')
sink()

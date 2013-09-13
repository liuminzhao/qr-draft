#!/bin/Rscript
##' Time-stamp: <liuminzhao 09/13/2013 15:22:38>
##' 2013/08/31 simulation M1
##' 2013/09/03 new

sink('sim-0913.txt')
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


###############
## PARAMETERS
###############
n <- 500
mcmc <- list(nburn=20000, nskip=1, nsave=20000, ndisp=30000, arate=0.4)
b1 <- 1
quan <- c(0.5,  0.875)
###############
## SIMULATION
###############

r3MN <- function(n){
  ind <- sample(4, n, replace = TRUE)
  return(I(ind == 1)*rnorm(n, -4, 1) + I(ind == 2)*rnorm(n) + I(ind == 3)*rnorm(n) + I(ind == 4)*rnorm(n, 4, 1))
}

p3MN <- function(x, tau){
  0.25*pnorm(x, -4, 1) + 0.5*pnorm(x) + 0.25*pnorm(x, 4, 1) - tau
}

q5 <- uniroot(p3MN, c(-10, 10), tau = 0.5)$root
q875 <- uniroot(p3MN, c(-10, 10), tau = 0.875)$root

boot <- 100
start <- proc.time()[3]

result <- foreach(icount(boot), .combine=rbind) %dopar% {

  x1 <- runif(n, max = 4)
  e1 <- r3MN(n)

  y1 <- 1 + x1*b1 + e1

  X <- cbind(1,x1)

  ## rq
  modrq5 <- rq(y1 ~ x1, 0.5)
  modrq9 <- rq(y1 ~ x1, 0.875)

  ## bqr
  modbqr5 <- BayesQReg(y1, X, 0.5)
  modbqr9 <- BayesQReg(y1, X, 0.875)

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

write.table(result, file="sim-result-0913.txt", row.names = F, col.names = F)
sendEmail(subject = "simulation", text = "done", address = "liuminzhao@gmail.com")



###############
## TRUE VALUE
###############
result <- read.table('sim-result-0913.txt')
truebetatau5 <- c(1 + q5, 1)
truebetatau9 <- c(1 + q9, 1)
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

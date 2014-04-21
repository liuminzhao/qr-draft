#!/bin/Rscript
##' Time-stamp: <liuminzhao 04/20/2014 18:59:47>
##' 2013/08/31 simulation M1
##' 2013/09/03 new

sink('sim-m5-0420.txt')
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
n <- 200
tuneinit <- c(0.3, 0.3, 1, 0.3, 0.3, 0.3)
mcmc <- list(nburn = 10000, nskip=5, nsave=30000, ndisp=30000, arate=0.2, tuneinit = tuneinit)
b1 <- 1
quan <- c(0.5,  0.9)
p <- 0.5
alpha <- 0

###############
## SIMULATION
###############

boot <- 100
start <- proc.time()[3]

result <- foreach(icount(boot), .combine=rbind) %dopar% {
  R <- rbinom(n, 1, p)
  x1 <- runif(n, min = -1, max = 1)
  y1 <- rep(n, 0)
  for (i in 1:n){
    if (R[i] == 1){
      y1[i] <- 2 + x1[i] + (1 + alpha*x1[i])*rnorm(1)
    } else {
      y1[i] <- -2 - x1[i] + (1 + alpha*x1[i])*rnorm(1)
    }
  }

  X <- cbind(1,x1)

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

write.table(result, file="sim-m5-result-0420.txt", row.names = F, col.names = F)
sendEmail(subject = "simulation", text = "done", address = "liuminzhao@gmail.com")



###############
## TRUE VALUE
###############
result <- read.table('sim-m5-result-0420.txt')

quan1 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 2+x, 1 + alpha*x) - .5*pnorm(y, -2-x, 1 + alpha*x))
}

SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}

xsim <- seq(-1, 1, len = 100)
y11 <- sapply(xsim, function(x) SolveQuan1(x, 0.1))
y13 <- sapply(xsim, function(x) SolveQuan1(x, 0.3))
y15 <- sapply(xsim, function(x) SolveQuan1(x, 0.5))
y17 <- sapply(xsim, function(x) SolveQuan1(x, 0.7))
y19 <- sapply(xsim, function(x) SolveQuan1(x, 0.9))

q11 <- lm(y11~xsim)$coef
q13 <- lm(y13~xsim)$coef
q15 <- lm(y15~xsim)$coef
q17 <- lm(y17~xsim)$coef
q19 <- lm(y19~xsim)$coef

truebetatau5 <- q15
truebetatau9 <- q19
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

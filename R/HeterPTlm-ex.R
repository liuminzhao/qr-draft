#!/bin/Rscript
rm(list=ls(all=TRUE))

set.seed(1)
source('HeterPTlm.R')
rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}

n <- 100
p <- 3

x1 <- rnorm(n)
x2 <- x1^2
X <- cbind(1,x1,x2)

e1 <- rnorm(n)
e2 <- rgamma(n,3,1)
e3 <- rMN3(n)

## plot(density(e1))
## plot(density(e2))
## plot(density(e3))

y1 <- X[,1]+X[,2]+X[,3]+e1
y2 <- X[,1]+X[,2]+X[,3]+e2
y3 <- X[,1]+X[,2]+X[,3]+e3
y4 <- X[,1]+X[,2]+X[,3]+(1-0.5*X[,2]+0.5*X[,3])*e3
y5 <- X[,1]+X[,2]+X[,3]+(1-0.5*X[,2]+0.5*X[,3])*e2

mod3 <- lm(y3 ~ X[,2] + X[,3])

betapm <- coef(mod3)
betapv <- vcov(mod3)
gammapm <- rep(0, p)
gammapv <- diag(p)*100
sigmap <- c(2,1/2)
#sigmap <- c(1,1) # for normal candidate and gamma prior for sigma and alpha
a0b0 <- c(1,1)

# change prior, reduce noise

prior <- list(betapm=betapm, betapv=betapv,
              gammapm=gammapm, gammapv=gammapv,
              sigmap=sigmap, a0b0=a0b0)

# mcmc <- list(nburn=20000, nskip=20, nsave=10000, ndisp=1000, arate=0.44)
mcmc <- list(nburn=20000, nskip=10, nsave=10000, ndisp=100, arate=0.44)

quan <- c(0.5, 0.9)

## foo1 <- HeterPTlm(y1, X, mcmc, prior, quan)
## foo1median <- HeterPTlm(y1, X, mcmc, prior, quan, method = 'median')
foo1ss <- HeterPTlm(y3, X, mcmc, prior, quan, method = 'ss')
## foo2 <- HeterPTlm(y2, X, mcmc, prior, quan)
## foo3 <- HeterPTlm(y3, X, mcmc, prior, quan)
## foo4 <- HeterPTlm(y4, X, mcmc, prior, quan)
## foo5 <- HeterPTlm(y5, X, mcmc, prior, quan)

## foo1.2 <- HeterPTlm(y1, X, mcmc, prior, quan)
## foo2.2 <- HeterPTlm(y2, X, mcmc, prior, quan)
## foo3.2 <- HeterPTlm(y3, X, mcmc, prior, quan)
## foo4.2 <- HeterPTlm(y4, X, mcmc, prior, quan)
## foo5.2 <- HeterPTlm(y5, X, mcmc, prior, quan)

print(summary(foo1ss))

## summary(foo3)
## summary(foo4)

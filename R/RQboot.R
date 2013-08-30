######################################################
# 2012/02/26 RQ boot, w/ betatau, credible interval,
# length of interval, coverage rate
# quantreg, rq function
# M1-M5
######################################################
library(quantreg)

RQBootSummary <- function(foo, truebetatau){
  # est
  beta <- foo$coef[,1]

  # lbd
  lbd <- foo$coef[,2]
  
  # ubd
  ubd <- foo$coef[,3]
  
  # length

  len <- ubd-lbd

  # mse

  mse <- beta - truebetatau

  # cover
  p <- length(truebetatau)
  for (j in 2:p)   cover <- prod((truebetatau[j]>lbd[j] && truebetatau[j]<ubd[j]))

  return(list(beta=beta, lbd=lbd, ubd=ubd, length=len, mse=mse, cover=cover))
}

rASL <- function(n,lambda,tau){
  posneg<-rbinom(n,1,tau)
  (1-posneg)*rexp(n,tau*lambda)-posneg*rexp(n,(1-tau)*lambda)
}

rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}


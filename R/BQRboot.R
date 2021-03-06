####################################################
# Time-stamp: <liuminzhao 09/02/2013 00:06:14>
# 2012/02/27 Reich simulation BQRiid on M1.5 M1.9
####################################################

source('../R/BQRiid.R')
BQR.Summary <- function(foo, truebetatau){

  burn <- foo$burn
  runs <- foo$runs
  beta.coef <- apply(foo$beta[burn:runs,],2,median)
  lbd <- apply(foo$beta[burn:runs,], 2, function(x) quantile(x, 0.025))
  ubd <- apply(foo$beta[burn:runs,], 2, function(x) quantile(x, 0.975))
  len <- ubd-lbd
  mse <- beta.coef-truebetatau
  p <- length(truebetatau)
  for (j in 2:p)   cover <- prod((truebetatau[j]>lbd[j] && truebetatau[j]<ubd[j]))

  return(list(beta=beta.coef, lbd=lbd, ubd=ubd, length=len, mse=mse, cover=cover))

}

BQRCoef <- function(foo){
  burn <- foo$burn
  runs <- foo$runs
  return(apply(foo$beta[burn:runs, ], 2, median))
}

BQRCIlbd <- function(foo){
  burn <- foo$burn
  runs <- foo$runs
  return(apply(foo$beta[burn:runs, ], 2, function(x) quantile(x, probs = c(0.025))))
}

BQRCIubd <- function(foo){
  burn <- foo$burn
  runs <- foo$runs
  return(apply(foo$beta[burn:runs, ], 2, function(x) quantile(x, probs = c(0.975))))
}

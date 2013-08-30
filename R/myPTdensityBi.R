rm(list=ls())
myPTdensityBi <- function(y, mcmc,   maxm = floor(log(length(y))/log(2))){
  dyn.load("myptdensitybi.so")
  n <- length(y)
  maxm <- maxm
  nsub <- dim(y)[1]
  q <- dim(y)[2]
  mdzero <- 1
  ngrid <- 200
  f <- matrix(0, ngrid, 2)
  left <- min(c(y)) - 0.5*sd(c(y))
  right <- max(c(y)) + 0.5*sd(c(y))
  grid <- seq(left, right, length=ngrid)
  arate <- 0.25
  mu <- rep(0, 2)
  Sigma <- diag(2)
  alpha <- 1
  nburn <- mcmc[1]
  nsave <- mcmc[2]
  nskip <- mcmc[3]
  sigmasave <- matrix(0, nsave, 3)
  alphasave <- rep(0, nsave)
  musave <- matrix(0, nsave, 2)
  fsave <- array(0, c(nsave, ngrid, 2))
  ratesave <- tunesave <- matrix(0, nburn, 7)
  
  foo <- .Fortran("myptdensitybi",
                  maxm  = as.integer(maxm),
                  nsub  = as.integer(nsub),
                  q     = as.integer(q),
                  y     = as.double(y),
                  mu    = as.double(mu),
                  Sigma = as.double(Sigma),
                  nburn = as.integer(nburn),
                  nsave = as.integer(nsave),
                  nskip = as.integer(nskip),
                  arate = as.double(arate),
                  grid  = as.double(grid),
                  ngrid = as.integer(ngrid),
                  f     = as.double(f),
                  alpha = as.double(alpha),
                  musave = as.double(musave),
                  sigmasave = as.double(sigmasave),
                  alphasave = as.double(alphasave),
                  fsave = as.double(fsave),
                  ratesave = as.double(ratesave),
                  tunesave = as.double(tunesave)
                  )

  musave <- matrix(foo$musave, nsave, 2)
  sigmasave <- matrix(foo$sigmasave, nsave, 3)

  ratesave <- matrix(foo$ratesave, nburn, 7)
  tunesave <- matrix(foo$tunesave, nburn, 7)
  dens <- matrix(foo$f, ngrid, 2)

  return(list(foo=foo, ratesave=ratesave, tunesave=tunesave, dens=dens,
              musave = musave, sigmasave = sigmasave,
              alphasave = foo$alphasave))
}

rho <- 0.5 

nrec <- 200
nsub <- 100
p <- 3
q <- 2


rho <- -0.5
Sigma <- matrix(c(1, rho, rho,1),2,2)

e1 <- matrix(rnorm(nsub*q), nsub, q)
e2.matrix <- kronecker(diag(nsub),t(chol(Sigma)))%*%rnorm(nrec)
e2 <- matrix(e2.matrix, nsub, q, byrow=T)
e3.matrix <- matrix(0, q, nsub)
mypi <- rbinom(nsub, size=1, prob=0.5)
mnSigma <- matrix(c(0.15, 0.02, 0.02, 0.04),2,2)
mn1 <- matrix(c(-1.35, 0.28),2,1)
mn2 <- matrix(c(1.35, 0.28),2,1)
for (i in 1:nsub){
  if (mypi[i]==1)
    e3.matrix[,i] <- mn1+t(chol(mnSigma))%*%rnorm(q)
  else
    e3.matrix[,i] <- mn2+t(chol(mnSigma))%*%rnorm(q)
}
e3 <- t(e3.matrix)

mcmc <- c(20000, 10000, 10)

mcmc <- c(20000, 20000, 20)

a <- myPTdensityBi(e1, mcmc)
a2 <- myPTdensityBi(e2, mcmc)
a3 <- myPTdensityBi(e3, mcmc)

plot(a$foo$grid, a$dens[,1] , 'l')
plot(a$foo$grid, a$dens[,2] , 'l')

plot(a2$foo$grid, a2$dens[,1] , 'l')
plot(a2$foo$grid, a2$dens[,2] , 'l')

plot(a3$foo$grid, a3$dens[,1] , 'l')
plot(a3$foo$grid, a3$dens[,2] , 'l')

Diagnose <- function(mod, ask = FALSE){
  par(mfrow = c(2, 2), ask = ask)
  for (i in 1:7) {
    plot(mod$ratesave[, i])
    plot(mod$tunesave[, i])
  }
  plot(ts(mod$musave[, 1]))
  plot(ts(mod$musave[, 2]))
  plot(ts(mod$sigmasave[, 1]))
  plot(ts(mod$sigmasave[, 2]))
  plot(ts(mod$sigmasave[, 3]))
  plot(ts(mod$alphasave))

  plot(mod$foo$grid, mod$dens[,1] , 'l')
  plot(mod$foo$grid, mod$dens[,2] , 'l')


  acf(mod$musave)
  acf(mod$sigmasave)
  acf(mod$alphasave)

  
}

Diagnose(a, ask = TRUE)
Diagnose(a2, ask = FALSE)
Diagnose(a3, ask = TRUE)

acf(a3$musave)
acf(a3$sigmasave)
acf(a3$alphasave)

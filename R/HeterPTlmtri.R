##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multivariate Bayesian Linear quantile regression with
##' Mixture of Polya tree priors , here d=3
##' @param y 
##' @param X 
##' @param nsub 
##' @param mcmc 
##' @param prior 
##' @param quan 
##' @return class HeterPTlmtri
##' @author Minzhao Liu
HeterPTlmtri <- function(y, X, nsub, mcmc, prior, quan){

  dyn.load('heterptlmtri.so')
  
  # DATA
  nrec <- prod(dim(y))
  q <- nrec/nsub
  p <- dim(X)[2]

  # PT
  if (is.null(prior$maxm)) {
    maxm <- floor(log(nsub)/log(2)) }  else maxm <- prior$maxm
  if (is.null(prior$mdzero)) {mdzero <- 1}  else  mdzero <- prior$mdzero

  # PRIOR
  betapm <- prior$betapm
  betapv <- prior$betapv
  gammapm <- prior$gammapm
  gammapv <- prior$gammapv
  tau <- prior$tau
  a0b0 <- prior$a0b0

  # MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp
  arate <- mcmc$arate

  # QUAN
  quan <- quan
  nquan <- length(quan)
  
  # SAVE
  betasave<- array(0, c(p,q, nsave))
  gammasave <- array(0, c( p, q, nsave))
  sigmasave <- matrix(0, nsave, 6) # sigma1, sigma2, sigma3, rho1, rho2, rho3
  alphasave <- rep(0, nsave)
  quansave <- matrix(0, nsave, nquan*q)

  tunesave <- aratesave <- matrix(0, nburn/100, 2*q+8)
  
  # GRID
  ngrid <- 200
  grid <- seq(-5,5, length=ngrid)
  f <- matrix(0, ngrid, q)
  
  # INITIAL
  beta <- solve(t(X)%*%X)%*%t(X)%*%y
  gamma <- matrix(0, p, q)
  gamma[1,] <- 1
  sigmavec <- c(1, 1, 1,0,0,0)
  Sigma <- diag(3)
  alpha <- 1
  b <- (y-X%*%beta)/(X%*%gamma)
  propv <- solve(t(X)%*%X)
  # WORKING
  whicho <- whichn <- rep(0, nsub)
  
  # DEBUG, TUNING
   # not used yet

  ####################################
  # MCMC USING FORTRAN

  foo <- .Fortran("heterptlmtri",
                  maxm=as.integer(maxm),
                  mdzero=as.integer(mdzero),
                  nrec=as.integer(nrec),
                  nsub=as.integer(nsub),
                  p=as.integer(p),
                  q=as.integer(q),
                  x=as.double(X),
                  y=as.double(y),
                  betapm=as.double(betapm),
                  betapv=as.double(betapv),
                  tau=as.double(tau),
                  a0b0=as.double(a0b0),
                  gammapm=as.double(gammapm),
                  gammapv=as.double(gammapv),
                  nburn=as.integer(nburn),
                  nskip=as.integer(nskip),
                  nsave=as.integer(nsave),
                  ndisp=as.integer(ndisp),
                  betasave=as.double(betasave),
                  gammasave=as.double(gammasave),
                  sigmasave=as.double(sigmasave),
                  alphasave=as.double(alphasave),
                  beta=as.double(beta),
                  gamma=as.double(gamma),
                  alpha=as.double(alpha),
                  Sigma=as.double(Sigma),
                  propv=as.double(propv),
                  arate=as.double(arate),
                  ngrid=as.integer(ngrid),
                  grid=as.double(grid),
                  f=as.double(f),
                  nquan=as.integer(nquan),
                  qtile=as.double(quan),
                  quansave=as.double(quansave)
#                  aratesave=as.double(aratesave),
#                  tunesave=as.double(tunesave)
              )

  ####################################

  betasave <- array(foo$betasave, dim=c(p, q, nsave))
  gammasave <- array(foo$gammasave, c(p, q, nsave))
  alphasave <- foo$alphasave
  sigmasave <- matrix(foo$sigmasave, nsave, 6)
  quansave <- matrix(foo$quansave, nsave, q*nquan)
  dens <- matrix(foo$f, ngrid, q)

#  aratesave <- matrix(foo$aratesave, nburn/100, 2*q + 8)
#  tunesave <- matrix(foo$tunesave, nburn/100, 2*q + 8)
  
  betatau <- matrix(0, nquan*q, p)
  for (j in 1:q){
    for (i in 1:nquan){
      tmp <- betasave[,j,]+gammasave[,j,]*as.numeric(quansave[,(j-1)*nquan+i])
      betatau[(j-1)*nquan+i,] <- apply(tmp, 1, mean)
    }
  }
  
  coef <- list(beta=apply(betasave, c(1,2), mean),
               gamma=apply(gammasave,c(1,2),mean),
               alpha=mean(alphasave),
               sigma=apply(sigmasave, 2, mean),
               quan=apply(quansave, 2, mean),
               betatau=betatau
               )

  z <- list(coef=coef,
            betasave=betasave,
            gammasave=gammasave,
            alphasave=alphasave,
            sigmasave=sigmasave,
            quansave=quansave,
            p=p,
            quan=quan,
            n=nsub,
            q=q,
            x=X,
            y=y,
            mcmc=mcmc,
            prior=prior,
            dens=dens,
            grid=grid
            )

  class(z) <- "HeterPTlmtri"

  return(z)
  
}

###############################################

summary.HeterPTlmtri <- function(obj){
  list(coef=obj$coef, n=obj$n, p=obj$p, q=obj$q, quan=obj$quan)
}

#################################################
plot.HeterPTlmtri <- function(obj, ask=FALSE){
  par(mfrow=c(obj$p, 2),ask=ask)
  for (j in 1:obj$q){
    for (i in 1:obj$p){
      title1 <- paste("Trace of beta" , i-1, sep=" ")
      title2 <- paste("Density of beta", i-1, sep=" ")
      plot(obj$betasave[i,j,], type='l', main=title1, xlab="MCMC scan", ylab=" ")
      plot(density(obj$betasave[i,j,]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
    }
  }

  for (j in 1:obj$q){
    for (i in 2:obj$p){
      title1 <- paste("Trace of gamma" , i-1, sep=" ")
      title2 <- paste("Density of gamma", i-1, sep=" ")
      plot(obj$gammasave[i,j,], type='l', main=title1, xlab="MCMC scan", ylab=" ")
      plot(density(obj$gammasave[i,j,]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
    }
  }

  for (i in 1:6){
    title1 <- paste("Trace of sigma2", i, sep=" ")
    title2 <- paste("Density of sigma2", i, sep= " ")
    plot(obj$sigmasave[,i], typ='l', main=title1, xlab="MCMC scan", ylab=" ")
    plot(density(obj$sigmasave[,i]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
  }

  title1 <- "Trace of alpha"
  title2 <- "Density of alpha"
  plot(obj$alphasave, typ='l', main=title1, xlab="MCMC scan", ylab=" ")
  plot(density(obj$alphasave), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')

  for (i in 1:q){
    title1 <- "Predictive Error Density"
    plot(obj$grid, obj$dens[,i], ylab="density", main=title1, type='l', lwd=2, xlab="values")
  }
}

######## boot summary ##############
bootsummary.HeterPTlmtri <- function(obj, truebetatau){
  mse <- matrix(0, 6, p)
  ## for (i in 1:6){
  ##   mse[i] <- mean((obj$coef$betatau[i,-1]-truebetatau[i, -1])^2)
  ## }

  mse <- (obj$coef$betatau - truebetatau)^2
  mseabs <- obj$coef$betatau - truebetatau
  mse <- as.vector(t(mse))
  mseabs <- as.vector(t(mseabs))
  
  return(list(mse=mse, mseabs=mseabs))
}


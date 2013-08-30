##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multivariate Bayesian Linear quantile regression with
##' Mixture of Polya tree priors , here d=2 
##' @param y 
##' @param X 
##' @param nsub 
##' @param mcmc 
##' @param prior 
##' @param quan 
##' @return class HeterPTlm
##' @author Minzhao Liu
HeterPTlmm <- function(y, X, nsub, mcmc, prior, quan){

  dyn.load('heterptlmm.so')
  
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
  sigmasave <- matrix(0, nsave, 3) # sigma1, sigma2, rho
  alphasave <- rep(0, nsave)
  quansave <- matrix(0, nsave, nquan*q)
  tunesave <- matrix(0, nburn, p*q+q+5)
  ratesave <- matrix(0, nburn/50, p*q+q+5)
  
  # GRID
  ngrid <- 200
  grid <- seq(-5,5, length=ngrid)
  f <- matrix(0, ngrid, q)
  
  # INITIAL
  beta <- solve(t(X)%*%X)%*%t(X)%*%y
  gamma <- matrix(0, p, q)
  gamma[1,] <- 1
  sigmavec <- c(1, 1, 0)
  Sigma <- diag(2)
  alpha <- 1
  b <- (y-X%*%beta)/(X%*%gamma)
  propv <- solve(t(X)%*%X)
  # WORKING
  whicho <- whichn <- rep(0, nsub)
  
  # DEBUG, TUNING
   # not used yet

  ####################################
  # MCMC USING FORTRAN

  foo <- .Fortran("heterptlmm",
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
                  quansave=as.double(quansave),
                  tunesave=as.double(tunesave),
                  ratesave=as.double(ratesave)                  
                  )

  ####################################

  betasave <- array(foo$betasave, dim=c(p, q, nsave))
  gammasave <- array(foo$gammasave, c(p, q, nsave))
  alphasave <- foo$alphasave
  sigmasave <- matrix(foo$sigmasave, nsave, 3)
  quansave <- matrix(foo$quansave, nsave, q*nquan)
  dens <- matrix(foo$f, ngrid, q)
  tunesave <- matrix(foo$tunesave, nburn, p*q+q+5)
  ratesave <- matrix(foo$ratesave, nburn/50, p*q+q+5)

  
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
            grid=grid,
            ratesave=ratesave,
            tunesave=tunesave
            )

  class(z) <- "HeterPTlmm"

  return(z)
  
}

###############################################

summary.HeterPTlmm <- function(obj){
  list(coef=obj$coef, n=obj$n, p=obj$p, q=obj$q, quan=obj$quan)
}

#################################################
plot.HeterPTlmm <- function(obj, ask=FALSE){
  par(mfrow=c(obj$p, 2),ask=ask)
  p <- obj$p
  q <- obj$q
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

  for (i in 1:3){
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


  apply(obj$betasave, c(1,2) , acf)
  apply(obj$gammasave[-1,,], c(1,2) , acf)
  acf(obj$sigmasave)
  acf(obj$alphasave)

  for (i in 1:(p*q+q+5)){
    plot(ts(obj$tunesave[,i]),cex=0.5)
    plot(ts(obj$ratesave[,i]),cex=0.5)
  }

}

######## boot summary ##############
bootsummary.HeterPTlmm <- function(obj, truebetatau){

  mse <- c(t(obj$coef$betatau - truebetatau))

  return(list(mse = mse))
}


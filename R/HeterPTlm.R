####################################################################
##' Time-stamp: <liuminzhao 06/24/2013 14:50:44>
##'
##' 2012/03/29 wrap heterptlm.f,
####################################################################
##' .. content for \description{} (no empty lines) ..
##' Bayesian Quantile Regression with Polya Tree
##' for univariate outcome. Unify three methods:
##' mean: Bayesian posterior mean
##' median: Bayesian posterior median
##' ss: Spike-stab prior
##' .. content for \details{} ..
##' @title
##' @param y
##' @param x
##' @param mcmc
##' @param prior
##' @param quan
##' @param method
##' @return
##' @author Minzhao Liu, Mike Daniels
HeterPTlm <- function(y, x, mcmc, prior, quan = 0.5, method = "ss"){

  ## LOAD DYN
  if (method == 'mean' | method == 'median'){
    if (!is.loaded('heterptlm')) {
      if (file.exists('heterptlm.so')) {
        dyn.load('heterptlm.so')
      } else {
        if (file.exists('~/Documents/bqrpt/code/heterptlm.so')) {
          dyn.load("~/Documents/bqrpt/code/heterptlm.so")
        } else stop('no shared library found.')
      }
    }
  }
  if (method == 'ss') {
    if (!is.loaded('heterptlmss')) {
      if (file.exists('heterptlm-ss.so')) {
        dyn.load('heterptlm-ss.so')
      } else {
        if (file.exists('~/Documents/bqrpt/code/heterptlm-ss.so')) {
          dyn.load("~/Documents/bqrpt/code/heterptlm-ss.so")
        } else stop('no shared library found.')
      }
    }
  }

  ## DATA
  nrec <- length(y)
  p <- dim(x)[2]

  ## PT
  if (is.null(prior$maxm)) {
    maxm <- floor(log(nrec)/log(2)) }  else maxm <- prior$maxm
  if (is.null(prior$mdzero)) {mdzero <- 0}  else  mdzero <- prior$mdzero

  ## PRIOR
  if (is.null(prior$betapm)){
    betapm <- solve(t(x)%*%x)%*%t(x)%*%y
    res <- y - x%*%betapm
    betapv <- solve(t(x)%*%x)*sum(res^2)/(nrec - p)
    gammapm <- rep(0, p)
    gammapv <- diag(p)*100
    sigmap <- c(2, 1/2)
    a0b0 <- c(1, 1)
  } else {
    betapm <- prior$betapm
    betapv <- prior$betapv
    gammapm <- prior$gammapm
    gammapv <- prior$gammapv
    sigmap <- prior$sigmap
    a0b0 <- prior$a0b0
  }

  ## MCMC
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nburn <- mcmc$nburn
  ndisp <- mcmc$ndisp
  arate <- mcmc$arate
  mcmc <- c(nburn, nskip, nsave, ndisp)

  ## QUAN
  quan <- quan
  nquan <- length(quan)

  ## SAVE
  betasave <- gammasave <- matrix(0, nsave, p)
  sigmasave <- alphasave <- rep(0, nsave)
  quansave <- matrix(0, nsave, nquan)

  ## grid
  ngrid <- 200
  f <- rep(0, ngrid)

  ## INITIAL
  beta <- as.vector(solve(t(x)%*%x)%*%t(x)%*%y)
  gamma <- c(1, rep(0, p-1))
  sigma2 <- 1
  alpha <- 1
  v <- as.vector((y-x%*%beta)/(x%*%gamma))

  ## new grid

  left <- min(v) - 0.5*sd(v)
  right <- max(v) + 0.5*sd(v)
  grid <- seq(left, right, length=ngrid)

  ## WORKING
  whicho <- whichn <- rep(0, nrec)

  ## DEBUG
  ratesave <- matrix(0, nburn/50, 2*p+2)
  tunesave <- matrix(0, nburn, 2*p+2)
#  hetersave <- rep(0, nburn)
  hetersave <- matrix(0, nburn, p)
  propv <- solve(t(x)%*%x)
  ###############################################

  if (!method == 'ss') {
    foo <- .Fortran("heterptlm",
                    maxm=as.integer(maxm),
                    mdzero=as.integer(mdzero),
                    nrec=as.integer(nrec),
                    p = as.integer(p),
                    x = as.double(x),
                    y = as.double(y),
                    betapm= as.double(betapm),
                    betapv = as.double(betapv),
                    sigmap = as.double(sigmap),
                    betasave= as.double(betasave),
                    gammasave= as.double(gammasave),
                    gammapm = as.double(gammapm),
                    gammapv=as.double(gammapv),
                    alpha = as.double(alpha),
                    beta = as.double(beta),
                    gamma = as.double(gamma),
                    nsave = as.integer(nsave),
                    sigma2 = as.double(sigma2),
                    v = as.double(v),
                    a0b0=as.double(a0b0),
                    mcmc = as.integer(mcmc),
                    whicho=as.integer(whicho),
                    whichn=as.integer(whichn),
                    sigmasave=as.double(sigmasave),
                    alphasave=as.double(alphasave),
                    f=as.double(f),
                    ngrid=as.integer(ngrid),
                    grid=as.double(grid),
                    quan=as.double(quan),
                    quansave=as.double(quansave),
                    nquan=as.integer(nquan),
                    ratesave=as.double(ratesave),
                    tunesave=as.double(tunesave),
                    hetersave=as.double(hetersave),
                    propv=as.double(propv),
                    arate=as.double(arate)
                    )
  } else {
      foo <- .Fortran("heterptlmss",
                  maxm=as.integer(maxm),
                  mdzero=as.integer(mdzero),
                  nrec=as.integer(nrec),
                  p = as.integer(p),
                  x = as.double(x),
                  y = as.double(y),
                  betapm= as.double(betapm),
                  betapv = as.double(betapv),
                  sigmap = as.double(sigmap),
                  betasave= as.double(betasave),
                  gammasave= as.double(gammasave),
                  gammapm = as.double(gammapm),
                  gammapv=as.double(gammapv),
                  alpha = as.double(alpha),
                  beta = as.double(beta),
                  gamma = as.double(gamma),
                  nsave = as.integer(nsave),
                  sigma2 = as.double(sigma2),
                  v = as.double(v),
                  a0b0=as.double(a0b0),
                  mcmc = as.integer(mcmc),
                  whicho=as.integer(whicho),
                  whichn=as.integer(whichn),
                  sigmasave=as.double(sigmasave),
                  alphasave=as.double(alphasave),
                  f=as.double(f),
                  ngrid=as.integer(ngrid),
                  grid=as.double(grid),
                  quan=as.double(quan),
                  quansave=as.double(quansave),
                  nquan=as.integer(nquan),
                  ratesave=as.double(ratesave),
                  tunesave=as.double(tunesave),
                  hetersave=as.double(hetersave),
                  propv=as.double(propv),
                  arate=as.double(arate)
                  )
    }
    #################################################
  betasave <- matrix(foo$betasave, nsave, p)
  gammasave <- matrix(foo$gammasave, nsave, p)
  quansave <- matrix(foo$quansave, nsave,nquan)

  if (method == 'median') {
    coef.beta <- apply(betasave, 2, median)
    coef.gamma <- apply(gammasave, 2, median)
    coef.quan <- apply(quansave,2,median)
    coef.betatau <- matrix(0, nquan, p)
    rownames(coef.betatau) <- quan
    for (i in 1:nquan) {
      tmp <- betasave + gammasave*as.numeric(quansave[,i])
      coef.betatau[i, ] <- apply(tmp, 2, median)
    }
  } else if (method == 'mean' | method == 'ss') {
    coef.beta <- apply(betasave, 2, mean)
    coef.gamma <- apply(gammasave, 2, mean)
    coef.quan <- apply(quansave,2,mean)
    coef.betatau <- matrix(0, nquan, p)
    rownames(coef.betatau) <- quan
    for (i in 1:nquan) {
      tmp <- betasave + gammasave*as.numeric(quansave[,i])
      coef.betatau[i, ] <- apply(tmp, 2, mean)
    }
  }

  coef <- list(beta=coef.beta, gamma=coef.gamma, quan=coef.quan,betatau=coef.betatau)

  ## DEBUG
  ratesave <- matrix(foo$ratesave, nburn/50, 2*p+2)
  tunesave <- matrix(foo$tunesave, nburn, 2*p+2)
  hetersave <- matrix(foo$hetersave, nburn, p)

  z <- list(coef=coef,
            betasave=betasave,
            gammasave=gammasave,
            quansave=quansave,
            sigmasave=foo$sigmasave,
            dens=foo$f,
            grid=grid,
            mcmc=mcmc,
            prior=prior,
            n=nrec,
            p=p,
            quan=quan,
            y=y,
            X=x,
            ratesave=ratesave,
            tunesave=tunesave,
            hetersave=hetersave,
            alphasave=foo$alphasave,
            method = method
            )

  class(z) <- "HeterPTlm"

  return(z)
}


############################################################
plot.HeterPTlm <- function(obj, ask=FALSE){
  par(mfrow=c(obj$p, 2),ask=ask)
  for (i in 1:obj$p){
    title1 <- paste("Trace of beta" , i-1, sep=" ")
    title2 <- paste("Density of beta", i-1, sep=" ")
    plot(obj$betasave[,i], type='l', main=title1, xlab="MCMC scan", ylab=" ")
    plot(density(obj$betasave[,i]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
  }

  for (i in 2:obj$p){
    title1 <- paste("Trace of gamma" , i-1, sep=" ")
    title2 <- paste("Density of gamma", i-1, sep=" ")
    plot(obj$gammasave[,i], type='l', main=title1, xlab="MCMC scan", ylab=" ")
    plot(density(obj$gammasave[,i]), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')
  }

  title1 <- "Trace of sigma2"
  title2 <- "Density of sigma2"
  plot(obj$sigmasave, typ='l', main=title1, xlab="MCMC scan", ylab=" ")
  plot(density(obj$sigmasave), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')

  title1 <- "Trace of alpha"
  title2 <- "Density of alpha"
  plot(obj$alphasave, typ='l', main=title1, xlab="MCMC scan", ylab=" ")
  plot(density(obj$alphasave), lwd=1.2, main=title2, xlab="values", ylab="density", col='red')

  title1 <- "Predictive Error Density"
  plot(obj$grid, obj$dens, ylab="density", main=title1, type='l', lwd=2, xlab="values")
}

############################################################
summary.HeterPTlm <- function(obj){
#  list(coef=obj$coef,n=obj$n, p=obj$p, quan=obj$quan, mcmc=obj$mcmc, prior=obj$prior)
  list(coef=obj$coef,n=obj$n, p=obj$p, quan=obj$quan)
}
###########################################################
bootsummary.HeterPTlm <- function(obj, truebetatau){
  nquan <- dim(truebetatau)[1]
  p <- obj$p
  quan <- obj$quan
  runs <- obj$mcmc$nsave

  betatau1 <- obj$betasave + obj$gammasave*as.numeric(obj$quansave[,1])
  betatau2 <- obj$betasave + obj$gammasave*as.numeric(obj$quansave[,2])

  betatau.coef1 <- apply(betatau1, 2, median)
  betatau.coef2 <- apply(betatau2, 2, median)

  lbd1 <- apply(betatau1, 2, function(x) quantile(x, 0.025))
  lbd2 <- apply(betatau2, 2, function(x) quantile(x, 0.025))

  ubd1 <- apply(betatau1, 2, function(x) quantile(x, 0.975))
  ubd2 <- apply(betatau2, 2, function(x) quantile(x, 0.975))

  len1 <- ubd1-lbd1
  len2 <- ubd2-lbd2

  mse1 <- betatau.coef1-truebetatau[1,]
  mse2 <- betatau.coef2-truebetatau[2,]

  for (j in 2:p)   {
    cover1 <- prod((truebetatau[1,j]>lbd1[j] && truebetatau[1,j]<ubd1[j]))
    cover2 <- prod((truebetatau[2,j]>lbd2[j] && truebetatau[2,j]<ubd2[j]))
  }

  return(list(beta1=betatau.coef1, lbd1=lbd1, ubd1=ubd1, length1=len1, mse1=mse1, cover1=cover1,
              beta2=betatau.coef2, lbd2=lbd2, ubd2=ubd2, length2=len2, mse2=mse2, cover2=cover2))


}
#############################################################
Diagnose.HeterPTlm <- function(obj, ask=FALSE){
  par(mfrow=c(4,4))
  for (i in 1:8){
    plot(obj$tunesave[,i],cex=0.5)
    plot(obj$ratesave[,i],cex=0.5)
  }

  par(mfrow=c(2,2))
  plot(obj$hetersave[,2],type='l')
  plot(obj$hetersave[,3],type='l')
  plot(obj$sigmasave,type='l')
  plot(obj$alphasave,type='l')

  # auto correlation
  par(mfrow=c(2,3))
  for (i in 1:p){
    acf(obj$betasave[,i])
  }
  for (i in 2:p){
    acf(obj$gammasave[,i])
  }
  acf(obj$sigmasave)

}

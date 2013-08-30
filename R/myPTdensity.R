rm(list=ls())
n <- 100
y <- rnorm(n)
maxm <- floor(log(n)/log(2))
mdzero <- 0
ngrid <- 200
f <- rep(0, ngrid)
left <- min(y) - 0.5*sd(y)
right <- max(y) + 0.5*sd(y)
grid <- seq(left, right, length=ngrid)
arate <- 0.25
mu <- 0
sigma <- 1
alpha <- 1
whicho <- whichn <- rep(0, n)
logprioro <- 0

dyn.load("myptdensity.so")

GetLogDensity <- function(grid){
  tmp <- .Fortran('gridupptprior',
                  theta=as.double(grid),
                  maxm=as.integer(maxm),
                  mdzero=as.integer(0),
                  nsubject=as.integer(n),
                  alpha=as.double(alpha),
                  mu=as.double(mu),
                  sigma=as.double(sigma),
                  b=as.double(y),
                  whicho=as.integer(whicho),
                  whichn=as.integer(whichn),
                  ll=as.double(logprioro)
                  )
  tmp$ll
}

y <- rMN3(n)
for (i in 1:ngrid) 
  f[i] <- exp(GetLogDensity(grid[i]))

plot(grid,f, 'l')

myPTdensity <- function(y, mcmc,   maxm = floor(log(length(y))/log(2))){
  dyn.load("myptdensity.so")
  n <- length(y)
  maxm <- maxm
  mdzero <- 0
  ngrid <- 200
  f <- rep(0, ngrid)
  left <- min(y) - 0.5*sd(y)
  right <- max(y) + 0.5*sd(y)
  grid <- seq(left, right, length=ngrid)
  arate <- 0.25
  mu <- 0
  sigma <- 1
  alpha <- 1
  nburn <- mcmc[1]
  nsave <- mcmc[2]
  nskip <- mcmc[3]
  sigmasave <- alphasave <- musave <- rep(0, nsave)
  fsave <- matrix(0, nsave, ngrid)
  ratesave <- tunesave <- matrix(0, nburn, 3)
  
  foo <- .Fortran("myptdensity1",
                  maxm  = as.integer(maxm),
                  n     = as.integer(n),
                  y     = as.double(y),
                  mu    = as.double(mu),
                  sigma = as.double(sigma),
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

  ratesave <- matrix(foo$ratesave, nburn, 3)
  tunesave <- matrix(foo$tunesave, nburn, 3)
  dens <- foo$f
  plot(grid, dens, 'l')
  return(list(foo=foo, ratesave=ratesave, tunesave=tunesave))
}

y <- rnorm(100)
rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}
y <- rMN3(100) + 5
mcmc <- c(10000, 1000, 10)
mcmc <- c(20000, 10000, 20)
a <- myPTdensity(y, mcmc)
a <- myPTdensity(y, mcmc, maxm = 4)
a2 <- myPTdensity(y, mcmc)
a3 <- myPTdensity(y, mcmc, maxm = 10)

plot(a$ratesave[,1])
plot(a$ratesave[,2])
plot(a$ratesave[,3])

plot(a$tunesave[,1])
plot(a$tunesave[,2])
plot(a$tunesave[,3])

png('u3-mixing.png')
par(mfrow = c(4, 2))
plot(ts(a$foo$musave))
plot(ts(a$foo$sigmasave))
plot(ts(a$foo$alphasave))
acf(a$foo$musave)
acf(a$foo$sigmasave)
acf(a$foo$alphasave)
plot(a$foo$grid, a$foo$f, 'l')
dev.off()

plot(a$foo$grid, a$foo$f, 'l', col='black')
lines(a2$foo$grid, a2$foo$f, 'l', col='red', lty=2)
lines(a3$foo$grid, a3$foo$f, 'l', col='blue', lty=3)
legend("topright",c("M=4", "M=6", "M =10"), col=c('black', 'red', 'blue'), lty=1:3)

####################

myPTdensity2 <- function(y, mcmc, maxm = floor(log(length(y))/log(2))){
  dyn.load("myptdensity2.so")
  n <- length(y)
  maxm <- maxm
  mdzero <- 0
  ngrid <- 200
  f <- rep(0, ngrid)
  left <- min(y) - 0.5*sd(y)
  right <- max(y) + 0.5*sd(y)
  grid <- seq(left, right, length=ngrid)
  arate <- 0.25
  mu <- 0
  sigma <- 1
  alpha <- 1
  nburn <- mcmc[1]
  nsave <- mcmc[2]
  nskip <- mcmc[3]
  sigmasave <- alphasave <- musave <- rep(0, nsave)
  fsave <- matrix(0, nsave, ngrid)
  ratesave <- tunesave <- matrix(0, nburn, 3)
  
  foo <- .Fortran("myptdensity2",
                  maxm  = as.integer(maxm),
                  n     = as.integer(n),
                  y     = as.double(y),
                  mu    = as.double(mu),
                  sigma = as.double(sigma),
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

  ratesave <- matrix(foo$ratesave, nburn, 3)
  tunesave <- matrix(foo$tunesave, nburn, 3)
  dens <- foo$f
  plot(grid, dens, 'l')
  return(list(foo=foo, ratesave=ratesave, tunesave=tunesave))
}

y <- rnorm(100)
rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}
y <- rMN3(500)
mcmc <- c(10000, 1000, 10)
b <- myPTdensity2(y, mcmc)

plot(b$ratesave[,1])
plot(b$ratesave[,2])
plot(b$ratesave[,3])

plot(b$tunesave[,1])
plot(b$tunesave[,2])
plot(b$tunesave[,3])

plot(b$foo$grid, b$foo$f, 'l')
plot(b$foo$grid, b$foo$f, 'l', col='red', lty=2)
lines(a$foo$grid, a$foo$f, 'l', col='blue', lty=1)

# blue is more smooth 

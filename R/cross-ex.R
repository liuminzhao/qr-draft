rm(list = ls())
source('BQRboot.R')
library(bqrpt)
set.seed(1)
n <- 200
tuneinit <- c(0.3, 0.3, 0.3, 1, 0.3, 0.3, 0.04, 0.1)
mcmc <- list(nburn=30000, nskip=5, nsave=30000, ndisp=10000, arate=0.2, tuneinit = tuneinit)
b1 <- 1
quan <- c(0.25, 0.5, 0.75, 0.85, 0.9, 0.91)
x1 <- runif(n, max = 4)
e1 <- rnorm(n)

y1 <- 1 + x1*b1 + e1

X <- cbind(1,x1)
g1 <- 0.2

x1 <- runif(n, max = 4)
e1 <- rt(n, df = 3)

X <- cbind(1,x1)

y1 <- 1 + x1*b1 + (1 + x1*g1)*e1


## bqrpt
modpt <- HeterPTlmMH(y1, X, mcmc = mcmc, quan = quan, prior = list(maxm = 6))
yfit <- newdata %*% t(coef(modpt)$betatauMedian)
matplot(yfit, type = 'l')

## bqr
modbqr1 <- BayesQReg(y1, X, 0.25)
modbqr2 <- BayesQReg(y1, X, 0.5)
modbqr3 <- BayesQReg(y1, X, 0.75)
modbqr4 <- BayesQReg(y1, X, 0.85)
modbqr5 <- BayesQReg(y1, X, 0.9)
modbqr6 <- BayesQReg(y1, X, 0.95)
modbqr6 <- BayesQReg(y1, X, 0.91)
coefbqr1 <- BQRCoef(modbqr1)
coefbqr2 <- BQRCoef(modbqr2)
coefbqr3 <- BQRCoef(modbqr3)
coefbqr4 <- BQRCoef(modbqr4)
coefbqr5 <- BQRCoef(modbqr5)
coefbqr6 <- BQRCoef(modbqr6)

## plot

plot(x1*10, y1)
newdata <- cbind(1, x, xsq)


yfit1 <- newdata %*% coefbqr1
yfit2 <- newdata %*% coefbqr2
yfit3 <- newdata %*% coefbqr3
yfit4 <- newdata %*% coefbqr4
yfit5 <- newdata %*% coefbqr5
yfit6 <- newdata %*% coefbqr6
lines(x*10, yfit1, col = 1)
lines(x*10, yfit2, col = 2)
lines(x*10, yfit3, col = 3)
lines(x*10, yfit4, col = 4)
lines(x*10, yfit5, col = 5)
lines(x*10, yfit6, col = 6)



abline(coefbqr1, col = 1)
abline(coefbqr2, col = 2)
abline(coefbqr3, col = 3)
abline(coefbqr4, col = 4)
abline(coefbqr5, col = 5)
abline(coefbqr6, col = 6)


### cars
x1 <- cars$speed/10
y1 <- cars$dist
x2 <- x1^2
X <- cbind(1, x1, x2)

library(quantreg)
mod <- rq(y1 ~ x1 + x2, tau = c(0.25, 0.5, 0.75, 0.85, 0.9, 0.91))

x <- seq(2, 25, by=0.1)/10
xsq <- x^2
new <- data.frame(x1=x, x2=xsq)
pred <- predict(mod, newdata=new)
matplot(x, pred, lty = 1, type = 'l')

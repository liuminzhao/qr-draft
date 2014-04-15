rm(list = ls())
source('BQRboot.R')
set.seed(1)
n <- 200

## sim - m1

x1 <- runif(n, max = 4)
e1 <- rnorm(n)

y1 <- 1 + x1*b1 + e1

X <- cbind(1,x1)


## sim m2
  x1 <- runif(n, max = 4)
  e1 <- rt(n, df = 3)

  y1 <- 1 + x1*b1 + e1

  X <- cbind(1,x1)

## m3
rMN3 <- function(n){
  posneg <- rbinom(n,1,0.5)
  (1-posneg)*rnorm(n, -2, 1) + posneg*rnorm(n, 2,1)
}

x1 <- runif(n, max = 4)
  e1 <- rMN3(n)
  y1 <- 1 + x1*b1 + e1
  X <- cbind(1,x1)


## m4

rMN <-function(n){
  posneg<-rbinom(n,1,0.8)
  posneg*rnorm(n) + (1-posneg)*rnorm(n, 3, sqrt(3))
}

  x1 <- runif(n, max = 4)
  e1 <- rMN(n)

  X <- cbind(1,x1)

  y1 <- 1 + x1*b1 + e1


modbqr1 <- BayesQReg(y1, X, 0.1)
modbqr3 <- BayesQReg(y1, X, 0.3)
modbqr5 <- BayesQReg(y1, X, 0.5)
modbqr7 <- BayesQReg(y1, X, 0.7)
modbqr9 <- BayesQReg(y1, X, 0.9)
modbqr75 <- BayesQReg(y1, X, 0.75)
modbqr91 <- BayesQReg(y1, X, 0.91)
modbqr95 <- BayesQReg(y1, X, 0.95)

coefmat <- rbind(BQRCoef(modbqr1), BQRCoef(modbqr3), BQRCoef(modbqr5), BQRCoef(modbqr7), BQRCoef(modbqr75), BQRCoef(modbqr9), BQRCoef(modbqr91), BQRCoef(modbqr95))
rownames(coefmat) <- c('10', '30', '50', '70', '75','90', '91','95')
colnames(coefmat) <- c('int', 'x')

## plot race

png('m4-cross.png')
plot(c(0, 4), c(-2, 10), type = 'n', xlab = 'x', ylab  ='y')
for (i in 1:8) abline(coefmat[i, ], col = 9-i)
legend('bottomright', c('10', '30', '50', '70', '75', '90','91', '95'), col = 8:1, lty = 1)
dev.off()

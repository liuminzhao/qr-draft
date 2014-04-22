##' 2014-04-20 to summarise simulation result into latex table

library(xtable)

boot <- 100

result1 <- read.table('sim-m1-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1+qnorm(0.9),1)
t1 <- rep(c(truebetatau5, truebetatau9), 4)

g1 <- 0.2
result1h <- read.table('sim-m1h-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1,1) + c(1 , g1)*(qnorm(0.9))
t1h <- rep(c(truebetatau5, truebetatau9), 4)

result2 <- read.table('sim-m2-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1+qt(0.9, df = 3),1)
t2 <- rep(c(truebetatau5, truebetatau9), 4)

result2h <- read.table('sim-m2h-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1,1) + c(1 , g1)*(qt(0.9, df = 3))
t2h <- rep(c(truebetatau5, truebetatau9), 4)

result3 <- read.table('sim-m3-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1+2 + qnorm(0.8),1)
t3 <- rep(c(truebetatau5, truebetatau9), 4)

result3h <- read.table('sim-m3h-result-0420.txt')
truebetatau5 <- c(1,1)
truebetatau9 <- c(1,1) + c(1 , g1)*(2 + qnorm(0.8))
t3h <- rep(c(truebetatau5, truebetatau9), 4)

result4 <- read.table('sim-m4-result-0420.txt')
pMN.5 <- function(x){
  0.8*pnorm(x)+0.2*pnorm(x,3,sqrt(3))-0.5
}
pMN.9 <- function(x){
  0.8*pnorm(x)+0.2*pnorm(x,3,sqrt(3))-0.9
}
q5 <- uniroot(pMN.5, c(-5,5))$root
q9 <- uniroot(pMN.9, c(-5,5))$root
truebetatau5 <- c(1,1) + c(1, 0)*q5
truebetatau9 <- c(1,1) + c(1, 0)*q9
t4 <- rep(c(truebetatau5, truebetatau9), 4)

result4h <- read.table('sim-m4h-result-0420.txt')
truebetatau5 <- c(1,1) + c(1, g1)*q5
truebetatau9 <- c(1,1) + c(1, g1)*q9
t4h <- rep(c(truebetatau5, truebetatau9), 4)

alpha <- 0
result5 <- read.table('sim-m5-result-0420.txt')
quan1 <- function(y, x, tau){
  return(tau - .5*pnorm(y, 2+x, 1 + alpha*x) - .5*pnorm(y, -2-x, 1 + alpha*x))
}
SolveQuan1 <- function(x, tau){
  uniroot(quan1, c(-30, 30), x = x, tau = tau)$root
}
xsim <- seq(-1, 1, len = 100)
y11 <- sapply(xsim, function(x) SolveQuan1(x, 0.1))
y13 <- sapply(xsim, function(x) SolveQuan1(x, 0.3))
y15 <- sapply(xsim, function(x) SolveQuan1(x, 0.5))
y17 <- sapply(xsim, function(x) SolveQuan1(x, 0.7))
y19 <- sapply(xsim, function(x) SolveQuan1(x, 0.9))
q11 <- lm(y11~xsim)$coef
q13 <- lm(y13~xsim)$coef
q15 <- lm(y15~xsim)$coef
q17 <- lm(y17~xsim)$coef
q19 <- lm(y19~xsim)$coef
truebetatau5 <- q15
truebetatau9 <- q19
t5 <- rep(c(truebetatau5, truebetatau9), 4)

## MSE

getmsemcse <- function(result, true){
    mse <- rep(0, 16)
    for (i in 1:16){
        mse[i] <- mean((result[,i] - true[i])^2)*100
    }
    mse <- matrix(mse, 4, 4)
    colnames(mse) <- c('RQ', 'BQR', 'PT', 'PTSS')

    ## MCSE
    mcse <- rep(0, 16)
    for (i in 1:16){
        mcse[i] <- sd((result[,i] - true[i])^2)/sqrt(boot) * 100
    }
    mcse <- matrix(mcse, 4, 4)
    colnames(mcse) <- c('RQ', 'BQR', 'PT', 'PTSS')

    ## combine mse and MCSE
    msemcse <- matrix(0, 4, 8)
    for (i in 1:4){
        msemcse[, i*2 - 1] <- mse[, i]
        msemcse[, i*2] <- mcse[, i]
    }

    return(msemcse)
}

msemcse1 <- getmsemcse(result1, t1)
msemcse1h <- getmsemcse(result1, t1h)
msemcse2 <- getmsemcse(result2, t2)
msemcse2h <- getmsemcse(result2h, t2h)
msemcse3 <- getmsemcse(result3, t3)
msemcse3h <- getmsemcse(result3h, t3h)
msemcse4 <- getmsemcse(result4, t4)
msemcse4h <- getmsemcse(result4h, t4h)
msemcse5 <- getmsemcse(result5, t5)
void <- matrix(0, 4, 8)

table1 <- cbind(msemcse1, msemcse1h)
table2 <- cbind(msemcse2, msemcse2h)
table3 <- cbind(msemcse3, msemcse3h)
table4 <- cbind(msemcse4, msemcse4h)
table5 <- cbind(msemcse5, void)

whole <- rbind(table1, table2, table3, table4, table5)

new <- matrix(0, nrow(whole), ncol(whole)/2)

for (i in 1:ncol(new)){
    new[, i] <- paste0(round(whole[, 2*i - 1], 2), "(", round(whole[, 2*i], 2), ")")
}

colnames(new) <- rep(c("RQ", "FBQR", "PT", "PTSS"), 2)
rownames(new) <- rep(c('\beta0', '\beta1'), 10)

xtable(new)

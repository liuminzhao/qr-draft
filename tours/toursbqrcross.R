modbqr1 <- BayesQReg(y, X, 0.1)
modbqr15 <- BayesQReg(y, X, 0.15)
modbqr3 <- BayesQReg(y, X, 0.3)
modbqr5 <- BayesQReg(y, X, 0.5)
modbqr7 <- BayesQReg(y, X, 0.7)
modbqr85 <- BayesQReg(y, X, 0.85)
modbqr9 <- BayesQReg(y, X, 0.9)
modbqr91 <- BayesQReg(y, X, 0.91)
modbqr95 <- BayesQReg(y, X, 0.95)

## new 10, 15, 30, 50, 70, 85, 90, 91 from BQR for tours

coefmat <- rbind(BQRCoef(modbqr1), BQRCoef(modbqr15), BQRCoef(modbqr3), BQRCoef(modbqr5), BQRCoef(modbqr7), BQRCoef(modbqr85), BQRCoef(modbqr9), BQRCoef(modbqr91), BQRCoef(modbqr95))
rownames(coefmat) <- c('10', '15', '30', '50', '70','85', '90','91', '95')
colnames(coefmat) <- c('int', 'age', 'race')

## plot race
png('ToursBQRRace.png')
plot(c(0, 1), c(0, 20), type = 'n', xlab = 'race', ylab  ='weight loss')
for (i in 1:9) abline(coefmat[i, c(1, 3)], col = 10-i)
legend('bottomright', c('10', '15', '30', '50', '70','85', '90','91', '95'), col = 9:1, lty = 1)
dev.off()


png('ToursBQRAge.png')
plot(c(0, 5), c(0, 20), type = 'n', xlab = 'Age', ylab  ='weight loss')
for (i in 1:9) abline(coefmat[i, c(1, 2)], col = 10-i)
legend('topright', c('10', '15', '30', '50', '70','85', '90','91', '95'), col = 9:1, lty = 1)
dev.off()

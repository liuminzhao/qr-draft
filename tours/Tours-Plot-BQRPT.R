rm(list=ls())
set.seed(1)

TOURS <- read.csv('~/Documents/qr-draft/tours/tours.csv')

TOURS <- subset(TOURS, RACE==1 | RACE==3)

weight1 <- TOURS$wtkg1
weight2 <- TOURS$wtkg2
weight3 <- TOURS$wtkg3
age <- TOURS$AGE
trt <- TOURS$TREATMENT
age_center <- (age-50)/25
race3 <- as.numeric(TOURS$RACE == 3)

y <- weight1 - weight2
n <- length(y)
X <- matrix(0, n, 3)
X[,1] <- 1
X[,2] <- age_center
X[,3] <- race3

dat <- data.frame(loss = y, weight2, weight3, trt, age_center, age=TOURS$AGE, race = factor(TOURS$RACE), base = weight1)

###############
## PLOT
###############
library(ggplot2)
library(gridExtra)

box1 <- ggplot(data = dat, aes(x = age, y = loss)) + geom_point() + ylab('Weight Loss (Kg)') + xlab('Age') + ggtitle("Weight Loss vs Age")

box2 <- ggplot(data = dat, aes(x = race, y = loss)) + geom_boxplot() + scale_x_discrete(labels=c("Black", "White")) + ylab('Weight Loss (Kg)') + xlab('Race') + ggtitle("Weight Loss vs Race")

box3 <- ggplot(data = dat, aes(x = age, y = loss, color = race)) + geom_point() + ylab('Weight Loss (Kg)') + xlab('Age') + ggtitle("Weight Loss vs Age") + ggtitle('Weight Loss vs Age and Race')

pdf('../image/weight-age-race.pdf', width = 20, height = 7)
sds <- grid.arrange(box1, box2, box3, nrow = 1, ncol = 3)
dev.off()

library(pcalg)

##################################################
## Generiere Modell
##################################################
set.seed(123)
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

## only 2 Variables
b0 <- 0
b1 <- 1
n <- 1000
x1 <- sample(c(0,1),n,replace=TRUE)
p <- expit(b0 + b1*x1)
x2 <- rep(NA,length(x1))
for (i in 1:n) {
  x2[i] <- sample(c(0,1),1,prob=c(1-p[i],p[i]))
}

res1 <- fisher.test(x1,x2)$p.value
res2 <- gSquareBin(1,2,NULL,cbind(x1,x2))
ok1 <- round(res1,15) == round(res2,15)

##################################################
## chain of 3 variables: x1 -> x2 -> x3
##################################################
b0 <- 0
b1 <- 1
b2 <- 1
n <- 10000
x1 <- sample(c(0,1),n,replace=TRUE)

p2 <- expit(b0 + b1*x1)
x2 <- rep(NA,length(x1))
for (i in 1:n) {
  x2[i] <- sample(c(0,1),1,prob=c(1-p2[i],p2[i]))
}
p3 <- expit(b0 + b2*x2)
x3 <- rep(NA,length(x2))
for (i in 1:n) {
  x3[i] <- sample(c(0,1),1,prob=c(1-p3[i],p3[i]))
}

dat <- cbind(x1,x2,x3)

## should all be TRUE
ok2 <- all(c(gSquareBin(3,1,2,dat)>0.05,
             gSquareBin(1,3,2,dat)>0.05,
             gSquareBin(1,2,3,dat)<0.05,
             gSquareBin(3,2,1,dat)<0.05))

##################################################
## collider: x1 -> x3 <- x2
##################################################
b0 <- 0
b1 <- 1.3
b2 <- 1.7
n <- 10000
x1 <- sample(c(0,1),n,replace=TRUE)
x2 <- sample(c(0,1),n,replace=TRUE)

p3 <- expit(b0 + b1*x1 + b2*x2)
x3 <- rep(NA,n)
for (i in 1:n) {
  x3[i] <- sample(c(0,1),1,prob=c(1-p3[i],p3[i]))
}

dat <- cbind(x1,x2,x3)

## should all be TRUE
ok3 <- all(gSquareBin(1,2,3,dat)<0.05,
           gSquareBin(1,2,NULL,dat)>0.05,
           gSquareBin(2,1,3,dat)<0.05,
           gSquareBin(1,3,2,dat)<0.05)

##################################################
## add noise variables to collider model
##################################################
x4 <- sample(c(0,1),n,replace=TRUE)
x5 <- sample(c(0,1),n,replace=TRUE)
x6 <- sample(c(0,1),n,replace=TRUE)

dat <- cbind(x1,x2,x3,x4,x5,x6)

ok4 <- all(gSquareBin(1,2,c(3,4,5,6),dat)<0.05,
           gSquareBin(4,2,c(3,5,6),dat)>0.05)

##################################################
## rectangle model
##################################################
b0 <- 0
b1 <- 1
b2 <- 1
b3 <- 0.5
b4 <- 0.5
n <- 10000
x1 <- sample(c(0,1),n,replace=TRUE)

p2 <- expit(b0 + b1*x1)
p3 <- expit(b0 + b2*x1)
x2 <- rep(NA,n)
x3 <- rep(NA,n)
for (i in 1:n) {
  x2[i] <- sample(c(0,1),1,prob=c(1-p2[i],p2[i]))
  x3[i] <- sample(c(0,1),1,prob=c(1-p3[i],p3[i]))
}

p4 <- expit(b0 + b3*x2 + b4*x3)
x4 <- rep(NA,n)
for (i in 1:n) {
  x4[i] <- sample(c(0,1),1,prob=c(1-p4[i],p4[i]))
}

dat <- cbind(x1,x2,x3,x4)

ok5 <- all(gSquareBin(1,4,c(2,3),dat)>0.05,
           gSquareBin(1,4,2,dat)<0.05,
           gSquareBin(1,4,3,dat)<0.05,
           gSquareBin(1,4,NULL,dat)<0.05,
           gSquareBin(2,3,4,dat)<0.05,
           gSquareBin(2,3,1,dat)>0.05)

ok <- all(ok1,ok2,ok3,ok4,ok5)

if (!ok) stop("Test gSquareBin wrong: Some dependence was not estimated correctly!")

library(pcalg)

set.seed(100)

wmat <- t(cbind(c(0,1,0,0,0),c(0,0,0,1,0),c(0,0,0,1,0),c(0,0,0,0,1),
                c(0,0,0,0,0)))
colnames(wmat) <- rownames(wmat) <- c("1","2","3","4","5")
g <- as(wmat,"graphNEL")

e.true <- 0
var.true <- 5

dat <- rmvDAG(1000,g)
x5 <- dat[,5]

## test mean
if (t.test(x5,alternative="two.sided")$p.value<0.05) {
  stop("Test of rmvDAG: Mean not correct!")
}

## test variance
if (var.test(x5,rnorm(1000,0,sqrt(5)),ratio=1,
             alternative="two.sided")$p.value<0.05) {
  stop("Test of rmvDAG: Variance not correct!")
}
  

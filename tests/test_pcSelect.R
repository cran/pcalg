library(pcalg)

p <- 10
n <- 1000

set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
res1 <- pcSelect(d.mat[,10],d.mat[,-10],alpha=0.05)
resTrue <- rep(FALSE,p-1)
resTrue[c(4,5,6)] <- TRUE
if (!all(res1$G==resTrue)) {
  stop("Test of pcSelect: Consistency problem 101")
}

set.seed(456)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
res1 <- pcSelect(d.mat[,10],d.mat[,-10],alpha=0.05)
resTrue <- rep(FALSE,p-1)
resTrue[c(5,8,9)] <- TRUE
if (!all(res1$G==resTrue)) {
  stop("Test of pcSelect: Consistency problem 101")
}

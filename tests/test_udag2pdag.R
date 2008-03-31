library(pcalg)

## acyclic graphs
nreps <- 30
p <- 8
n <- 1000

my.seed <- cyc.res <- rep(0,nreps)
for (i in 1:nreps) {
  set.seed(i)
  myDAG <- randomDAG(p, prob = 0.2)
  d.mat <- rmvDAG(n, myDAG, errDist = "normal")
  res <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard",directed=TRUE)
  res.amat <- wgtMatrix(res@graph)
  res.amat[res.amat!=0] <- 1
  undir.amat <- res.amat+t(res.amat)
  undir.amat[undir.amat==1] <- 0
  undir.amat[undir.amat==2] <- 1
  res.dir <- res.amat-undir.amat

  cyc.res[i] <- isAcyclic(res.dir)
}

if (any(!cyc.res)) {
  stop("Test of udag2pdag: Cyclic part in PDAG!")
}

## find collider correctly
set.seed(123)
myDAG <- randomDAG(3, prob = 0.5)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
res <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard",directed=TRUE)
gEst <- wgtMatrix(res@graph)
gTrue <- cbind(c(0,0,1),c(0,0,1),c(0,0,0))
if (!all(gEst==gTrue)) {
  stop("Test of udag2pdag: Problem finding a collider!")
}

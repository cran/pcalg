library(pcalg)

set.seed(234)
p <- 10
nreps <- 20
resG <- resP <- rep(FALSE,nreps)
for (i in 1:nreps) {
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.3)
  mycor <- cov2cor(trueCov(myDAG))
  amat <- wgtMatrix(myDAG)
  amat[amat!=0] <- 1
  amat <- amat + t(amat)
  amat[amat!=0] <- 1

  ## Gaussian
  suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
  indepTest <- gaussCItest 

  resU <- skeleton(suffStat, indepTest, p, 0.99)

  resG[i] <- all(as(resU@graph,"matrix") == amat)
  resP[i] <- all(resU@pMax[as(resU@graph,"matrix") == TRUE] < 0.99)
}
if (!all(resG)) stop("Test skeleton wrong: Some skeleton was not found correctly!")
if (!all(resP)) stop("Test skeleton wrong: There was an inconsistency with an entry in pMax!")

library(pcalg)

nreps <- 10
check.res <- rep(FALSE,nreps)
set.seed(234)
for (ii in 1:nreps) {
  p <- 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.3)
  myCPDAG <- dag2cpdag(myDAG)

  suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
  indepTest <- gaussCItest 

  res <- pc(suffStat, indepTest, p, 0.99)

  check.res[ii] <- (shd(res,myCPDAG)==0)
}

if (!all(check.res)) stop("Test pc wrong: Some CPDAG was not found correctly!")

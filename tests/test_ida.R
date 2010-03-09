library(pcalg)

set.seed(123)
nreps <- 100
res <- rep(FALSE,nreps)
all.eff.true <- res
for (i in 1:nreps) {
  p <- 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.2)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)

  x <- sample(1:10,1)
  y <- sample(setdiff(1:10,x),1)

  ## plot(myCPDAG)
  eff.true <- round(causalEffect(myDAG, y, x),14)
  all.eff.true[i] <- eff.true
  ## cat("x=",x," y=",y," eff=",eff.true,"\n")

  eff.est <- round(ida(x,y,mcov,myCPDAG,method="local",verbose=FALSE),14)

  res[i] <- (eff.true %in% eff.est)
}

if (!all(res)) stop("Test ida: True effects were not recovered!")

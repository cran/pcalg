library(pcalg)

p <- 10
set.seed(21)
reps <- 10
res <- rep(FALSE,reps)
stopifnot(require("ggm"))# e.g. isAcyclic() below
for (i in 1:reps) {
  amat <- matrix(sample(c(0,1),p*p,replace=TRUE),p,p)
  diag(amat) <- rep(0,p)
  skelT <- amat+t(amat)
  skelT[skelT!=0] <- 1

  ## same skeleton?
  my.dag <- amat2dag(amat)
  skelDAG <- my.dag+t(my.dag)
  res1 <- all(skelDAG==skelT)

  ## acyclic?
  res2 <- isAcyclic(my.dag)

  res[i] <- (res1&res2)
}

if(!all(res)) {
  stop("Test amat2dag: Problem")
}

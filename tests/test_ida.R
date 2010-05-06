library(pcalg)

set.seed(123)
nreps <- 100
res <- rep(FALSE,nreps)
all.eff.true <- res
Rnd <- function(e) round(e, 14)
for (i in 1:nreps) {
  p <- 2 + rpois(1, lambda = 8) # ==>  p >= 2, E[p] = 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.2)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)

  x <- sample(1:p, 1)
  y <- sample(setdiff(1:p, x),1)

  ## plot(myCPDAG)
  eff.true <- Rnd(causalEffect(myDAG, y, x))
  all.eff.true[i] <- eff.true
  ## cat("x=",x," y=",y," eff=",eff.true,"\n")

  eff.est <- Rnd(ida(x,y, mcov, myCPDAG, method="local"))
  res[i] <- (eff.true %in% eff.est)
}
cat('Time elapsed: ', (.pt <- proc.time()),"\n")

stem(all.eff.true)
if (!all(res)) stop("Test ida: True effects were not recovered!")

## *one* test for  method="global" :
eff.g.est <- Rnd(ida(x,y, mcov, myCPDAG, method="global", verbose=TRUE))
stopifnot(identical(eff.est, eff.g.est))

cat('Time elapsed additionally: ', proc.time() - .pt,"\n")

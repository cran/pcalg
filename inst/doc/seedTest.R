##################################################
## Gauss: Seed = 123
##################################################
nreps <- 10000
gut <- rep(FALSE, nreps)
for (i in 1:nreps) {
  cat("i=",i,"\n")
  i <- 40
  set.seed(i)
  p <- 8
  n <- 5000
  gGtrue <- randomDAG(p, prob = 0.3)
  datG <- rmvDAG(n, gGtrue)

  indepTest <- gaussCItest 
  suffStat <- list(C = cor(datG), n = 5000)
  pc.fit <- pc(suffStat, indepTest, p = 8, alpha = 0.01)

      par(mfrow = c(1,2))
      plot(gGtrue, main = "True DAG")
      plot(pc.fit, main = "Estimated DAG")

  causalEffect(gGtrue, 6, 1)
  ida(1, 6, cov(datG), pc.fit@graph)
  idaFast(1, c(4,5,6), cov(datG), pc.fit@graph)
  
  if ((length(res) == 2) & (all(abs(res)>0.1)) & (abs(causalEffect(gGtrue, 6,1)) > 0.1)) {
##    gut[i] <- TRUE
  }
}



##################################################
## IDA: Seed = 123
##################################################
seedMax <- 1000
res <- rep(FALSE, seedMax)

for (i in 1:seedMax) {
  cat("i=",i,"\n")
  set.seed(i)
  p <- 7
  myDAG <- randomDAG(p, prob = 0.2) ## true DAG

  n <- 10000
  dat <- rmvDAG(n, myDAG)

  if (all(datI[1,] == dat[1,])) res[i] <- TRUE
}



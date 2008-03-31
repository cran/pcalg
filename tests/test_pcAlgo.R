library(pcalg)

seeds <- c(123,444,324,543,764,2316,3341,934,22335,343)
n <- 100000
p <- 6
a <- 0.001
en <- 3
correctEst <- rep(FALSE,length(seeds))

s <- en/(p-1)

for (i in 1:length(seeds)) {
  cat("i=",i,"\n")
  set.seed(seeds[i])
  myDAG <- randomDAG(p, prob = s)
  true.cpdag <- dag2cpdag(myDAG)
  d.mat <- rmvDAG(n, myDAG, errDist = "normal")
  res <- pcAlgo(d.mat, alpha = a, corMethod = "standard",directed=TRUE,verbose=0)
  compare.res <- compareGraphs(res@graph,true.cpdag)

  if ((compare.res["tpr"]==1) & (compare.res["fpr"]==0)) {
    correctEst[i] <- TRUE
  }
}

if (!all(correctEst)) stop("Test pcAlgo: Consistency wrong!")

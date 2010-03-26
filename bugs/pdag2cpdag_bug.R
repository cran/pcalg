library(pcalg)
load("/u/kalisch/R/bugs/pdag2cpdag_bug.dat")

## gPDAG1
par(mfrow=c(1,3))
plot(gPDAG1)
res1 <- pdag2dag(gPDAG1@graph)
plot(res1)
res2 <- dag2cpdag(res1)
plot(res2)

## gPDAG2
par(mfrow=c(1,3))
plot(gPDAG2)
res1 <- pdag2dag(gPDAG2@graph)
plot(res1)
res2 <- dag2cpdag(res1)
plot(res2)

## gPDAG3
par(mfrow=c(1,3))
plot(gPDAG3)
res1 <- pdag2dag(gPDAG3@graph)
plot(res1)
res2 <- dag2cpdag(res1)
plot(res2)


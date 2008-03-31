library(pcalg)

set.seed(123)

g <- randomDAG(5,0.3)
pdag <- dag2cpdag(g)
amatP <- wgtMatrix(pdag)
pdagTrue <- matrix(c(0,0,0,0,1, 0,0,1,0,1, 0,1,0,1,0, 0,0,1,0,0, 0,0,0,0,0),5,5)
if (!all(amatP==pdagTrue)) {
  stop("Test of dag2cpdag: True CPDAG not found!")
}

set.seed(567)

g <- randomDAG(5,0.3)
pdag <- dag2cpdag(g)
amatP <- wgtMatrix(pdag)
pdagTrue <- matrix(c(0,0,1,0,1, 0,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0, 1,0,0,0,0),5,5)
if (!all(amatP==pdagTrue)) {
  stop("Test of dag2cpdag: True CPDAG not found!")
}


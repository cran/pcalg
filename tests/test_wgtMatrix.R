library(pcalg)

set.seed(42)

## wmat_ij is edge from j to i
g <- randomDAG(3,0.4)
wmat <- wgtMatrix(g)
if (!(wmat[2,1]!=0 & wmat[1,2]==0 & wmat[3,1]!=0 & wmat[3,2]==0)) {
  stop("Test of wgtMatrix: Something with orientation of edges is wrong!")
}

## test weird parameters
g <- randomDAG(3,0)
wmat <- wgtMatrix(g)
if (!all(wmat==matrix(0,3,3))) {
  stop("Test of wgtMatrix: Problem when used on empty graph!")
}

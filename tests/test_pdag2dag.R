library(pcalg)

## expl 1
amat1 <- cbind(c(0,1),c(1,0))
g1 <- as(amat1,"graphNEL")
d1 <- pdag2dag(g1)
res1 <- all(cbind(c(0,0),c(1,0))==wgtMatrix(d1$graph))

## expl2
amat2 <- cbind(c(0,0,1),c(1,0,1),c(0,1,0))
g2 <- as(amat2,"graphNEL")
d2 <- pdag2dag(g2)
amatT2 <- cbind(c(0,1,0),c(0,0,0),c(1,1,0))
res2 <- all(amatT2==wgtMatrix(d2$graph))

## expl3 - not extendable
amat3 <- cbind(c(0,1,0,1),c(1,0,1,0),c(0,1,0,1),c(1,0,1,0))
g3 <- as(amat3,"graphNEL")
d3 <- pdag2dag(g3)
res3 <- !d3$success

if (!all(c(res1,res2,res3))) {
  stop("Test of pdag2dag: Problem when extending PDAG!")
}

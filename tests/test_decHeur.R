library(pcalg)

set.seed(123)

p <- 8
n <- 100

g1 <- randomDAG(p,0.2)
dag <- rmvDAG(n,g1,errDist="normal")
okN <- !decHeur(dag)$use.rob

dag <- rmvDAG(n,g1,errDist="mixt3")
okT <- !decHeur(dag)$use.rob

dag <- rmvDAG(n,g1,errDist="mix")
okC <- decHeur(dag)$use.rob

if (!all(c(okN,okT,okC))) {
  stop("Test of decHeur: Wrong decision!")
}

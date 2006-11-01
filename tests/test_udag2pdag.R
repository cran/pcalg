library(pcalg)

p <- 10 # number of random variables
n <- 10000 # number of samples
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
d <- rmvDAG(n,g) # generate random samples

gSkel <- 
  pcAlgo(d,alpha=0.05) # estimate of the skeleton
as(gSkel@graph,"matrix")

gPDAG <- udag2pdag(gSkel)
as(gPDAG@graph,"matrix")

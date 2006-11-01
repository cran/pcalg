library(pcalg)

p <- 7 # number of random variables
n <- 10000 # number of samples
s <- 0.3 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
d <- rmvDAG(n,g) # generate random samples

gSkel <- 
  pcAlgo(d,alpha=0.05) # estimate of the skeleton
(as(gSkel@graph,"matrix"))
gCPDAG <- 
  udag2cpdag(gSkel) # transform skeleton to cpdag
(as(gCPDAG@graph,"matrix"))


library(pcalg)

p <- 10 # number of random variables
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
res <- dag2cpdag(g)
as(res,"matrix")

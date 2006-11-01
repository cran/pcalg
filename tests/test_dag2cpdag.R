library(pcalg)

p <- 10 # number of random variables
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG

## transform graph to adjacency matrix
amat <- as(g,"matrix")
amat[amat!=0] <- 1
res <- dag2cpdag(amat)
as(res,"matrix")

library(pcalg)

n.perm <- 5

set.seed(123)

# A: adjacency matrix of DAG;
# B: adjacency matrix of CPDAG
# Setting 3 by courtesy of Jonas Peters: in pcalg <= 2.0.8,
# setting i = 3, k = 3 failed
A <- list(
  matrix(c(0,0,0,0,1, 0,0,1,0,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,0,0,0, 0,0,0,1,0, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,0, 1,1,1,0), 4, 4),
  matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,0, 1,1,1,0), 4, 4),
  matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,0, 1,1,1,0), 4, 4))
B <- list(
  matrix(c(0,0,0,0,1, 0,0,1,0,1, 0,1,0,1,0, 0,0,1,0,0, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,0,0,0, 1,0,0,1,0, 0,0,0,1,0, 0,0,0,0,1, 0,0,0,0,0), 5, 5, byrow = TRUE),
  matrix(c(0,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0), 4, 4),
  matrix(c(0,1,1,0, 1,0,1,0, 1,1,0,0, 1,1,1,0), 4, 4),
  matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,1, 1,1,1,0), 4, 4))
targets <- list(
  list(integer(0)),
  list(integer(0)),
  list(integer(0)),
  list(integer(0), 4),
  list(integer(0), 2))

for (i in 1:length(A)) {
  for (k in 1:n.perm) {
    p <- nrow(A[[i]])
    
    ind <- if(k == 1) 1:p else sample.int(p)
    permTargets <- lapply(targets[[i]], function(v) match(v, ind))
    
    # Test functionality with a matrix
    B.hat <- dag2essgraph(A[[i]][ind, ind], targets = permTargets)
    if (!all(B.hat == B[[i]][ind, ind])) {
      stop(sprintf("True CPDAG not found! (setting: matrix, i = %d, k = %d)", i, k))
    }
    
    # Test functionality with grephNEL objects
    g <- as(A[[i]][ind, ind], "graphNEL")
    pdag <- dag2essgraph(g, targets = permTargets)
    B.hat <- as(pdag, "matrix")
    if (!all(B.hat == B[[i]][ind, ind])) {
      stop(sprintf("True CPDAG not found! (setting:  graphNEL, i = %d, k = %d)", i, k))
    }
    
    # Test functionality with ParDAG/EssGraph objects
    g <- as(A[[i]][ind, ind], "GaussParDAG")
    pdag <- dag2essgraph(g, targets = permTargets)
    B.hat <- as(pdag, "matrix")
    if (!all(B.hat == B[[i]][ind, ind])) {
      stop(sprintf("True CPDAG not found! (setting:  ParDAG, graphNEL, i = %d, k = %d)", i, k))
    }
    
    # par(mfrow = c(1, 2))
    # plot(g)
    # plot(pdag)
  }
}

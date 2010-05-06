library(pcalg)

p <- 6 # number of random variables
n <- 10000 # number of samples
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
g2 <- dag2cpdag(g)
plot(g2)

plot(pdag2dag(g2)$graph)

am2 <- as(g2,"matrix")
g2neu <- as(am2,"graphNEL")

rownames(am2) <- colnames(am2) <- c("a","b","c","A","B","C")
g3 <- as(am2,"graphNEL")
g3neu <- pdag2dag(g3)
plot(g3neu$graph)

pdag2dag <- function(g) {
  ## Purpose: Generate a consistent extension of a PDAG to a DAG; if this
  ## is not possible, a random extension of the skeleton is returned and
  ## a warning is issued.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g: PDAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:21

  if (numEdges(g)==0) {
    res <- g
  } else {
    gm <- wgtMatrix(g) ## gm_i_j is edge from j to i
    gm[which(gm>0 & gm!=1)] <- 1
    p <- dim(gm)[1]

    ## temporarily change rownames
    rowNames <- rownames(gm)
    rownames(gm) <- colnames(gm) <- as.character(1:p)
    
    gm2 <- gm
    a <- gm
    go.on <- TRUE
    go.on2 <- FALSE
    while(length(a)>1 & sum(a)>0 & go.on) {
      go.on <- FALSE
      go.on2 <- TRUE
      sinks <- find.sink(a)
      if (length(sinks)>0) {
        counter <- 1
        while(counter<=length(sinks) & go.on2==TRUE) {
          x <- sinks[counter]
          if (adj.check(a,x)) {
            go.on2 <- FALSE
            ## orient edges
            inc.to.x <- which(a[,x]==1 & a[x,]==1) ## undirected
            if (length(inc.to.x)>0) {
              real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
              real.x <- as.numeric(row.names(a)[x])
              gm2[real.x,real.inc.to.x] <- rep(1,length(inc.to.x))
              gm2[real.inc.to.x,real.x] <- rep(0,length(inc.to.x))
            }
            ## remove x and all edges connected to it
            a <- a[-x,-x]
          }
          counter <- counter+1
        }
      }
      go.on <- !go.on2
    }
    if (go.on2==TRUE) {
      rownames(gm) <- colnames(gm) <- rowNames
      res <- as(amat2dag(gm),"graphNEL")
      warning("PDAG not extendible: Random DAG on skeleton drawn")
      succ <- FALSE
    } else {
      rownames(gm2) <- colnames(gm2) <- rowNames
      res <- as(t(gm2),"graphNEL")
      succ <- TRUE
    }
  }
  list(graph=res,success=succ)
}


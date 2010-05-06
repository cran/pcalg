library(pcalg)

## Sarah's Beispiel (8 Knoten)
# load graphNEL object (pdag G.test)
load('/u/gerster/Diplomarbeit/PCAlg/PDAG2DAG_bug/pdag.dat')

# pdag
par(mfrow=c(1,2))
plot(G.test)

# transform pdag to dag
G.dag <- pdag2dag(G.test)

## resulting 'DAG' is not acyclic...
## no warning is printed...
plot(G.dag)

## Einfacheres Beispiel (5 Knoten)
a.mat2 <- t(matrix(c(0,1,0,0,0, 0,0,0,1,0, 0,1,0,0,1, 0,1,0,0,1, 0,1,0,0,0),5,5))
colnames(a.mat2) <- rownames(a.mat2) <- as.character(1:5)
G.test2 <- as(a.mat2,"graphNEL")
par(mfrow=c(1,2))
plot(G.test2)
G.dag2 <- pdag2dag(G.test2)
plot(G.dag2)

## Noch einfacheres Beispiel (3 Knoten)
a.mat3 <- matrix(c(0,1,1,1,0,0,0,1,0),3,3)
colnames(a.mat3) <- rownames(a.mat3) <- as.character(1:3)
G.test3 <- as(a.mat3,"graphNEL")
plot(G.test3)
G.dag3 <- pdag2dag(G.test3)
plot(G.dag3)

##################################################
## Teste mit simulierten Daten
##################################################
p <- 8
## generate and draw random DAG :
## set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
## plot(myDAG, main = "randomDAG(10, prob = 0.2)")

## generate 1000 samples of DAG using standard normal error distribution
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton and CPDAG of given data
resD <- pcAlgo(d.mat, alpha = 0.05, corMethod =
               "standard",directed=TRUE)
par(mfrow=c(1,2))
plot(resD)
resD.dag <- pdag2dag(resD@graph)
plot(resD.dag)
##################################################
##################################################

find.sink <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none
  ## a sink my have incident indirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:28

  ## treat undirected edges
  undir.nbrs <- which((gm==t(gm) & gm==1),arr.ind=TRUE)
  gm[undir.nbrs] <- 0
  ## treat directed edges
  res <- which(apply(gm,1,sum)==0)
  res
}

adj.check <- function(gm,x) {
  ## Purpose:  Return "TRUE", if:
  ## For every vertex y, adj to x, with (x,y) undirected, y is adjacent to
  ## all the other vertices which are adjacent to x
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: adjacency matrix of graph
  ## - x: node number (number)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:27

  res <- TRUE
  ## find undirected neighbors of x
  un <- which(gm[x,]==1 & gm[,x]==1)
  if (length(un)>0) {
    for (i in 1:length(un)) {
      y <- un[i]
      adj.x <- setdiff(which(gm[x,]==1 | gm[,x]==1),y)
      adj.y <- setdiff(which(gm[y,]==1 | gm[,y]==1),x)
      axINay <- all(adj.x %in% adj.y)
      res <- (res & axINay)
    }
  }
  res
}

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
    gm <- as(g,"matrix")
    gm[which(gm>0 & gm!=1)] <- 1
    p <- dim(gm)[1]
    
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
              gm2[real.inc.to.x,real.x] <- rep(1,length(inc.to.x))
              gm2[real.x,real.inc.to.x] <- rep(0,length(inc.to.x))
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
      res <- as(amat2dag(gm),"graphNEL")
      warning("PDAG not extendible: Random DAG on skeleton drawn")
    } else {
      res <- as(gm2,"graphNEL")
    }
  }
  res
}

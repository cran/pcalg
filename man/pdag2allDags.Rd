\name{pdag2allDags}
\alias{pdag2allDags}
\title{Enumerate All DAGs in a Markov Equivalence Class}
\description{
  \code{pdag2allDags} computes all DAGs in the Markov Equivalence Class
  Represented by a Given Partially Directed Acyclic Graph (PDAG).
}
\usage{
pdag2allDags(gm, verbose = FALSE)
}
\arguments{
  \item{gm}{adjacency matrix of type \link{amat.cpdag}}
  \item{verbose}{logical; if true, some output is produced during
    computation}
}
\details{
  All DAGs extending the given PDAG are computing while avoiding new
  v-structures and cycles. If no DAG is found, the function returns \code{NULL}.
}
\value{
  List with two elements:
\item{dags:}{Matrix; every row corresponds to a DAG; every column
  corresponds to an entry in the adjacency matrix of this DAG. Thus, the
adjacency matrix (of type \link{amat.cpdag}) contained in the i-th row
of matrix \code{dags} can be obtained by calling
\code{matrix(dags[i,],p,p, byrow = TRUE)} (assuming the input PDAG has
\code{p} nodes).} 
\item{nodeNms}{Node labels of the input PDAG.}
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\examples{
## Example 1
gm <- rbind(c(0,1),
            c(1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:2]
res1 <- pdag2allDags(gm)
## adjacency matrix of first DAG in output
amat1 <- matrix(res1$dags[1,],2,2, byrow = TRUE)
colnames(amat1) <- rownames(amat1) <- res1$nodeNms
amat1 ## A --> B

## Example 2
gm <- rbind(c(0,1,1),
            c(1,0,1),
            c(1,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res2 <- pdag2allDags(gm)
## adjacency matrix of first DAG in output
amat2 <- matrix(res2$dags[1,],3,3, byrow = TRUE)
colnames(amat2) <- rownames(amat2) <- res2$nodeNms
amat2

## Example 3
gm <- rbind(c(0,1,1,0,0),
            c(1,0,0,0,0),
            c(1,0,0,0,0),
            c(0,1,1,0,1),
            c(0,0,0,1,0))
colnames(gm) <- rownames(gm) <- LETTERS[1:ncol(gm)]
res3 <- pdag2allDags(gm)
## adjacency matrix of first DAG in output
amat3 <- matrix(res3$dags[1,],5,5, byrow = TRUE)
colnames(amat3) <- rownames(amat3) <- res3$nodeNms
amat3

if (require(Rgraphviz)) {
## for convenience a simple plotting function
## for the function output
plotAllDags <- function(res) {
    require(graph)
    p <- sqrt(ncol(res$dags))
    nDags <- ceiling(sqrt(nrow(res$dags)))
    par(mfrow = c(nDags, nDags))
    for (i in 1:nrow(res$dags)) {
        tmp <- matrix(res$dags[i,],p,p)
        colnames(tmp) <- rownames(tmp) <- res$nodeNms
        plot(as(tmp, "graphNEL"))
    }
}
plotAllDags(res1)
amat1 ## adj.matrix corresponding to the first plot for expl 1
plotAllDags(res2)
amat2 ## adj.matrix corresponding to the first plot for expl 2
plotAllDags(res3)
amat3 ## adj.matrix corresponding to the first plot for expl 3
}
}

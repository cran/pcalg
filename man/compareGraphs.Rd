\name{compareGraphs}
\alias{compareGraphs}
\title{Compare two graphs in terms of TPR, FPR and TDR}
\description{
  Compares the true undirected graph with an estimated
  undirected graph in terms of True Positive Rate (TPR), False Positive
  Rate (FPR) and True Discovery Rate (TDR).
}
\usage{
compareGraphs(gl, gt)
}
\arguments{
  \item{gl}{estimated graph (graph object)}
  \item{gt}{true graph (graph object)}
}
\value{a named numeric vector with three numbers
  \item{tpr}{True Positive Rate: Number of correctly found edges (in
    estimated graph) divided by number of true edges (in true graph)}
  \item{fpr}{False Positive Rate: Number of incorrectly found edges
    divided by number of true gaps}
  \item{tdr}{True Discovery Rate: Number of correctly found edges divided
    by number of found edges (both in estimated graph)}
}
\details{
  If the input graph is directed, the directions are omitted.
  Special cases:
  \itemize{
    \item If the true graph contains no edges, the tpr is defined to be
    zero.
    \item Similarly, if the true graph contains no gaps, the fpr is
    defined to be one.
    \item
    If there are no edges in the true graph and there are none in
    the estimated graph, tdr is one. If there are none in the true graph but
    there are some in the estimated graph, tdr is zero.
  }
}
\seealso{\code{\link{randomDAG}} for generating a random DAG.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler
}
\examples{
## generate a graph with 4 nodes
V <- LETTERS[1:4]
edL2 <- vector("list", length=4)
names(edL2) <- V
edL2[[1]] <- list(edges= 2)
edL2[[2]] <- list(edges= c(1,3,4))
edL2[[3]] <- list(edges= c(2,4))
edL2[[4]] <- list(edges= c(2,3))
gt <- new("graphNEL", nodes=V, edgeL=edL2, edgemode="undirected")

## change graph
gl <- addEdge("A","C", gt,1)

## compare the two graphs
par(mfrow=c(2,1))
plot(gt) ; title("True graph")
plot(gl) ; title("Estimated graph")
(cg <- compareGraphs(gl,gt))
\dontshow{
stopifnot(cg == c(1, 1/2, 4/5))
}
}
\keyword{graphs}

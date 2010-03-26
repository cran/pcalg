\name{plotSG}
\alias{plotSG}
\title{Plot the subgraph around a Specific Node in a Graph Object}
\description{
  Plots a subgraph for a specified starting node and a given graph.  The
  subgraph consists of those nodes that can be reached from the starting
  node by passing no more than a specified number of edges.
}
\usage{
plotSG(graphObj, y, dist, amat = NA, directed = TRUE, main = )
}
\arguments{
  \item{graphObj}{an \R object of class \code{\link[graph:graph-class]{graph}}.}
  \item{y}{starting node.}
  \item{dist}{distance of nodes included in subgraph from starting node \code{y}.}
  \item{amat}{adjacency matrix of skeleton graph (optional).}
  \item{directed}{\code{\link{logical}} indicating if the subgraph
    should be directed.}
  \item{main}{title to be used, with a sensible default; see \code{\link{title}}.}
}
\value{the desired subgraph is also plotted and returned (made
  \code{\link{invisible}}).
}

\details{
  Commencing at the starting point \code{y} the function looks for the
  neighbouring nodes.  Beginning with direct parents and children it
  will continue hierarchically through the distances to \code{y}.  If
  \code{directed} is true (as per default), the orientation of the edges
  is taken from the initial graph.

  The package \pkg{Rgraphviz} must be installed, and is used for the
  plotting.
}

\author{
  Daniel Stekhoven (\email{hoven@stat.math.ethz.ch})
}
\examples{
if (require(Rgraphviz)) {
## generate a random DAG:
p <- 10
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

## plot whole the DAG
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

op <- par(mfrow = c(3,2))
## plot the neighbours of node number 8 up to distance 1
plotSG(myDAG, 8, 1, directed = TRUE)
plotSG(myDAG, 8, 1, directed = FALSE)

## plot the neighbours of node number 8 up to distance 2
plotSG(myDAG, 8, 2, directed = TRUE)
plotSG(myDAG, 8, 2, directed = FALSE)

## plot the neighbours of node number 8 up to distance 3
plotSG(myDAG, 8, 3, directed = TRUE)
plotSG(myDAG, 8, 3, directed = FALSE)

## Note that the layout of the subgraph might be different than in the
## original graph, but the graph structure is identical
par(op)
}% only when Rgraphviz is present
}
\keyword{graphs}

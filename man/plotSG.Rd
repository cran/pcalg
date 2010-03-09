\name{plotSG}
\alias{plotSG}
\title{Plot the neighbourhood of a specific node in a graph object}
\description{
  This function plots a subgraph for a specified starting node and a
  given graph. The subgraph consists of those nodes, that can be reached
  from the starting node by passing no more than a specified number of edges.
}
\usage{
plotSG( graphObj, y, dist, amat = NA, directed = TRUE )
}
\arguments{
  \item{graphObj}{Graph object.}
  \item{y}{Starting node.}
  \item{dist}{Distance of nodes included in subgraph from starting node y.}
  \item{amat}{Adjacency matrix of skeleton graph (optional).}
  \item{directed}{Boolean; should the plotted subgraph be directed?}
}
\value{A graph object containing the subgraph}

\details{
  Commencing at the starting point 'y' the function looks for the
  neighbouring nodes. Beginning with direct parents and children it
  will continue hierarchically through the distances to 'y'. If
  'directed' is set to TRUE the orientation of the edges is taken from
  the initial graph.
}

\author{
  Daniel Stekhoven (\email{hoven@stat.math.ethz.ch})
}
\examples{
## generate a random DAG:
p <- 10
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

## plot whole the DAG
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

x11()
par(mfrow = c(3,2))
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
## orignial graph, but the graph structure is identical
}
\keyword{graphs}

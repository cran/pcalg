\name{plotSG}
\alias{plotSG}
\title{Plot the subgraph around a Specific Node in a Graph Object}
\description{
  Plots a subgraph for a specified starting node and a given graph.  The
  subgraph consists of those nodes that can be reached from the starting
  node by passing no more than a specified number of edges.
}
\usage{
plotSG(graphObj, y, dist, amat = NA, directed = TRUE,
       plot = requireNamespace("Rgraphviz"), main = ,
       cex.main = 1.25, font.main = par("font.main"), col.main=par("col.main"),
       \dots)
}
\arguments{
  \item{graphObj}{An \R object of class \code{\link[graph:graph-class]{graph}}.}
  \item{y}{(integer) position of the starting node in the adjacency matrix.}
  \item{dist}{Distance of nodes included in subgraph from starting node \code{y}.}
  \item{amat}{Precomputed adjacency matrix of type \link{amat.cpdag} (optional)}
  \item{directed}{\code{\link{logical}} indicating if the subgraph
    should be directed.}
  \item{plot}{logical indicating if the subgraph should be plotted (or
    just returned).  Defaults to true when \pkg{Rgraphviz} is installed.}
  \item{main}{title to be used, with a sensible default; see
    \code{\link{title}}.}
  \item{cex.main, font.main, col.main}{optional settings for the
    \code{main} title; see \code{\link{title}}.}
  \item{\dots}{optional arguments passed to the \code{plot()}
    from package \pkg{Rgraphviz}.}
}
\value{the desired subgraph is returned; invisibly, i.e., via
  \code{\link{invisible}}, if \code{plot} is true.
}

\details{
  Commencing at the starting point \code{y} the function looks for the
  neighbouring nodes.  Beginning with direct parents and children it
  will continue hierarchically through the distances to \code{y}.  Note
  that the neighbourhood does \emph{not} depend on edge directions.  If
  \code{directed} is true (as per default), the orientation of the edges
  is taken from the initial graph.

  For the plotting, the package \pkg{Rgraphviz} must be installed.
}

\author{
  Daniel Stekhoven, then Martin Maechler.
}
\examples{
## generate a random DAG:
p <- 10
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

if(requireNamespace("Rgraphviz")) {
  ## plot whole the DAG
  plot(myDAG, main = "randomDAG(10, prob = 0.3)")

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
} else { ## without 'Rgraphviz'
  sg2d <- plotSG(myDAG, 8, 2, directed = TRUE,  plot=FALSE)
  sg2u <- plotSG(myDAG, 8, 2, directed = FALSE, plot=FALSE)
}
}
\keyword{graphs}

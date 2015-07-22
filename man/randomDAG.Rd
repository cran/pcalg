\name{randomDAG}
\alias{randomDAG}
\title{Generate a Directed Acyclic Graph (DAG) randomly}
\description{
  Generate a random Directed Acyclic Graph (DAG).  The resulting graph
  is topologically ordered from low to high node numbers.
}
\usage{
randomDAG(n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
}
\arguments{
  \item{n}{Number of nodes, \eqn{n \ge 2}{N >= 2}.}
  \item{prob}{Probability of connecting a node to another node with
    higher topological ordering.}
  \item{lB, uB}{Lower and upper limit of edge weights, chosen uniformly
    at random, i.e., by \code{\link{runif}(., min=lB, max=uB)}.}
  \item{V}{\code{\link{character}} vector length \code{n} of node names.}
}
\value{
  An object of class \code{"graphNEL"}, see
  \code{\link[graph]{graph-class}} from package \pkg{graph}, with \code{n}
  named ("1" to "n") nodes and directed edges.  The graph is
  topologically ordered.
  Each edge has a weight between \code{lB} and \code{uB}.
}
\details{
  The \code{n} nodes are ordered.  Start with first node.  Let the
  number of nodes with higher order be k.  Then, the number of
  neighbouring nodes is drawn as Bin(k, \code{prob}).  The neighbours are
  then drawn without replacement from the nodes with higher order.  For
  each node, a weight is uniformly sampled from \code{lB} to \code{uB}.
  This procedure is repeated for the next node in the original ordering
  and so on.
}
\seealso{\code{\link{randDAG}} for a more elaborate version of this
  function; \code{\link{rmvDAG}} for generating data according to a
  DAG; \code{\link{compareGraphs}} for comparing the skeleton of a DAG
  with some other undirected graph (in terms of TPR, FPR and TDR).
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler
}
\examples{
set.seed(101)
myDAG <- randomDAG(n = 20, prob= 0.2, lB = 0.1, uB = 1)
## require(Rgraphviz)
plot(myDAG)
}
\keyword{graphs}
\keyword{datagen}

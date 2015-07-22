\name{wgtMatrix}
\alias{wgtMatrix}
\title{Weight Matrix of a Graph, e.g., a simulated DAG}
\description{
  Given a \code{\link[graph:graph-class]{graph}} object \code{g}, as
  generated e.g., by \code{\link{randomDAG}}, return the matrix of its
  edge weights, the \dQuote{weight matrix}.
}
\usage{
wgtMatrix(g, transpose = TRUE)
}
\arguments{
  \item{g}{\code{\link[graph:graph-class]{graph}} object (package
    \pkg{graph}) of, say, \eqn{p} nodes, e.g. containing a DAG.}
  \item{transpose}{logical indicating if the weight matrix should be
    transposed (\code{\link{t}(.)}, see details).}
}
\details{
  When generating a DAG (e.g. using \code{\link{randomDAG}}), a graph
  object is usually generated and edge weights are usually specified.
  This function extracts the edge weights and arranges them in a matrix
  \eqn{M}.

  If \code{transpose} is \code{TRUE} (default), \code{M[i,j]} is the
  weight of the edge from j to i. If \code{transpose} is false, M[i,j]
  is the weight of the edge from i to j.

  Nowadays, this is a trivial wrapper around \code{as(g, "matrix")}
  using the (\code{\link{coerce}}) method provided by the \pkg{graph}
  package.
}
\value{
  The \eqn{p \times p}{p * p} weight matrix \code{M}.
}
\author{
Markus Kalisch
}
\note{
  This function can \emph{not} be used to estimate the edge weights in
  an estimated DAG / CPDAG.
}

\seealso{
  \code{\link{randomDAG}} for generating a random DAG;
  \code{\link{rmvDAG}} for simulating data from a generated DAG.
}
\examples{
set.seed(123)
g <- randomDAG(n = 5, prob = 0.3) ## generate random DAG
if(require(Rgraphviz)) {
  plot(g)
}

## edge weights as matrix
wgtMatrix(g)

## for comparison: edge weights in graph object
g@edgeData@data
}
\keyword{misc}

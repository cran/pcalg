\name{getGraph}
\alias{getGraph}
\alias{getGraph-methods}
\alias{getGraph,ANY-method}
\alias{getGraph,fciAlgo-method}
\alias{getGraph,matrix-method}
\alias{getGraph,pcAlgo-method}
\title{Get the "graph" Part or Aspect of R Object}
\description{
  Get the \code{\link[graph]{graph-class}}  part or \dQuote{aspect} of an \R
  object, notably from our \code{\link{pc}()}, \code{\link{skeleton}()},
  \code{\link{fci}()}, etc, results.
}
\usage{
getGraph(x)
}
\arguments{
  \item{x}{potentially any \R object which can be interpreted as a graph
    (with nodes and edges).}
}
\value{
  a \code{\link[graph]{graph-class}} object, i.e., one inheriting from (the
  virtual) class \code{"graph"}, package \pkg{graph}.
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "ANY")}}{the default method just tries
      \code{as(x, "graph")},  so works when a \code{\link{coerce}} (S4)
      method is defined for \code{x}.}
    \item{\code{signature(x = "pcAlgo")}}{and}
    \item{\code{signature(x = "fciAlgo")}}{extract the graph part
      explicitly.}
    \item{\code{signature(x = "matrix")}}{interpret \code{x} as adjacency
      matrix and return the corresponding \code{"\link[graph]{graphAM}"}
      object.}
  }
  For sparseMatrix methods, see the \sQuote{Note}.
}
\note{For large graphs, it may be attractive to work with \bold{sparse
    matrices} from the \pkg{Matrix} package.
  If desired, you can activate this by \preformatted{%
 require(Matrix)
 setMethod("getGraph", "sparseMatrix", function(x) as(x, "graphNEL"))
 setMethod("getGraph", "Matrix",       function(x) as(x, "graphAM"))
  }
}
\author{Martin Maechler}
\seealso{
  \code{\link{fci}}, etc.
  The \code{\link[graph]{graph-class}} class description in package \pkg{graph}.
}
\examples{
A <- rbind(c(0,1,0,0,1),
           c(0,0,0,1,1),
           c(1,0,0,1,0),
           c(1,0,0,0,1),
           c(0,0,0,1,0))
sum(A) # 9
getGraph(A) ## a graph with 5 nodes and  'sum(A)' edges
}
\keyword{methods}
\keyword{graph}

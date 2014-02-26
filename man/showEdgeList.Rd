\name{showEdgeList}% FIXME: "too much" in parallel with showAmat.Rd % (notably example)
\alias{showEdgeList}
\title{Show Edge List of pcAlgo object}
\description{
  Show the list of edges (of the graph) of a pcAlgo object; this is
  intended to be an alternative if \pkg{Rgraphviz} does not work.
}
\usage{
showEdgeList(object, labels = NULL)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{pcAlgo}}, as returned
    from \code{\link{skeleton}()} or \code{\link{pc}()}.}
  \item{labels}{optional labels for nodes; by default, the labels from
    the \code{object} are used.}
}
\value{none; the purpose is in (the side effect of) printing
  the edge list.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\note{
  This is not quite ok for \code{"\linkS4class{fciAlgo}"} objects, yet.% FIXME: do the work
}
\seealso{
  \code{\link{showAmat}} for the adjacency matrix of a \code{\linkS4class{pcAlgo}} object.
  \code{\link{iplotPC}} for plotting a pcAlgo object using the package
  \pkg{igraph}, also for an example of \code{showEdgeList()}.
}
\keyword{graphs}

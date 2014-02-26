\name{showAmat}
\alias{showAmat}
\title{Show Adjacency Matrix of pcAlgo object}
\description{
  Show the  adjacency matrix of a \code{"pcAlgo"} object; this is
  intended to be an alternative if the \pkg{Rgraphviz} package does not work.
}
\usage{
showAmat(object)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{pcAlgo}}, as returned
    from \code{\link{skeleton}()} or \code{\link{pc}()}.}
}
\value{
  The adjacency matrix.
}
\author{
Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\note{
  For \code{"\linkS4class{fciAlgo}"} objects, the show method produces a similar result.
}
\seealso{
  \code{\link{showEdgeList}} for showing the edge list of a \code{\linkS4class{pcAlgo}} object.
  \code{\link{iplotPC}} for plotting a \code{"pcAlgo"} object using the package
  \pkg{igraph} also for an example of \code{showAmat()}.
}
\keyword{hplot}

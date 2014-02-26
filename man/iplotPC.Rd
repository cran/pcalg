\name{iplotPC}
\alias{iplotPC}
\title{Plotting a pcAlgo object using the package igraph}
\description{
  Notably, when the \pkg{Rgraphviz} package is not easily available,
  \code{iplotPC()} is an alternative for plotting a \code{"\linkS4class{pcAlgo}"} object,
  making use of package \pkg{igraph}.

  It extracts the adjacency matrix and converts it into an
  object from package \pkg{igraph} which is then plotted.
}
\usage{
iplotPC(pc.fit, labels = NULL)
}
\arguments{
  \item{pc.fit}{an \R object of class \code{\linkS4class{pcAlgo}}, as returned
    from \code{\link{skeleton}()} or \code{\link{pc}()}.}
  \item{labels}{optional labels for nodes; by default, the labels from
    the \code{pc.fit} object are used.}
}
\value{
  Nothing.  As side effect, the plot of pcAlgo object \code{pc.fit}.
}
\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch}
}
\note{
  Note that this function does not work on \code{\linkS4class{fciAlgo}}
  objects, as those need different edge marks.
}
\seealso{
  \code{\link{showEdgeList}} for printing the edge list of a
  \code{\linkS4class{pcAlgo}} object; \code{\link{showAmat}} for
  printing the adjacency matrix of a pcAlgo object.
}
\examples{
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
pc.fit <- pc(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V, verbose = TRUE)

## Edge list
showEdgeList(pc.fit)

## Adjacency matrix
showAmat(pc.fit)

## Plot using package igraph; show estimated CPDAG:
iplotPC(pc.fit)
}
\keyword{graphs}

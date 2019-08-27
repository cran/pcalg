\name{pcAlgo-class}
\title{Class "pcAlgo" of PC Algorithm Results, incl. Skeleton}
\docType{class}
\alias{pcAlgo-class}
\alias{plot,pcAlgo,ANY-method}
\alias{show,pcAlgo-method}
\alias{summary,pcAlgo-method}
\alias{print.pcAlgo}
\description{This class of objects is returned by the functions
  \code{\link{skeleton}} and \code{\link{pc}} to represent the
  (skeleton) of an estimated CPDAG.
  Objects of this class have methods for the functions plot, show and
  summary.
}
\usage{% usage ..for methods with "surprising arguments":
\S4method{plot}{pcAlgo,ANY}(x, y, main = NULL,
     zvalue.lwd = FALSE, lwd.max = 7, labels = NULL, \dots)
\S3method{print}{pcAlgo}(x, amat = FALSE, zero.print = ".", \dots)

\S4method{summary}{pcAlgo}(object, amat = TRUE, zero.print = ".", \dots)
\S4method{show}{pcAlgo}(object)
}
\arguments{
  \item{x, object}{a \code{"pcAlgo"} object.}
  \item{y}{(generic \code{plot()} argument; unused).}
  \item{main}{main title for the plot (with an automatic default).}
  \item{zvalue.lwd}{\code{\link{logical}} indicating if the line width
    (\code{lwd}) of the edges should be made proportional to the entries
    of matrix \code{zMin} (originally) or derived from matrix \code{pMax}.}
  \item{lwd.max}{maximal \code{lwd} to be used, if \code{zvalue.lwd} is true.}
  \item{labels}{if non-\code{NULL}, these are used to define node attributes
    \code{nodeAttrs} and \code{attrs}, passed to \code{\link[Rgraphviz]{agopen}()}
    from package \pkg{Rgraphviz}.}
  \item{amat}{\code{\link{logical}} indicating if the adjacency matrix
    should be shown (printed) as well.}
  \item{zero.print}{string for printing \code{0} (\sQuote{zero}) entries
    in the adjacency matrix.}
  \item{\dots}{optional further arguments (passed from and to methods).}
}
\section{Creation of objects}{
  Objects are typically created as result from
  \code{\link{skeleton}()} or \code{\link{pc}()}, but could be
  be created by calls of the form \code{new("pcAlgo", ...)}.
}
\section{Slots}{
  The slots \code{call}, \code{n}, \code{max.ord}, \code{n.edgetests},
  \code{sepset}, and \code{pMax} are inherited from class
  \code{"\linkS4class{gAlgo}"}, see there.

  In addition, \code{"pcAlgo"} has slots
  \describe{
    \item{\code{graph}:}{Object of class \code{"\linkS4class{graph}"}:
      the undirected or partially directed graph that was estimated.}

    \item{\code{zMin}:}{Deprecated.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{gAlgo}"}.
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "pcAlgo")}: Plot the resulting
      graph. If argument \code{"zvalue.lwd"} is true, the
      linewidth an edge reflects \code{zMin}, so that
      thicker lines indicate more reliable dependencies. The argument
      \code{"lwd.max"} controls the maximum linewidth.}
    \item{show}{\code{signature(object = "pcAlgo")}: Show basic properties of
      the fitted object}
    \item{summary}{\code{signature(object = "pcAlgo")}: Show details of
      the fitted object }
  }
}
\author{Markus Kalisch and Martin Maechler}

\seealso{
  \code{\link{pc}}, \code{\link{skeleton}},
  \code{\linkS4class{fciAlgo}}
}

\examples{
  showClass("pcAlgo")

## generate a pcAlgo object
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
suffStat <- list(C = cor(d.mat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p = p)

## use methods of class pcAlgo
show(pc.fit)
if(require(Rgraphviz))
  plot(pc.fit, main = "Fitted graph")
summary(pc.fit)

## access slots of this object
(g  <- pc.fit@graph)
str(ss <- pc.fit@sepset, max=1)
}
\keyword{classes}

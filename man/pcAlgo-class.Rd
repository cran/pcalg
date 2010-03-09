\name{pcAlgo-class}
\docType{class}
\alias{pcAlgo-class}
\alias{plot,pcAlgo-method}
\alias{show,pcAlgo-method}
\alias{summary,pcAlgo-method}

\title{Class "pcAlgo"}
\description{This class of objects is returned by the functions 
  \code{\link{skeleton}} and \code{\link{pc}} to represent the 
  (skeleton) of an estimated CPDAG. 
  Objects of this class have methods for the functions plot, show and  
  summary.}
\section{Creation of objects}{
Objects can be created by calls of the form \code{new("pcAlgo", ...)} or
by function call to \code{\link{skeleton}} or \code{\link{pc}}.
}
\section{Slots}{
	 \describe{
    \item{\code{graph}:}{Object of the class \code{"graph"}:  
      the undirected or partially directed graph that was estimated.}
    \item{\code{call}:}{Object of the class \code{"call"}: the
      original function call.}
    \item{\code{n}:}{Object of the class \code{"integer"}: the
      sample size used to estimate the graph.}
    \item{\code{max.ord}:}{Object of class \code{"integer"}: the
      maximum size of the conditioning set used in the conditional 
      independence tests of the algorithm.}
    \item{\code{n.edgetests}:}{Object of class \code{"numeric"}: the number
      of conditional independence tests performed by the algorithm.}
    \item{\code{sepset}:}{Object of class \code{"list"}: the
      conditioning sets that led to edge deletions. The set that led to
    the removal of the edge \code{i-j} is saved in either
    sepset[[i]][[j]] or in sepset[[j]][[i]].}
    \item{\code{pMax}:}{Object of class \code{"matrix"}: the (i,j)th
      entry of the matrix contains the maximum p-value of all 
      conditional independence tests for edge i-j.}
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
	\code{\linkS4class{fci}}
      }

\examples{
  showClass("pcAlgo")

## generate a pcAlgo object
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
indepTest <- gaussCItest 
suffStat <- list(C = cor(d.mat), n = n)
alpha <- 0.01
pc.fit <- pc(suffStat, indepTest, p, alpha)

## access slots of this object
g <- pc.fit@graph
ss <- pc.fit@sepset
}

\keyword{classes}

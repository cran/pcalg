\name{fciAlgo-class}
\title{Class "fciAlgo"}
\docType{class}
\alias{fciAlgo-class}
\alias{plot,fciAlgo,ANY-method}
\alias{show,fciAlgo-method}
\alias{summary,fciAlgo-method}
\description{This class of objects is returned by the function
  \code{\link{fci}} to represent the estimated PAG.
  Objects of this class have methods for the functions
  \code{\link{plot}}, \code{\link{show}} and \code{\link{summary}}.
}
\section{Creation of objects}{
  Objects can be created by calls of the form \code{new("fciAlgo",
    ...)}, but are typically the result of \code{\link{fci}(..)}.
}
\section{Slots}{
  \describe{
    \item{\code{amat}:}{Object of class \code{"matrix"}: The
      the estimated graph, represented by its adjacency matrix.
      The edge marks are encoded
      by numbers: 0 = no edge, 1 = circle, 2 = arrowhead, 3 =
      tail. If \code{amat[i,j] = 1} and \code{amat[j,i] = 2},
      this represents the edge \code{i <-o j}.}
    \item{\code{call}:}{Object of class \code{"call"}: the
    original function call}
    \item{\code{n}:}{Object of class \code{"integer"}: the sample size
      used to estimate the graph.}
    \item{\code{max.ord}:}{Object of class \code{"integer"}: the maximum
      size of the conditioning sets used in the conditional indepedence tests
      in the first part of the algorithm (i.e., in the function
      \code{\link{skeleton}}).}
    \item{\code{n.edgetests}:}{Object of class \code{"numeric"}: the number
      of conditional independence tests performed in the first part of
      the algorithm.}
    \item{\code{sepset}:}{Object of class \code{"list"}: the
      conditioning sets that led to edge deletions. The set that led to
      the removal of the edge \code{i-j} is saved in either
      sepset[[i]][[j]] or in sepset[[j]][[i]].}
    \item{\code{pMax}:}{Object of class \code{"matrix"}:  the (i,j)th
      entry of the matrix contains the maximum p-value of all
      conditional independence tests for edge i-j.}
    \item{\code{allPdsep}}{Object of class \code{"list"}: the ith entry of
      this list contains Possible D-SEP of node \code{i}.}
    \item{\code{n.edgetestsPDSEP}}{Object of class \code{"numeric"}: the number
      of new conditional independence tests (i.e., tests that were not
      done in the first part of the algorithm) that were performed
      while checking subsets of Possible D-SEP.}
    \item{\code{max.ordPDSEP}}{Object of class \code{"integer"}: the maximum
      size of the conditioning sets used in the new conditional independence
      that were performed when checking subsets of Possible D-SEP.}
  }
}

\section{Extends}{
  Class \code{"\linkS4class{gAlgo}"}.
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "fciAlgo")}: Plot the resulting graph}
    \item{show}{\code{signature(object = "fciAlgo")}: Show basic properties of
      the fitted object}
    \item{summary}{\code{signature(object = "fciAlgo")}: Show details of
      the fitted object }
  }
}
\author{Markus Kalisch and Martin Maechler}

\seealso{
  \code{\link{fci}}, \code{\linkS4class{pcAlgo}}
}

\examples{
## look at slots of the class
showClass("fciAlgo")
\dontrun{
## Suppose, fciObj is an object of class fciAlgo
## access slots by using the @ symbol
fciObj@amat   ## adjacency matrix
fciObj@sepset ## separation sets

## use show, summary and plot method
show(fciObj)
summary(fciObj)
plot(fciObj)
}

## Also look at the extensive examples in  ?fci !
}
\keyword{classes}

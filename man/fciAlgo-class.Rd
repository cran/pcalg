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
  The slots \code{call}, \code{n}, \code{max.ord}, \code{n.edgetests},
  \code{sepset}, and \code{pMax} are inherited from class
  \code{"\linkS4class{gAlgo}"}, see there.
  
  In addition, \code{"fciAlgo"} has slots
  \describe{
    \item{\code{amat}:}{a \code{\link{matrix}}: The
      the estimated graph, represented by its adjacency matrix.
      The edge marks are encoded
      by numbers: 0 = no edge, 1 = circle, 2 = arrowhead, 3 =
      tail. If \code{amat[i,j] = 1} and \code{amat[j,i] = 2},
      this represents the edge \code{i <-o j}.}
    \item{\code{allPdsep}}{a \code{\link{list}}: the ith entry of
      this list contains Possible D-SEP of node number \code{i}.}
    \item{\code{n.edgetestsPDSEP}}{the number of new conditional
      independence tests (i.e., tests that were not done in the first
      part of the algorithm) that were performed while checking subsets
      of Possible D-SEP.}
    \item{\code{max.ordPDSEP}}{an \code{\link{integer}}: the maximum
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

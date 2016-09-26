\name{fciAlgo-class}
\title{Class "fciAlgo" of FCI Algorithm Results}
\docType{class}
\alias{fciAlgo-class}
\alias{plot,fciAlgo,ANY-method}
\alias{show,fciAlgo-method}
\alias{summary,fciAlgo-method}
\alias{print.fciAlgo}
\description{This class of objects is returned by functions
  \code{\link{fci}()}, \code{\link{rfci}()}, \code{\link{fciPlus}}, and
  \code{\link{dag2pag}} and represent the estimated PAG (and sometimes
  properties of the algorithm).
  Objects of this class have methods for the functions
  \code{\link{plot}}, \code{\link{show}} and \code{\link{summary}}.
}
%% Not interesting:
%% \section{Creation of objects}{
%%   Objects can be created by calls of the form \code{new("fciAlgo",
%%     ...)}, but are typically the result of \code{\link{fci}(..)}.
%% }
\section{Slots}{
  The slots \code{call}, \code{n}, \code{max.ord}, \code{n.edgetests},
  \code{sepset}, and \code{pMax} are inherited from class
  \code{"\linkS4class{gAlgo}"}, see there.

  In addition, \code{"fciAlgo"} has slots
  \describe{
    \item{\code{amat}:}{adjacency matrix; for the coding of the
  adjacency matrix see \link{amatType}}
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
\usage{% usage ..for methods with "surprising arguments":
\S4method{show}{fciAlgo}(object)
\S3method{print}{fciAlgo}(x, amat = FALSE, zero.print = ".", \dots)

\S4method{summary}{fciAlgo}(object, amat = TRUE, zero.print = ".", \dots)
\S4method{plot}{fciAlgo,ANY}(x, y, main = NULL, \dots)
}
\arguments{
  \item{x, object}{a \code{"fciAlgo"} object.}
  \item{amat}{\code{\link{logical}} indicating if the adjacency matrix
    should be shown (printed) as well.}
  \item{zero.print}{string for printing \code{0} (\sQuote{zero}) entries
    in the adjacency matrix.}
  \item{y}{(generic \code{plot()} argument; unused).}
  \item{main}{main title, not yet supported.}% <-> Rgraphviz support
  \item{\dots}{optional further arguments (passed from and to methods).}
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
  \code{\link{fci}}, \code{\link{fciPlus}}, etc (see above);
  \code{\linkS4class{pcAlgo}}
}

\examples{
## look at slots of the class
showClass("fciAlgo")

## Also look at the extensive examples in ?fci , ?fciPlus, etc !

\dontrun{
## Suppose, fciObj is an object of class fciAlgo
## access slots by using the @ symbol
fciObj@amat   ## adjacency matrix
fciObj@sepset ## separation sets

## use show, summary and plot method
fciObj ## same as  show(fciObj)
show(fciObj)
summary(fciObj)
plot(fciObj)
}%dont
}
\keyword{classes}

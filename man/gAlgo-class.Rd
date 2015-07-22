\name{gAlgo-class}
\docType{class}
\alias{gAlgo-class}
\title{Class \code{"gAlgo"}}
\description{
  \code{"gAlgo"} is a \code{"VIRTUAL"} class, the common basis of classes
  \code{"\linkS4class{pcAlgo}"} and \code{"\linkS4class{fciAlgo}"}.

  We describe the common slots here; for more see the help pages of the
  specific classes.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{a \code{\link{call}} object: the original function call.}
    \item{\code{n}:}{an \code{"integer"}, the sample size used to estimate the graph.}
    \item{\code{max.ord}:}{an \code{\link{integer}}, the maximum size of
      the conditioning set used in the conditional independence tests of
      the (first part of the algorithm), in function \code{\link{skeleton}}.}
    \item{\code{n.edgetests}:}{the number of conditional independence
      tests performed by the (first part of the) algorithm.}
    \item{\code{sepset}:}{a \code{\link{list}}, the conditioning sets
      that led to edge deletions.  The set that led to the removal of
      the edge \code{i -- j} is saved in either \code{sepset[[i]][[j]]} or
      in \code{sepset[[j]][[i]]}.}
    \item{\code{pMax}:}{a numeric square \code{\link{matrix}}, where the
      \eqn{(i,j)}th entry contains the maximal p-value of all conditional
      independence tests for edge \eqn{i -- j}.}
  }
}
%% \section{Methods}{
%% No methods defined with class "gAlgo" in the signature.
%% }
\author{Martin Maechler}
\seealso{\code{"\linkS4class{pcAlgo}"} and \code{"\linkS4class{fciAlgo}"}.
}
\examples{
showClass("gAlgo")
}
\keyword{classes}

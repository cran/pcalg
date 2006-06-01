\name{pcAlgo-class}
\docType{class}
\alias{pcAlgo-class}
\alias{plot,pcAlgo-method}
\alias{show,pcAlgo-method}
\alias{summary,pcAlgo-method}

\title{Class "pcAlgo"}
\description{The object contains the results of a call of the function
  \code{\link{pcAlgo}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pcAlgo", ...)} or
by a call of the function \code{\link{pcAlgo}}.
}
\section{Slots}{
	 \describe{
    \item{\code{graph}:}{Object of class \code{"graph"} which contains
      the undirected graph that was estimated }
    \item{\code{call}:}{Object of class \code{"call"} which contains the
    original call}
    \item{\code{n}:}{Object of class \code{"integer"}; number of samples
    that were used to estimate the graph}
    \item{\code{max.ord}:}{Object of class \code{"integer"}; maximum
      size of neighbourhood set that was conditioned on during the algorithm }
    \item{\code{n.edgetests}:}{Object of class \code{"numeric"}; number
      of edge tests that were performed }
    \item{\code{sepset}:}{Object of class \code{"list"}; the
      conditioning sets that made an edge drop out during the run of the
    algorithm}
    \item{\code{zMin}:}{Object of class \code{"matrix"}; contains at
      entry ij the minimum z-Value of all conditional independence tests
    of edge ij}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "pcAlgo")}: Plot the resulting
      undirected graph; if argument \code{"zvalue.lwd"} is true, the
      linewidth of the edges reflects the minimum absolut z-Value that
      was computed during all tests for the particular edge, i.e.,
      thicker lines show more reliable dependencies; the argument
      \code{"lwd.max"} controls the maximal linewidth.}
    \item{show}{\code{signature(object = "pcAlgo")}: Show properties of
      the fitted object}
    \item{summary}{\code{signature(object = "pcAlgo")}: Show details of
      the fitted object }
	 }
}
\author{Markus Kalisch and Martin Maechler}

\seealso{
	\code{\link{pcAlgo}}
}
% \examples{
% ## generate a pcAlgo-object using the pcAlgo function:
% p <- 10
% ## generate and draw random DAG :
% set.seed(101)
% class(myDAG <- randomDAG(p, prob = 0.2))
% plot(myDAG, main = "randomDAG(10, prob = 0.2)")
% ## generate 1000 samples of DAG using standard normal error distribution
% n <- 1000
% d.mat <- rmvDAG(n, myDAG, errDist = "normal")

% ## estimate skeleton given data
% res <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard")
% res # this is a pcAlgo object
% show(res)
% summary(res)
% plot(res)# << using the plot() method for 'pcAlgo' objects!
% plot(res, zvalue.lwd=TRUE) # indicate reliable edges by linewidth
% plot(res, zvalue.lwd=TRUE, lwd.max=20) # very thick lines
% }
\keyword{classes}

\name{pcAlgo}
\alias{pcAlgo}
\title{PC-Algorithm: Estimate the Underlying Graph (Skeleton) or
  Equivalence Class (CPDAG) of a DAG}
\description{
  This function is DEPRECATED! Use \code{\link{skeleton}}, \code{\link{pc}} or
  \code{\link{fci}} instead.
}
\usage{
pcAlgo(dm = NA, C = NA, n=NA, alpha, corMethod = "standard",
       verbose=FALSE, directed=FALSE, G=NULL, datatype = "continuous",
       NAdelete=TRUE, m.max=Inf, u2pd = "rand", psepset=FALSE)
}
\arguments{
  \item{dm}{data matrix; rows correspond to samples, cols correspond to
    nodes.}
  \item{C}{correlation matrix; this is an alternative for specifying the
    data matrix.}
  \item{n}{sample size; this is only needed if the data matrix is not
    provided.}
  \item{alpha}{significance level for the individual partial correlation tests.}
  \item{corMethod}{a character string speciyfing the method for
    (partial) correlation estimation.
    "standard", "QnStable", "Qn" or "ogkQn" for standard and robust (based on
    the Qn scale estimator without and with OGK) correlation
    estimation. For robust estimation, we recommend QnStable.}
  \item{verbose}{0-no output, 1-small output, 2-details;using 1 and 2
  makes the function very much slower}
  \item{directed}{If \code{FALSE}, the underlying skeleton is computed;
    if \code{TRUE}, the underlying CPDAG is computed}
  \item{G}{the adjacency matrix of the graph from which the algorithm
    should start (logical)}
  \item{datatype}{distinguish between discrete and continuous data}
  \item{NAdelete}{delete edge if pval=NA (for discrete data)}
  \item{m.max}{maximal size of conditioning set}
  \item{u2pd}{Function used for converting skeleton to cpdag. "rand"
    (use udag2pdag); "relaxed" (use udag2pdagRelaxed); "retry" (use
    udag2pdagSpecial)}
  \item{psepset}{If true, also possible separation sets are tested.}
}
\value{
  An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an undirected graph
  (object of \code{\link{class}} \code{"graph"}, see
  \code{\link[graph]{graph-class}} from the package \pkg{graph})
  (without weigths) as estimate of the skeleton or the CPDAG of the
  underlying DAG.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  Kalisch M. and P. B\"uhlmann (2007)
  \emph{Estimating high-dimensional
    directed acyclic graphs with the PC-algorithm};
    JMLR, Vol. 8, 613-636, 2007.
}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler.
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}

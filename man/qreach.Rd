\name{qreach}
\alias{qreach}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute possible d-sep of node x}
\description{
  This function computes the Possible-D-Sep of node x in a partially
  directed acyclic graph (PDAG). It is intended
  for internal use only.
}
\usage{
qreach(x, amat, verbose = FALSE)
}
\arguments{
\item{x}{Column of node in adjacency matrix, of which possible d-sep is
  to be computed}
\item{amat}{Adjacency matrix as returned by
  amat[i,j] = 0 iff no edge btw i,j
  amat[i,j] = 1 iff i *-o j
  amat[i,j] = 2 iff i *-> j}
\item{verbose}{Details on output}
}

\value{
  Vector of column positions indicating the nodes in Possible-D-Sep of x.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{fci}} and \code{\link{pdsep}} which both use
  this function.}




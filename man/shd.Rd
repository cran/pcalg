\name{shd}
\alias{shd}
\title{Compute Structural Hamming Distance (SHD)}
\description{
  Compute the Structural Hamming Distance (SHD) between two graphs.  In
  simple terms, this is the number of edge insertions, deletions or
  flips in order to transform one graph to another graph.
}
\usage{
shd(g1,g2)
}
\arguments{
  \item{g1}{Graph object}
  \item{g2}{Graph object}
}
\value{
  The value of the SHD (numeric).
}
\references{
  I. Tsamardinos, L.E. Brown and C.F. Aliferis (2006).
  The Max-Min Hill-Climbing Bayesian Network Structure Learning
  Algorithm. \emph{JMLR} \bold{65}, 31--78.
}

\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch} and Martin Maechler
}
\examples{
## generate two graphs
g1 <- randomDAG(10, prob = 0.2)
g2 <- randomDAG(10, prob = 0.2)
## compute SHD
(shd.val <- shd(g1,g2))
}
\keyword{graphs}


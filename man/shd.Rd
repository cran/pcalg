\name{shd}
\alias{shd}
\title{Compute Structural Hamming Distance (SHD)}
\description{
  Computer the Structural Hamming Distance between two graphs. In simple
  terms, this is the number of edge instertion, deletions or flips in
  order to transform one graph to another graph.
}
\usage{
shd(g1,g2)
}
\arguments{
  \item{g1}{graph object}
  \item{g2}{graph object}
}
\value{
  The value of the SHD (numeric).
}
\details{
  The "standard" method envokes a standard correlation estimator. "Qn"
  envokes a robust, elementwise correlation estimator based on the Qn scale
  estimte. "QnStable" also uses the Qn scale estimator, but uses an
  improved way of transforming that into the correlation
  estimator. "ogkQn" envokes a correlation estimator based on Qn using
  OGK. 
}
\seealso{
  \link{pcAlgo}
}

\references{
I. Tsamardinos, L.E. Brown and C.F. Aliferis (2006), "The Max-Min
Hill-Climbing Bayesian Network Structure Learning Algorithm"; JMLR,
Vol. 65, 31-78.
}

\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch} and Martin Maechler
}
\examples{
## generate two graphs
g1 <- randomDAG(10, prob = 0.2)
g2 <- randomDAG(10, prob = 0.2)
## comute SHD
shd.val <- shd(g1,g2)
}
\keyword{graphs}


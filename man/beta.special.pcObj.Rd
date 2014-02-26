\name{beta.special.pcObj}
\alias{beta.special.pcObj}
\title{Compute set of intervention effects in a fast way}
\description{
  This function is DEPRECATED! Use \code{\link{ida}} or
  \code{\link{idaFast}} instead.
}
\usage{
beta.special.pcObj(x.pos, y.pos, pcObj, mcov=NA, amat=NA, amatSkel=NA, t.amat=NA)
}
\arguments{
\item{x.pos}{Column of x in dat}
\item{y.pos}{Column of y in dat}
\item{pcObj}{Precomputed pc-object}
\item{mcov}{Covariance that was used in the pc-object fit}
\item{amat,  amatSkel,  t.amat}{Matrices that can be precomputed,  if
  needed (see code for details on how to precompute)}
}
\value{
  estimates of intervention effects
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pcAlgo}}, \code{\link{dag2cpdag}},
  \code{\link{beta.special}}}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}



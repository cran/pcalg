\name{pcorOrder}
\alias{pcorOrder}
\title{Compute Partial Correlations}
\description{
  This function computes partial correlations given a correlation matrix
  using a recursive algorithm or matrix inversion.
 }
\usage{
pcorOrder(i,j, k, C, cut.at = 0.9999999)
}
\arguments{
  \item{i,j}{integer variable numbers to compute partial correlations of.}
  \item{k}{conditioning set for partial correlations (vector of integers).}
  \item{C}{correlation matrix (matrix)}
  \item{cut.at}{number slightly smaller than one; if \eqn{c} is
    \code{cut.at}, values outside of \eqn{[-c,c]} are set to \eqn{-c} or
    \eqn{c} respectively.}
}
\value{
  The partial correlation of i and j given the set k.
}
\details{If the size of the conditioning sets is less than four, the
  partial correlations are computed recursively starting from the
  correlation matrix. If the size of the conditioning sets is greater
  than or equal to four, the partial correlations are computed by matrix
  inversion. This chooses in both cases the fastest method.
  To avoid numeric problems in recursions, the
  partial correlation is set to \code{cut.at} (\eqn{= 0.9999999} by
  default; keep this value, unless you know better!), if it is above
  (e.g., if 1) and to \code{-cut.at} if it is below (e.g., if -1).
}
\seealso{
  \code{\link{condIndFisherZ}} for testing zero partial correlation.
}
\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch} and Martin Maechler
}
\examples{
## produce uncorrelated normal random variables
mat <- matrix(rnorm(3*20),20,3)
## compute partial correlation of var1 and var2 given var3
pcorOrder(1,2, 3, cor(mat))

## define graphical model, simulate data and compute
## partial correlation with bigger conditional set
genDAG <- randomDAG(20, prob = 0.2)
dat <- rmvDAG(1000, genDAG)
C <- cor(dat)
pcorOrder(2,5, k = c(3,7,8,14,19), C)
}

\keyword{multivariate}

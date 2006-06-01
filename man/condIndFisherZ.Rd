\name{condIndFisherZ}
\alias{condIndFisherZ}
\alias{zStat}
\title{Conditional Independence by Fisher's Z-Transformation}
\description{
  Using Fisher's z-transformation of the partial correlation,
  test for zero partial correlation for sets of normally distributed
  random variables.
}
\usage{
condIndFisherZ(x, y, S, C, n, cutoff,
               verbose= isTRUE(getOption("verbose.pcalg.condIFz")))
zStat(x, y, S, C, n)
}
\arguments{
  \item{x,y,S}{it is tested, whether \code{x} and \code{y} are conditionally
    independent given the subset \code{S} of the remaining nodes.
    \code{x, y, S} all are integers, corresponding to variable or node
    numbers.}
  \item{C}{correlation matrix of nodes}
  \item{n}{integer specifying the number of observations
    (\dQuote{samples}) used to estimate the correlation matrix \code{C}.}
  \item{cutoff}{numeric cutoff for significance level of individual partial
    correlation tests.  Must be set to \code{\link{qnorm}(1 - alpha/2)}
    for a test significance level of \code{alpha}.}
  \item{verbose}{logical indicating whether some intermediate output should be
    shown (WARNING: This decreases the performance dramatically!)}
}
\value{
  \code{zStat()} gives a number
  \deqn{Z = \sqrt{n - \left|S\right| - 3} \cdot \log((1+r)/(1-r))/2}{%
        Z = sqrt(n - |S| - 3) * log((1+r)/(1-r))/2,}
  which is asymptotically normally distributed under the null hypothesis
  of correlation 0.

  \code{condIndFisherZ()} returns a \code{\link{logical}}
  \eqn{L} indicating whether the
  \emph{\dQuote{partial correlation of x and y given S is zero}}
  could not be rejected on the given significance level.  More
  intuitively and for multivariate normal data, this means: If
  \code{TRUE} then it seems plausible, that x and y are conditionally
  independent given S.  If \code{FALSE} then there was strong evidence
  found against this conditional independence statement.
}
\details{
  For gaussian random variables and after performing Fisher's
  z-transformation of the partial correlation, the test statistic
  \code{zStat()} is (asymptotically for large enough \code{n}) standard
  normally distributed.
  Partial correlation is tested in a two-sided hypothesis test, i.e.,
  basically, \code{condIndFisherZ(*) == abs(zStat(*)) > qnorm(1 - alpha/2)}.
  In a multivariate normal distribution, zero partial correlation is
  equivalent to conditional independence.
}
\references{
  Markus Kalisch and Peter B\"uhlmann (2005)
  \emph{Estimating high-dimensional directed acyclic graphs with the
    PC-algorithm}; Research Report Nr.~130, ETH Zurich;\cr
  \url{http://stat.ethz.ch/research/research_reports/2005}
}
\seealso{\code{\link{pcorOrder}} for computing a partial correlation
  given the correlation matrix in a recursive way.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler
}
\examples{
set.seed(42)
## Generate four independent normal random variables
n <- 20
data <- matrix(rnorm(n*4),n,4)
## Compute corresponding correlation matrix
corMatrix <- cor(data)
## Test, whether variable 1 (col 1) and variable 2 (col 2) are
## independent given variable 3 (col 3) and variable 4 (col 4) on 0.05
## significance level
x <- 1
y <- 2
S <- c(3,4)
n <- 20
alpha <- 0.05
cutoff <- 1-qnorm(alpha/2)
(b1 <- condIndFisherZ(x,y,S,corMatrix,n,cutoff))
   # -> 1 and 2 seem to be conditionally independent given 3,4

## Now an example with conditional dependence
data <- matrix(rnorm(n*3),n,3)
data[,3] <- 2*data[,1]
corMatrix <- cor(data)
(b2 <- condIndFisherZ(1,3,2,corMatrix,n,cutoff))
   # -> 1 and 3 seem to be conditionally dependent given 2
}
\keyword{htest}
\keyword{multivariate}


\name{condIndFisherZ}
\title{Test Conditional Independence of Gaussians via Fisher's Z}
\alias{condIndFisherZ}
\alias{zStat}
\alias{gaussCItest}
\description{
  Using Fisher's z-transformation of the partial correlation,
  test for zero partial correlation of sets of normally / Gaussian
  distributed random variables.

  The \code{gaussCItest()} function, using \code{zStat()} to test for
  (conditional) independence between gaussian random variables, with an
  interface that can easily be used in \code{\link{skeleton}},
  \code{\link{pc}} and \code{\link{fci}}.
}
\usage{
condIndFisherZ(x, y, S, C, n, cutoff, verbose= )
zStat         (x, y, S, C, n)
gaussCItest   (x, y, S, suffStat)
}
\arguments{
  \item{x,y,S}{(integer) position of variable \eqn{X}, \eqn{Y} and set
  of variables \eqn{S}, respectively, in the adjacency matrix. It is
  tested, whether \code{X} and \code{Y} are conditionally 
    independent given the subset \code{S} of the remaining nodes.}
  \item{C}{Correlation matrix of nodes}
  \item{n}{Integer specifying the number of observations
    (\dQuote{samples}) used to estimate the correlation matrix \code{C}.}
  \item{cutoff}{Numeric cutoff for significance level of individual partial
    correlation tests.  Must be set to \code{\link{qnorm}(1 - alpha/2)}
    for a test significance level of \code{alpha}.}
  \item{verbose}{Logical indicating whether some intermediate output should be
    shown; currently not used.}% (WARNING: This decreases the performance dramatically!)}
  \item{suffStat}{A \code{\link{list}} with two elements, \code{"C"}
    and \code{"n"}, corresponding to the above arguments with the same name.}
}
\value{
  \code{zStat()} gives a number
  \deqn{Z = \sqrt{n - \left|S\right| - 3} \cdot \log((1+r)/(1-r))/2}{%
        Z = sqrt(n - |S| - 3) * log((1+r)/(1-r))/2,}
  which is asymptotically normally distributed under the null hypothesis
  of correlation 0.

  \code{condIndFisherZ()} returns a \code{\link{logical}}
  \eqn{L} indicating whether the \emph{\dQuote{partial correlation of x
      and y given S is zero}} could not be rejected on the given
  significance level.  More intuitively and for multivariate normal
  data, this means:  If \code{TRUE} then it seems plausible, that x and
  y are conditionally independent given S.  If \code{FALSE} then there
  was strong evidence found against this conditional independence
  statement.

  \code{gaussCItest()} returns the p-value of the test.
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
  M. Kalisch and P. Buehlmann (2007).
  Estimating high-dimensional directed acyclic graphs with the PC-algorithm.
  \emph{JMLR} \bold{8} 613-636.
}
\seealso{\code{\link{pcorOrder}} for computing a partial correlation
  given the correlation matrix in a recursive way.

  \code{\link{dsepTest}}, \code{\link{disCItest}} and
  \code{\link{binCItest}} for similar functions for a d-separation
  oracle, a conditional independence test for discrete variables and a
  conditional independence test for binary variables, respectively.
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
cutoff <- qnorm(1-alpha/2)
(b1 <- condIndFisherZ(x,y,S,corMatrix,n,cutoff))
   # -> 1 and 2 seem to be conditionally independent given 3,4

## Now an example with conditional dependence
data <- matrix(rnorm(n*3),n,3)
data[,3] <- 2*data[,1]
corMatrix <- cor(data)
(b2 <- condIndFisherZ(1,3,2,corMatrix,n,cutoff))
   # -> 1 and 3 seem to be conditionally dependent given 2

## simulate another dep.case: x -> y -> z
set.seed(29)
x <- rnorm(100)
y <- 3*x + rnorm(100)
z <- 2*y + rnorm(100)
dat <- cbind(x,y,z)

## analyze data
suffStat <- list(C = cor(dat), n = nrow(dat))
gaussCItest(1,3,NULL, suffStat) ## dependent [highly signif.]
gaussCItest(1,3,  2,  suffStat) ## independent | S
}
\keyword{htest}
\keyword{multivariate}


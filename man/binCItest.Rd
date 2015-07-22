\name{binCItest}
\alias{gSquareBin}
\alias{binCItest}
\title{G square Test for (Conditional) Independence of Binary Variables}
\description{
  \eqn{G^2} test for (conditional) independence of
  \emph{binary} variables \eqn{X} and \eqn{Y} given the (possibly empty)
  set of binary variables \eqn{S}.

  \code{binCItest()} is a wrapper of \code{gSquareBin()}, to be easily
  used in \code{\link{skeleton}}, \code{\link{pc}} and
  \code{\link{fci}}.
}
\usage{
gSquareBin(x, y, S, dm, adaptDF = FALSE, n.min = 10*df, verbose = FALSE)
binCItest (x, y, S, suffStat)
}
\arguments{
  \item{x,y}{(integer) position of variable \eqn{X} and \eqn{Y},
    respectively, in the adjacency matrix.}
  \item{S}{(integer) positions of zero or more conditioning variables in the
    adjacency matrix.}
  \item{dm}{data matrix (with \eqn{\{0,1\}} entries).}
  \item{adaptDF}{logical specifying if the degrees of freedom should be
    lowered by one for each zero count.  The value for the degrees of
    freedom cannot go below 1.}
  \item{n.min}{the smallest \eqn{n} (number of observations,
    \code{nrow(dm)}) for which the G^2 test is computed; for smaller
    \eqn{n}, independence is assumed (\eqn{G^2 := 1}) with a warning.  The
    default is \eqn{10 m}, where \eqn{m} is the degrees of freedom
    assuming no structural zeros, \eqn{2^{|S|}}{2^|S|}.}
  \item{verbose}{logical or integer indicating that increased diagnostic output is to
    be provided.}

  \item{suffStat}{a \code{\link{list}} with two elements,
    \code{"dm"}, and \code{"adaptDF"} corresponding to the above two
    arguments of \code{gSquareBin()}.}
}
\value{
  The p-value of the test.
}
\details{
  The \eqn{G^2} statistic is used to test for (conditional)
  independence of X and Y given a set S (can be \code{NULL}).  This
  function is a specialized version of \code{\link{gSquareDis}} which is
  for discrete variables with more than two levels.
}
\references{
  R.E. Neapolitan (2004).
  Learning Bayesian Networks. \emph{Prentice Hall Series in Artificial
  Intelligence}. Chapter 10.3.1
}
\seealso{
  \code{\link{gSquareDis}} for a (conditional) independence test
  for discrete variables with more than two levels.

  \code{\link{dsepTest}}, \code{\link{gaussCItest}} and
  \code{\link{disCItest}} for similar functions for a d-separation
  oracle, a conditional independence test for Gaussian variables and
  a conditional independence test for discrete variables, respectively.

  \code{\link{skeleton}}, \code{\link{pc}} or \code{\link{fci}} which
  need a testing function such as \code{binCItest}.
}
\author{
  Nicoletta Andri and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
n <- 100
set.seed(123)
## Simulate *independent data of {0,1}-variables:
x <- rbinom(n, 1, pr=1/2)
y <- rbinom(n, 1, pr=1/2)
z <- rbinom(n, 1, pr=1/2)
dat <- cbind(x,y,z)

binCItest(1,3,2, list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.
binCItest(1,3,2, list(dm = dat, adaptDF = TRUE )) # the same, here

## Simulate data from a chain of 3 variables: x1 -> x2 -> x3
set.seed(12)
b0 <- 0
b1 <- 1
b2 <- 1
n <- 10000
x1 <- rbinom(n, size=1, prob=1/2) ## = sample(c(0,1), n, replace=TRUE)

## NB:  plogis(u) := "expit(u)" := exp(u) / (1 + exp(u))
p2 <- plogis(b0 + b1*x1) ; x2 <- rbinom(n, 1, prob = p2) # {0,1}
p3 <- plogis(b0 + b2*x2) ; x3 <- rbinom(n, 1, prob = p2) # {0,1}

ftable(xtabs(~ x1+x2+x3))
dat <- cbind(x1,x2,x3)

## Test marginal and conditional independencies
gSquareBin(3,1,NULL,dat, verbose=TRUE)
gSquareBin(3,1, 2,  dat)
gSquareBin(1,3, 2,  dat) # the same
gSquareBin(1,3, 2,  dat, adaptDF=TRUE, verbose = 2)
\dontshow{%% NB: more in --->>  ../tests/test_gSquareBin.R
stopifnot(all.equal(gSquareBin(3,1, 2, dat),
                    gSquareBin(1,3, 2, dat)),
          all.equal(gSquareBin(3,1, 2, dat, adaptDF=TRUE),
                    gSquareBin(1,3, 2, dat, adaptDF=TRUE)))
}%dont
}

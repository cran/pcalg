\name{disCItest}
\alias{disCItest}
\alias{gSquareDis}
\title{G square Test for (conditional) Independence for Discrete Data}
\description{
  \eqn{G^2} test for (conditional) independence of \emph{discrete}
  variables \eqn{X} and \eqn{Y} given the (possibly empty) set of
  discrete variables \eqn{S}. 

  \code{disCItest()} is a wrapper of \code{gSquareDis()}, to be easily
  used in \code{\link{skeleton}}, \code{\link{pc}} and \code{\link{fci}}.
}
\usage{
gSquareDis(x, y, S, dm, nlev, adaptDF = FALSE, verbose = FALSE)
disCItest (x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of variable X in the adjacency matrix}
  \item{y}{Position of variable Y in the adjacency matrix}
  \item{S}{Position of conditioning variables in the adjacency matrix}
  \item{dm}{Data matrix (rows: samples, columns: variables) with numeric
    entries; the k levels for a given column are coded by the integers
  0,1,...,k-1. (see example)}
  \item{nlev}{Vector with numbers of levels for each variable}
  \item{adaptDF}{Lower the degrees of freedom by one for each zero
    count.  The value for the degrees of freedom cannot go below 1.}
  \item{verbose}{Logical indicating if detailed output is to be provided.}

  \item{suffStat}{A list with three elements: \code{"dm"},
    \code{"nlev"}, \code{"adaptDF"}; each corresponding to the above
    arguments of \code{gSquareDis()}.}
}
\value{
  The p-value of the test.
}
\details{
  The \eqn{G^2} statistic is used to test for (conditional) independence
  of X and Y given a set S (can be \code{NULL}).  If only binary
  variables are involved, \code{\link{gSquareBin}} is a specialized
  (a bit more efficient) alternative to \code{gSquareDis()}.
} 
\seealso{
  \code{\link{gSquareBin}} for a (conditional) independence test
  for binary variables. 

  \code{\link{dsepTest}}, \code{\link{gaussCItest}} and
  \code{\link{binCItest}} for similar functions for a d-separation
  oracle, a conditional independence test for gaussian variables and a
  conditional independence test for binary variables, respectively.
}
\references{
  R.E. Neapolitan (2004).
  Learning Bayesian Networks. \emph{Prentice Hall Series in Artificial
  Intelligence}. Chapter 10.3.1
}
\author{
  Nicoletta Andri and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\examples{
## Simulate data
n <- 100
set.seed(123)
x <- sample(0:2, n, TRUE) ## three levels
y <- sample(0:3, n, TRUE) ## four levels
z <- sample(0:1, n, TRUE) ## two levels
dat <- cbind(x,y,z)

## Analyze data
gSquareDis(1,3,2,dat,nlev = c(3,4,2))
gSquareDis(1,3,2,dat,nlev = c(3,4,2), verbose=TRUE, adaptDF=TRUE)

suffStat <- list(dm = dat, nlev = c(3,4,2), adaptDF = FALSE)
disCItest(1,3,2,suffStat)
}


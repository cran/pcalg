\name{gSquareDis}
\alias{gSquareDis}
\title{Test for (conditional) independence for discrete data}
\description{
  G^2 statistic to test for (conditional) independence of X and Y given
  set S.
}
\usage{
gSquareDis(x, y, S, dm, nlev, verbose = FALSE, adaptDF = FALSE)
}
\arguments{
  \item{x}{Position of variable X in the adjacency matrix}
  \item{y}{Position of variable Y in the adjacency matrix}
  \item{S}{Position of the conditioning variables in the adjacency set}
  \item{dm}{Data matrix (rows: samples, columns: variables) with binary
    entries}
  \item{nlev}{Vector with numbers of levels for each variable}
  \item{verbose}{If TRUE, detailed output is provided.}
  \item{adaptDF}{Lower the degrees of freedom by one for each zero
    count. The value for the degrees of freedom cannot go below 1.}
}
\value{The p-value of the test.
}

\details{The G^2 statistic is used to test for (conditional)
  independence of X and Y given a set S (cann be \code{NULL}). If only
  binary variables are involved, \code{\link{gSquareBin}} is a
  specialized alternative to this function.}

\references{
  R.E. Neapolitan (2004).
  Learning Bayesian Networks. \emph{Prentice Hall Series in Artificial
  Intelligence}. Chapter 10.3.1
}
\seealso{\code{\link{gSquareBin}} for a (conditional) independence test
  for binary variables. \code{\link{disCItest}} for a wrapper of this
  function that can be easily included in \code{\link{skeleton}},
  \code{\link{pc}} or \code{\link{fci}}. 
}
\author{
  Nicoletta Andri and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\examples{
## Simulate data
x <- sample(1:3,100,TRUE)
y <- sample(1:4,100,TRUE)
z <- sample(1:2,100,TRUE)
dat <- cbind(x,y,z)

## Analyze data
gSquareDis(1,3,2,dat,nlev = c(3,4,2))
}


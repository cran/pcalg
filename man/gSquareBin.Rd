\name{gSquareBin}
\alias{gSquareBin}
\title{Test for (conditional) independence for binary data}
\description{
  G^2 statistic to test for (conditional) independence of binary
  variables X and Y given set of binary variables S.
}
\usage{
gSquareBin(x, y, S, dm, verbose = FALSE, adaptDF = FALSE)
}
\arguments{
  \item{x}{Position of variable X in the adjacency matrix}
  \item{y}{Position of variable Y in the adjacency matrix}
  \item{S}{Position of the conditioning variables in the adjacency set.}
  \item{dm}{Data matrix (rows: samples, columns: variables) with binary
    entries}
  \item{verbose}{If TRUE, detailed output is provided.}
  \item{adaptDF}{Lower the degrees of freedom by one for each zero
    count. The value for the degrees of freedom cannot go below 1.}
}
\value{The p-value of the test.
}

\details{ The G^2 statistic is used to test for (conditional)
  independence of X and Y given a set S (cann be \code{NULL}). This
  function is a specialized version of \code{\link{gSquareDis}} which is
  for discrete variables with more than two levels.}

\references{
  R.E. Neapolitan (2004).
  Learning Bayesian Networks. \emph{Prentice Hall Series in Artificial
  Intelligence}. Chapter 10.3.1
}
\seealso{\code{\link{gSquareDis}} for a (conditional) independence test
  for discrete variables with more than two
  levels. \code{\link{binCItest}} for a wrapper of this function that
  can be easily included in \code{\link{skeleton}}, \code{\link{pc}} or
  \code{\link{fci}}. 
}
\author{
  Nicoletta Andri and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\examples{
## Simulate data from a chain of 3 variables: x1 -> x2 -> x3
set.seed(123)
expit <- function(x) { exp(x) / (1 + exp(x)) }
b0 <- 0
b1 <- 1
b2 <- 1
n <- 10000
x1 <- sample(c(0,1),n,replace=TRUE)

p2 <- expit(b0 + b1*x1)
x2 <- rep(NA,length(x1))
for (i in 1:n) {
  x2[i] <- sample(c(0,1),1,prob=c(1-p2[i],p2[i]))
}
p3 <- expit(b0 + b2*x2)
x3 <- rep(NA,length(x2))
for (i in 1:n) {
  x3[i] <- sample(c(0,1),1,prob=c(1-p3[i],p3[i]))
}

dat <- cbind(x1,x2,x3)

## Test marginal and conditional independencies
gSquareBin(3,1,NULL,dat)
gSquareBin(3,1,2,dat)
gSquareBin(1,3,2,dat)
}


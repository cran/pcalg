\name{gaussCItest}
\alias{gaussCItest}
\title{Test for (conditional) independence for gaussian data}
\description{
  This function tests for (condtitional) independence between gaussian
  random variables. The function was organized in a way, so that it can
  easily be used in \code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}.
}
\usage{
gaussCItest(x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of node X in adjacency matrix}
  \item{y}{Position of node Y in adjacency matrix}
  \item{S}{Position of set of nodes in adjacency matrix on which is
    conditioned}
  \item{suffStat}{A list with two elements: (1) Element "C" containing
    the correlation matrix of the data and (2) element "n" containing
    the sample size.}
}
\value{The p-value of the test is returned.
}

\details{ The Fisher z transformation is used. See the help on
  \code{\link{zStat}} for details.}

\seealso{\code{\link{dsepTest}}, \code{\link{disCItest}} and
  \code{\link{binCItest}} for similar functions for a d-separation
  oracle, conditional independence test on discrete variables and
  conditional independence test on binary variables, respectively.}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
## simulate data: x -> y -> z
set.seed(29)
x <- rnorm(100)
y <- 3*x + rnorm(100)
z <- 2*y + rnorm(100)
dat <- cbind(x,y,z)

## analyze data
suffStat <- list(C = cor(dat), n = nrow(dat))
gaussCItest(1,3,NULL,suffStat) ## dependent
gaussCItest(1,3,2,suffStat) ## independent
}


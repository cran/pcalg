\name{disCItest}
\alias{disCItest}
\title{Test for (conditional) independence for discrete data}
\description{
  This function tests for (condtitional) independence between discrete
  random variables. The function was organized in a way, so that it can
  easily be used in \code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}.
}
\usage{
disCItest(x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of node X in adjacency matrix}
  \item{y}{Position of node Y in adjacency matrix}
  \item{S}{Position of set of nodes in adjacency matrix on which is
    conditioned}
  \item{suffStat}{A list with three elements: (1) Element "dm" containing
    the data matrix (columns are variables, rows are samples), (2)
    element "nlev" containing a vector with numbers of levels for each
    variable and (3) element "adaptDF" as a boolean variable indicating
    whether to lower the degrees of freedom by one for each zero
    count. (The value for the degrees of freedom cannot go below 1.)}
}
\value{The p-value of the test is returned.
}

\details{This function is based on \code{\link{gSquareDis}}; see its
  help file for details.}

\seealso{\code{\link{dsepTest}}, \code{\link{gaussCItest}} and
  \code{\link{binCItest}} for similar functions for a d-separation
  oracle, conditional independence test on gaussian variables and
  conditional independence test on binary variables, respectively.}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
## Simulate data
set.seed(123)
x <- sample(1:3,100,TRUE)
y <- sample(1:4,100,TRUE)
z <- sample(1:2,100,TRUE)
dat <- cbind(x,y,z)

suffStat <- list(dm = dat, nlev = c(3,4,2), adaptDF = FALSE)
disCItest(1,3,2,suffStat)
}


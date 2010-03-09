\name{binCItest}
\alias{binCItest}
\title{Test for (conditional) independence for binary data}
\description{
  This function tests for (condtitional) independence between binary
  random variables. The function was organized in a way, so that it can
  easily be used in \code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}.
}
\usage{
binCItest(x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of node X in adjacency matrix}
  \item{y}{Position of node Y in adjacency matrix}
  \item{S}{Position of set of nodes in adjacency matrix on which is
    conditioned}
  \item{suffStat}{A list with two elements: (1) Element "dm" containing
    the data matrix (columns are variables, rows are samples) and (3)
    element "adaptDF" as a boolean variable indicating 
    whether to lower the degrees of freedom by one for each zero
    count. (The value for the degrees of freedom cannot go below 1.)}
}
\value{The p-value of the test is returned.
}

\details{This function is based on \code{\link{gSquareBin}}; see its
  help file for details.}

\seealso{\code{\link{dsepTest}}, \code{\link{gaussCItest}} and
  \code{\link{disCItest}} for similar functions for a d-separation
  oracle, conditional independence test on gaussian variables and
  conditional independence test on discrete variables, respectively.}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
## Simulate data
set.seed(123)
x <- sample(1:2,100,TRUE)
y <- sample(1:2,100,TRUE)
z <- sample(1:2,100,TRUE)
dat <- cbind(x,y,z)

suffStat <- list(dm = dat, adaptDF = FALSE)
binCItest(1,3,2,suffStat)
}


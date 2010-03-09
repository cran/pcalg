\name{dsepTest}
\alias{dsepTest}
\title{Test for d-separation in a DAG}
\description{
  This function tests for d-separation of nodes in a DAG. The function
  was organized in a way, so that it can 
  easily be used in \code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}.
}
\usage{
dsepTest(x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of node X in adjacency matrix}
  \item{y}{Position of node Y in adjacency matrix}
  \item{S}{Position of set of nodes in adjacency matrix on which is
    conditioned}
  \item{suffStat}{A list with two elements: (1) Element "g" containing
    the Directed Acyclic Graph (object of \code{\link{class}}
    \code{"graph"}, see \code{\link[graph]{graph-class}} from the
    package \pkg{graph}) and (2) element "jp" containing the shortest
    path distance matrix for all pairs of nodes as computed by
    \code{johnson.all.pairs.sp} in \pkg{graph}.} 
}
\value{If x and y are d-separated by S in DAG G the result is 1, otherwise it is
  0. This is analogous to the p-value of an ideal (without sampling
  error) conditional independence test on any distribution that is
  faithful to the DAG G.
}

\details{The function is based on \code{\link{dsep}}. For details on
d-separation see the reference.}

\references{
  Kalisch M. and P. B\"uhlmann (2007)
  \emph{Estimating high-dimensional
    directed acyclic graphs with the PC-algorithm};
    JMLR, Vol. 8, 613-636, 2007.
}

\seealso{\code{\link{gaussCItest}}, \code{\link{disCItest}} and
  \code{\link{binCItest}} for similar functions for a conditional
  independence test on gaussian variables, on discrete variables and
  on binary variables, respectively.}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

## plot the DAG
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

## define sufficient statistics (d-separation oracle)
suffStat <- list(g = myDAG, jp = johnson.all.pairs.sp(myDAG))

dsepTest(1,6,NULL,suffStat) ## not d-separated
dsepTest(1,6,3,suffStat) ## not d-separated by node 3
dsepTest(1,6,c(3,4),suffStat) ## d-separated by node 3 and 4
}


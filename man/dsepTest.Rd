\name{dsepTest}
\alias{dsepTest}
\title{Test for d-separation in a DAG}
\description{
  This function tests for d-separation of nodes in a DAG. The function
  is written, so that it can
  easily be used in \code{\link{skeleton}}, \code{\link{pc}}, \code{\link{fci}}.
}
\usage{
dsepTest(x, y, S, suffStat)
}
\arguments{
  \item{x}{Position of variable X in the adjacency matrix}
  \item{y}{Position of variable Y in the adjacency matrix}
  \item{S}{Positions of conditioning variables in the adjacency matrix}
  \item{suffStat}{A list with two elements: (1) Element "g" containing
    the Directed Acyclic Graph (object of \code{\link{class}}
    \code{"graph"}, see \code{\link[graph]{graph-class}} from the
    package \pkg{graph}) and (2) element "jp" containing the shortest
    path distance matrix for all pairs of nodes as computed by
    \code{\link[RBGL]{johnson.all.pairs.sp}} from package \pkg{RBGL}.}
}
\value{If x and y are d-separated by S in DAG G the result is 1, otherwise it is
  0. This is analogous to the p-value of an ideal (without sampling
  error) conditional independence test on any distribution that is
  faithful to the DAG G.
}

\details{The function is based on \code{\link{dsep}}. For details on
d-separation see the Lauritzen (2004).}

\references{
  S.L. Lauritzen (2004),
  Graphical Models, \emph{Oxford University Press}.

}

\seealso{\code{\link{gaussCItest}}, \code{\link{disCItest}} and
  \code{\link{binCItest}} for similar functions for a conditional
  independence test for gaussian, discrete and
  binary variables, respectively.}

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

if (require(Rgraphviz)) {
## plot the DAG
plot(myDAG, main = "randomDAG(10, prob = 0.2)")
}

## define sufficient statistics (d-separation oracle)
suffStat <- list(g = myDAG, jp = RBGL::johnson.all.pairs.sp(myDAG))

dsepTest(1,6, S= NULL,  suffStat) ## not d-separated
dsepTest(1,6, S= 3,     suffStat) ## not d-separated by node 3
dsepTest(1,6, S= c(3,4),suffStat) ## d-separated by node 3 and 4
}


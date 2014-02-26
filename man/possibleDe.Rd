\name{possibleDe}
\alias{possibleDe}
\title{
Find possible descendants on definite status paths.
}
\description{
In a DAG, CPDAG, MAG or PAG determine which nodes are possible descendants
of x on definite status paths. 
}
\usage{
possibleDe(amat, x)
}
\arguments{
  \item{amat}{
    Adjacency matrix of the DAG, CPDAG, MAG or PAG.
}
  \item{x}{
    Node of interest.
}
}
\details{
A non-endpoint vertex \code{X} on a path \code{p} in a partial mixed
graph is said to be of a \emph{definite status} if it is either a collider or a
definite non-collider on \code{p}. The path \code{p} is said to be of a
\emph{definite status} if all non-endpoint vertices on the path are of a
definite status (see e.g. Maathuis and Colombo (2013), Def. 3.4).

A possible descendent of x can be reached moving to adjacent nodes of x
but never going against an arrowhead. 
}
\value{
  Vector with possible descendents.
}
\references{
M.H. Maathuis and D. Colombo. A generalized backdoor criterion. (\href{http://arxiv.org/abs/1307.5636}{arXiv:1307.5636v2}) 
}
\author{
Diego Colombo
}
\seealso{
\code{\link{backdoor}}
}
\examples{
amat <- matrix( c(0,3,0,0,0,0, 2,0,2,0,0,0, 0,3,0,0,0,0, 0,0,0,0,1,0,
0,0,0,1,0,1, 0,0,0,0,1,0), 6,6)
colnames(amat) <- rownames(amat) <- letters[1:6]
if(require(Rgraphviz)) {
plotAG(amat)
}

possibleDe(amat, 1) ## a, b are poss. desc. of a
possibleDe(amat, 4) ## d, e, f are poss. desc. of d
}
\keyword{misc}


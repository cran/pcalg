\name{possibleDe}
\alias{possibleDe}
\title{
[DEPRECATED] Find possible descendants on definite status paths.
}
\description{
 This function is DEPRECATED! Use \code{\link{possDe}} instead.

  In a DAG, CPDAG, MAG or PAG determine which nodes are possible descendants
of x on definite status paths. 
}
\usage{
possibleDe(amat, x)
}
\arguments{
  \item{amat}{adjacency matrix of type \link{amat.pag}}
  \item{x}{(integer) position of node \eqn{x} (node of interest) in the
    adjacency matrix.}
  }
\details{
A non-endpoint vertex \code{X} on a path \code{p} in a partial mixed
graph is said to be of a \emph{definite status} if it is either a collider or a
definite non-collider on \code{p}. The path \code{p} is said to be of a
\emph{definite status} if all non-endpoint vertices on the path are of a
definite status (see e.g. Maathuis and Colombo (2015), Def. 3.4).

A possible descendent of x can be reached moving to adjacent nodes of x
but never going against an arrowhead. 
}
\value{
  Vector with possible descendents.
}
\references{
  M.H. Maathuis and D. Colombo (2015). A generalized back-door
  criterion. \emph{Annals of Statistics} \bold{43} 1060-1088.
}
\author{
Diego Colombo
}
\seealso{
\code{\link{backdoor}}, \link{amatType}
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


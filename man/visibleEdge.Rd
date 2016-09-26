\name{visibleEdge}
\alias{visibleEdge}
\title{
Check visible edge.
}
\description{
Check if the directed edge from x to z in a MAG or in a PAG is visible or not.
}
\usage{
visibleEdge(amat, x, z)
}
\arguments{
  \item{amat}{Adjacency matrix of type \link{amat.pag}}
  \item{x,z}{(integer) position of variable \eqn{x} and \eqn{z},
    respectively, in the adjacency matrix.}
}
\details{
All directed edges in DAGs and CPDAGs are said to be visible. Given a
MAG M / PAG P, a directed edge A -> B in M / P is visible if there is a
vertex C not adjacent to B, such that there is an edge between C and A
that is into A, or there is a collider path between C and A that is into
A and every non-endpoint vertex on the path is a parent of B. Otherwise
A -> B is said to be invisible. (see Maathuis and Colombo (2015), Def. 3.1)
}
\value{
\code{TRUE} if edge is visible, otherwise \code{FALSE}.
}
\references{
  M.H. Maathuis and D. Colombo (2015). A generalized backdoor
  criterion. Annals of Statistics 43 1060-1088.
}
\author{
  Diego Colombo
}
\seealso{
\code{\link{backdoor}}
}
\examples{
amat <- matrix(c(0,3,0,0, 2,0,2,3, 0,2,0,3, 0,2,2,0), 4,4)
colnames(amat) <- rownames(amat) <- letters[1:4]
if(require(Rgraphviz)) {
plotAG(amat)
}

visibleEdge(amat, 3, 4) ## visible
visibleEdge(amat, 2, 4) ## visible
visibleEdge(amat, 1, 2) ## invisible
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

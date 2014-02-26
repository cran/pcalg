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
  \item{amat}{Adjacency matrix (coding 0,1,2,3 for no edge, circle,
    arrowhead, tail; e.g., \code{amat[a,b] = 2} and \code{amat[b,a] = 3}
    implies a -> b).
}
  \item{x}{Node x.
}
  \item{z}{Node z.
}
}
\details{
All directed edges in DAGs and CPDAGs are said to be visible. Given a
MAG M / PAG P, a directed edge A -> B in M / P is visible if there is a
vertex C not adjacent to B, such that there is an edge between C and A
that is into A, or there is a collider path between C and A that is into
A and every non-endpoint vertex on the path is a parent of B. Otherwise
A -> B is said to be invisible. (see Maathuis and Colombo (2013), Def. 3.1)
}
\value{
\code{TRUE} if edge is visible, otherwise \code{FALSE}.
}
\references{
  M.H. Maathuis and D. Colombo (2013). A generalized backdoor
  criterion. arXiv preprint arXiv:1307.5636.
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

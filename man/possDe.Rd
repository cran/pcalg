\name{possDe}
\alias{possDe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find possible descendants of given node(s).
}
\description{
In a DAG, CPDAG, MAG or PAG determine which nodes are (possible) descendants of x on definite status or just any paths potentially avoiding given nodes on the paths.
}
\usage{
possDe(m, x, y = NULL, possible = TRUE, ds = TRUE,
type = c("cpdag", "pdag", "dag", "mag", "pag"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
  Adjacency matrix in coding according to type.
  }
  \item{x}{
  Node positions of starting nodes.
}
  \item{y}{
  Node positions of nodes through which a path must not go.
}
  \item{possible}{
  If \code{TRUE}, possible descendents are returned.
}
  \item{ds}{
  If \code{TRUE}, only paths of definite status are considered.
}
  \item{type}{
  Type of adjacency matrix in \code{m}. The coding is according to \code{\link{amatType}}.
}
}
\details{
  Not all possible combinations of the arguments are currently implemented and will issue an error if called.
}
\value{
  Vector of all node positions found as (possible) descendents of the nodes in \code{x}.
}
\author{
  Markus Kalisch
}

\seealso{
  \code{\link{amatType}}
}
\examples{
## a -> b -- c
amat <- matrix(c(0,1,0, 0,0,1, 0,1,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
if (require(Rgraphviz)) plot(as(t(amat), "graphNEL"))

possDe(m = amat, x = 1, possible = TRUE, ds = FALSE, type = "pdag") ## all nodes
possDe(m = amat, x = 1, possible = FALSE, ds = FALSE, type = "pdag") ## only nodes 1 and 2
possDe(m = amat, x = 1, y = 2, possible = TRUE, ds = FALSE, type = "pdag") ## only node 1
}
\keyword{graphs}% use one of  RShowDoc("KEYWORDS")
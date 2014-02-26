\name{legal.path}
\alias{legal.path}
\title{Check if a 3-node-path is Legal}
\description{
  Check if the path \eqn{a--b--c}{a-b-c} is legal.

  A 3-node path \eqn{a--b--c}{a-b-c}  is \dQuote{legal} iff either
  \eqn{b} is a collider or \eqn{a--b--c}{a-b-c} is a triangle.
}
\usage{
legal.path(a, b, c, amat)
}
\arguments{
  \item{a,b,c}{Integer positions in adjacency matrix of nodes \eqn{a},
    \eqn{b}, and \eqn{c}, respectively.}
  \item{amat}{Adjacency matrix (coding 0,1,2,3 for no edge, circle,
  arrowhead, tail; e.g., \code{amat[a,b] = 2} and \code{amat[b,a] = 3}
  implies a -> b)} 
}
\note{Prerequisite: \eqn{a--b--c}{a-b-c} \emph{must} be in a path (and
  this is \emph{not} checked by \code{legal.path()}).
}
\value{\code{TRUE} if path is legal, otherwise \code{FALSE}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
amat <- matrix( c(0,1,1,0,0, 2,0,1,0,0, 2,2,0,2,1,
                  0,0,1,0,0, 0,0,2,0,0), 5,5)
legal.path(1,3,5, amat)
legal.path(1,2,3, amat)
legal.path(2,3,4, amat)
}
\keyword{misc}

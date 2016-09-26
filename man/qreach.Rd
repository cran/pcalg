\name{qreach}
\alias{qreach}
\title{Compute Possible-D-SEP(x,G) of a node x in a PDAG G}
\description{
  Let G be a graph with the following edge types: o-o, o-> or <->, and
  let x be a vertex in the graph.  Then this function computes
  Possible-D-SEP(x,G), which is defined as follows: 

  v is in Possible-D-SEP(x,G) iff there is a path p between x and v in G
  such that for every subpath <s,t,u> of p, t is a collider on this
  subpath or <s,t,u> is a triangle in G.

  See Spirtes et al (2000) or Definition 3.3 of Colombo et al (2012).
}
\usage{
  qreach(x, amat, verbose = FALSE)
}
\arguments{
  \item{x}{(integer) position of vertex \eqn{x} in the adjacency matrix of
  which Possible-D-SEP set is to be computed.}
  \item{amat}{Adjacency matrix of type \link{amat.pag}.
  }
  \item{verbose}{Logical, asking for details on output}
}
\value{
  Vector of column positions indicating the nodes in Possible-D-SEP of x.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  D. Colombo, M.H. Maathuis, M. Kalisch, T.S. Richardson (2012).
  Learning high-dimensional directed acyclic graphs with latent and
  selection variables. 
  \emph{Annals of Statistics} \bold{40}, 294--321.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{fci}} and \code{\link{pdsep}} which both use
  this function.
}




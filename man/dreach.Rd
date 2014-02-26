\name{dreach}
\alias{dreach}
\title{Compute D-SEP(x,y,G)}
\description{
  Let x and y be two distinct vertices in a mixed graph G.  This
  function computes D-SEP(x,y,G), which is defined as follows:

  A node v is in D-SEP(x,y,G) iff v is not equal to x and there is a collider
  path between x and v in G such that every vertex on this path is an
  ancestor of x or y in G.

  See p.136 of Sprirtes et al (2000) or Definition 4.1 of Maathuis and
  Colombo (2013). 
}
\usage{
  dreach(x, y, amat, verbose = FALSE)
}
\arguments{
  \item{x}{First argument of D-SEP, given as the column number of the
    node in the adjacency matrix.} 
  \item{y}{Second argument of D-SEP, given as the column number of the
    node in the adjacency matrix (\code{y} must be different from \code{x}).} 
  \item{amat}{Adjacency matrix (coding 0,1,2,3 for no edge, circle,
  arrowhead, tail; e.g., \code{amat[a,b] = 2} and \code{amat[b,a] = 3}
  implies a -> b) 
  }
  \item{verbose}{Logical specifying details should be on output}
}
\value{
  Vector of column positions indicating the nodes in D-SEP(x,y,G).
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  M.H. Maathuis and D. Colombo (2013). A generalized backdoor criterion. arXiv:1307.5636.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{backdoor}} uses this function;
  \code{\link{pag2magAM}}.% FIXME --> \code{\link{pag2mag}} once we have that
}
%% FIXME:
%% \examples{
%% }




\name{find.unsh.triple}
\alias{find.unsh.triple}
\title{Find all Unshielded Triples in an Undirected Graph}
\description{
  Find all unshielded triples in an undirected graph, \eqn{q}, i.e.,
  the ordered (\eqn{(x,y,z)} with \eqn{x < z}) list of all the triples
  in the graph.
}
\usage{
find.unsh.triple(g, check=TRUE)
}
\arguments{
  \item{g}{adjacency matrix of type \link{amat.cpdag} representing the
    skeleton; since a skeleton consists only of undirected edges, 
    \code{g} must be symmetric.}
  \item{check}{logical indicating that the symmetry of \code{g} should
    be checked.}
}
\details{
  A triple of nodes \code{x}, \code{y} and \code{z} is
  \dQuote{unshielded}, if (all of these are true): \describe{
    \item{(i)}{\code{x} and \code{y} are connected;}
    \item{(ii)}{\code{y} and \code{z} are connected;}
    \item{(iii)}{\code{x} and \code{z} are \bold{not} connected.}
  }
}
\value{
  \item{unshTripl}{Matrix with 3 rows containing in each column
                   an unshielded triple}
  \item{unshVect}{Vector containing the unique number  for
    each column in unshTripl (for internal use only)}
}
\author{
 Diego Colombo, Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), and
 Martin Maechler
}
\examples{
data(gmG)
if (require(Rgraphviz)) {
  ## show graph
  plot(gmG$g, main = "True DAG")
}

## prepare skeleton use in example
g <- wgtMatrix(gmG$g) ## compute weight matrix
g <- 1*(g != 0) # wgts --> 0/1; still lower triangular
print.table(g, zero.print=".")
skel <- g + t(g) ## adjacency matrix of skeleton

## estimate unshielded triples -- there are  13 :
(uTr <- find.unsh.triple(skel))
\dontshow{stopifnot(ncol(uTr$unshTripl) == 13)}
}
\keyword{misc}

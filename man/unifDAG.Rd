\name{unifDAG}
\title{Uniform Sampling of Directed Acyclic Graphs (DAG)}
\alias{unifDAG}
\alias{unifDAG.approx}
\concept{exact enumeration}
\description{
  Uniform sampling of a labelled directed acyclic graph (DAG) with
  combinatorial enumeration.
}
\usage{
unifDAG (n, weighted=FALSE, wFUN=list(runif, min=0.1, max=1))
unifDAG.approx(n, n.exact=20, weighted=FALSE, wFUN=list(runif,min=0.1,max=1))
}
\arguments{
  \item{n}{integer larger than \code{1}, indicating the number of
    nodes in the DAG.  \code{unifDAG} can only be used for \code{n} up to
    \code{100}.  For larger \code{n}, use \code{unifDAG.approx}.}
  \item{weighted}{logical indicating if weights of the edges are
    computed according to \code{wFUN}.}
  \item{wFUN}{a \code{\link{function}} for computing the weights of the
    edges in the DAG.  It takes as first argument a number of edges
    \code{m} for which it returns a vector of length \code{m} containing
    the weights.  Alternatively, it could be a list consisting of the
    function in the first entry and of further arguments of the function
    in the additional entries.  The default (only if \code{weighted} is
    true)) is a uniform weight between \code{0.1} and \code{1}.  See the
    examples.}
  \item{n.exact }{an integer, at least \code{n} and between \code{2} and
    \code{100}, denoting the number of nodes up to which the exact
    method is used, followed by an approximation for larger numbers up
    to \code{n}.  See details on the quality of the approximation.}
}

\details{
  A (weighted) random graph with \code{n} nodes is uniformly drawn over
  the space of all labelled DAGs with \code{n} nodes.
  The main idea of these two methods is to first sample a random
  sequence of outpoints, that is, nodes without incoming edges.  This
  sequence is then used to construct an adjacency matrix, which is
  converted to the final DAG.  The presented methods differ only in the
  approach to find this sequence of outpoints.

  The method \code{unifDAG} builds the random sequence of outpoints
  based on precomputed enumeration tables.

  The method \code{unifDAG.approx} executes \code{unifDAG} up to the
  number \code{n.exact}, for larger number of nodes an approximation is
  used instead.  The default of \code{n.exact = 20} (\code{40}) should
  get the approximation within the uniformity limits of a 32 (64) bit
  integer sampler.  See reference for more details.
}

\value{A graph object of class \code{\link[graph]{graphNEL}}.}

\references{
  Jack Kuipers and Guisi Moffa (2015)
  Uniform random generation of large acyclic digraphs.
  \emph{Statistics and Computing} \bold{25(2)}, 227--242, Springer;
  \url{http://dx.doi.org/10.1007/s11222-013-9428-y}
}

\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Manuel Schuerch.}

\note{
  The main advantage of these algorithms is that they operate on the space
  of DAGs instead of the space of undirected graphs with an additional
  phase of orienting the edges.  With this approach the unintentional bias
  towards sparse graphs, for instance occurring when sampling adjacency
  matrices, can be eliminated.
}

\seealso{\code{\link{randDAG}} for generating different random DAGs,
  each intentionally shifted towards certain properties;
  \code{\link{randomDAG}} a limited and soon deprecated version of
  \code{\link{randDAG}} where the output is topologically sorted;
  \code{\link{rmvDAG}} for generating multivariate data according to a DAG.
}

\examples{
# require("Rgraphviz")
set.seed(123)
dag1 <- unifDAG(n=10)
dag2 <- unifDAG.approx(n=10, n.exact=5)

dag <- unifDAG(n=5)
plot(dag)
dag@edgeData   ## note the constant weights

dag <- unifDAG(n=5,weighted=TRUE)
plot(dag)
dag@edgeData   ## note the uniform weights between 0.1 and 1

wFUN <- function(m,lB,uB) { runif(m,lB,uB) }
dag <- unifDAG(n=5,weighted=TRUE,wFUN=list(wFUN,1,4))
dag@edgeData   ## note the uniform weights between 1 and 4
}
\keyword{graphs}
\keyword{models}

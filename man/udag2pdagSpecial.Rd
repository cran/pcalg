\name{udag2pdagSpecial}
\alias{udag2pdagSpecial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a PDAG using
  different methods if problems occur}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to a Partially
Directed Acyclic Graph (PDAG). The result is a pcAlgo-object 
as well.
}
\usage{
udag2pdagSpecial(gInput,verbose,n.max=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{pcAlgo-object containing skeleton and cond. ind. information}
  \item{verbose}{0: No output; 1: Details}
  \item{n.max}{Maximum number of tries for reontienting doubly visited edges.}
}
\details{
  The skeleton is extended to a PDAG using rules by Spirtes and Pearl
  (see References).
  If, after orienting the v-structures, the graph is not extendible to a
  DAG, the following mechanisms try to solve the problem:
  There might be edges, that were oriented twice while forming
  colliders. In this case, try all possible combinations of
  reorientations and check whether any reorientation is extendible. Take
  the first one that is extendible and make no more than n.max tries.
  If this fails, the original graph is extended arbitrarily to a DAG
  that fits on the skeleton. v-structures might have changed. The
  resulting DAG is then transformed to its CPDAG.
}
\value{
  \item{pcObj}{Oriented pc-Object}
  \item{evisit}{Matrix counting the number of orientation attemps per edge}
  \item{xtbl.orig}{Is original graph with v-structure extendable}
  \item{xtbl}{Is final graph with v-structure extendable}
  \item{amat0}{Adj.matrix of original graph with v-structures}
  \item{amat1}{Adj.matrix of graph with v-structures after reorienting
    edges from double edge visits}
  \item{status}{0: original try is extendable; 1: reorienting double edge
    visits helps; 2: orig. try is not extendable; reorienting double
    visits don't help; result is acyclic, has orig. v-structures, but perhaps
    additional v-structures}
  \item{counter}{Number of reorientation tries until success or max.tries}
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  J. Pearl (2000), \emph{Causality}, Cambridge University Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pdag2dag}}, \code{\link{dag2cpdag}},
  \code{\link{udag2pdag}}, \code{\link{udag2pdagRelaxed}}}
\examples{
## simulate data
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2)
trueCPDAG <- dag2cpdag(myDAG)
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton
resU <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard",directed=FALSE)
## orient edges using three different methods
resD1 <- udag2pdagRelaxed(resU,verbose=0)
resD2 <- udag2pdagSpecial(resU,verbose=0,n.max=100)
resD3 <- udag2pdag(resU,verbose=0)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



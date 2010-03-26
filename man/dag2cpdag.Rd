\name{dag2cpdag}
\alias{dag2cpdag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert a DAG to a CPDAG }
\description{
Convert a DAG to a Completed Partially Directed Acyclic Graph (CPDAG).
}
\usage{
dag2cpdag(dag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dag}{DAG (graph object)}
}
\details{
  This function converts a DAG (graph object) to its corresponding
  (unique) CPDAG (graph object), using the algorithm of
  Chickering (2002).
}
\value{
  A graph object containing the CPDAG.
}
\references{
D.M. Chickering (2002), Learning Equivalence Classes of Bayesian-Network
  Structures, \emph{Journal of Machine Learning Research} 2, 445-398.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{udag2pdag}}, \code{\link{pdag2dag}}}
\examples{
p <- 10 # number of random variables
s <- 0.4 # sparseness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG

res <- dag2cpdag(g)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


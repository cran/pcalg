\name{pdag2dag}
\alias{pdag2dag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a PDAG to a DAG}
\description{
This function extends a PDAG to a DAG, if this is possible.
}
\usage{
pdag2dag(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{input PDGA (graph object)}
}
\details{
  The PDAG is consistently extended using an algorithm by Dor and Tarsi
  (see References). If no extension is possible, a DAG corresponding to
  the skeleton of the PDAG is generated and a warning message is produced.
}
\value{
  consistent DAG extension (graph object)
}
\references{
  D.Dor, M.Tarsi, \emph{A simple algorithm to construct a
    consistent extension of a partially oriented graph, Technicial
    Report R-185, Cognitive Systems Laboratory, UCLA, 1992} 
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{udag2pdag}},\code{\link{pdag2dag}},\code{\link{udag2cpdag}}}
\examples{
p <- 10 # number of random variables
n <- 10000 # number of samples
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
d <- rmvDAG(n,g) # generate random samples

gSkel <- 
  pcAlgo(d,alpha=0.05) # estimate of the skeleton

gPDAG <- udag2pdag(gSkel)

gDAG <- pdag2dag(gPDAG@graph)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

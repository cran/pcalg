\name{pdag2dag}
\alias{pdag2dag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a Partially Directed Acyclic Graph (PDAG) to a DAG}
\description{
This function extends a PDAG to a DAG, if this is possible.
}
\usage{
pdag2dag(g,keepVstruct=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{input PDAG (graph object)}
  \item{keepVstruct}{If TRUE, the v-structures in g are kept. Otherwise they
    are ignored and an arbitrary extension is generated.}
}
\details{
  Direct undirected edges without creating directed cycles or additional
  v-structures. 
  The PDAG is consistently extended to a DAG using the algorithm by Dor and Tarsi
  (1992). If no extension is possible, a DAG corresponding to
  the skeleton of the PDAG is generated and a warning message is produced.
}
\value{
  List with two entries: "graph" contains a consistent DAG extension
  (graph object) and "success" is \code{TRUE} iff the extension was possible.
}
\references{
  D.Dor, M.Tarsi (1992). A simple algorithm to construct a
    consistent extension of a partially oriented graph. Technicial
    Report R-185, Cognitive Systems Laboratory, UCLA
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\examples{
p <- 10 # number of random variables
n <- 10000 # number of samples
s <- 0.4 # sparsness of the graph

## generate random data
set.seed(42)
g <- randomDAG(p,s) # generate a random DAG
d <- rmvDAG(n,g) # generate random samples

gSkel <- pcAlgo(d,alpha=0.05) # estimate of the skeleton

gPDAG <- udag2pdag(gSkel)

gDAG <- pdag2dag(gPDAG@graph)
}

\keyword{graphs}

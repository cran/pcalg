\name{udag2cpdag}
\alias{udag2cpdag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a CPDAG}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to a Completed
Partially Directed Acyclic Graph (CPDAG). The result is a pcAlgo-object
as well.
}
\usage{
udag2cpdag(pc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc}{pcAlgo-object containing skeleton and cond. ind. information}
}
\details{
  The skeleton is first extended to a PDAG using the function
  \code{udag2pdag}. Then, it is extended to a DAG using the functino
  \code{pdag2dag} and is then transferred to a CPDAG using the function
  \code{dag2cpdag}. 
}
\value{
  A pcAlgo-object containing the CPDAG.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{udag2pdag}},\code{\link{pdag2dag}},\code{\link{dag2cpdag}}}
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
gCPDAG <- 
  udag2cpdag(gSkel) # transform skeleton to cpdag

plot(g)
plot(gSkel,zvalue.lwd=TRUE)
plot(gCPDAG,zvalue.lwd=TRUE)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}


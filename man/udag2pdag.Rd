\name{udag2pdag}
\alias{udag2pdag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a PDAG}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to a Partially
Directed Acyclic Graph (PDAG). The result is a pcAlgo-object 
as well.

}
\usage{
udag2pdag(gInput,verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{pcAlgo-object containing skeleton and
    cond. ind. information}
  \item{verbose}{0: No output; 1: Details}
}
\details{
  The skeleton is extended to a PDAG using rules by Spirtes and Pearl
  (see References). 
}
\value{
  A pcAlgo-object containing the PDAG.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  J. Pearl (2000), \emph{Causality}, Cambridge University Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pdag2dag}},\code{\link{dag2cpdag}}}
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
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



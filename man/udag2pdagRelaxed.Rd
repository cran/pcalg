\name{udag2pdagRelaxed}
\alias{udag2pdagRelaxed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a PDAG}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to a Partially
Directed Acyclic Graph (PDAG). The result is a pcAlgo-object 
as well. There is no check whether the result is extendible to a DAG
}
\usage{
udag2pdagRelaxed(gInput, verbose, unfVect=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{pcAlgo-object containing skeleton and
    cond. ind. information}
  \item{verbose}{0: No output; 1: Details}
  \item{unfVect}{Vector containing numbers that encode the unfaithful
    triple (as returned by pc.cons.intern). This is needed in the conservative PC.}
}
\details{
  The skeleton is extended to a PDAG using rules by Spirtes and Pearl
  (see References). There is no test whether the result is really
  extendible. 
}
\value{
  \item{pcObj}{Oriented pc-Object}
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  J. Pearl (2000), \emph{Causality}, Cambridge University Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pdag2dag}}, \code{\link{dag2cpdag}},
  \code{\link{udag2pdag}}, \code{\link{udag2pdagRelaxed}}, \code{\link{udag2pdagSpecial}}}
\examples{
## simulate data
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2)
trueCPDAG <- dag2cpdag(myDAG)
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## define independence test (partial correlations)
indepTest <- gaussCItest 
## define sufficient statistics
suffStat <- list(C = cor(d.mat), n = n)
## estimate skeleton
resU <- skeleton(suffStat, indepTest, p, 0.05)

## orient edges using three different methods
resD1 <- udag2pdagRelaxed(resU,verbose=0)
resD2 <- udag2pdagSpecial(resU,verbose=0,n.max=100)
resD3 <- udag2pdag(resU,verbose=0)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



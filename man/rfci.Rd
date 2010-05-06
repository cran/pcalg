\name{rfci}
\alias{rfci}
\title{Estimate the equivalence class of a MAG (PAG) using the RFCI Algorithm}
\description{
Estimate the RFCI-PAG from observational data, using the RFCI-algorithm.
}

\usage{
rfci(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf)
}

\arguments{
  \item{suffStat}{Sufficient statistics: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{fixedGaps}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    removed before starting the algorithm. Therefore, this edge is
    guaranteed to be absent in the resulting graph.}
  \item{fixedEdges}{A logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    never considered for removal. Therefore, this edge is
    guaranteed to be present in the resulting graph.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{m.max}{Maximal size of the conditioning sets that are considered in the
    conditional independence tests.}
}
\value{An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}

\details{ This function is rather similar to \link{fci}. However, it
does not compute any possible d-sep sets and thus does not make tests
conditioning on them. This makes RFCI much faster than FCI. The
orientation rules for v-structures and rule 4 were modified in order to
produce a PAG which, in the oracle version, is guaranteed to have the
correct ancestral relationships.}

\references{
  Diego Colombo, Marloes H. Maathuis, Markus Kalisch, Thomas
  S. Richardson (2011). Learning high-dimensional DAGs with latent and
  selection variables. \emph{Arxiv}, http://arxiv.org/abs/1104.5617
}
\seealso{\code{\link{fci}} for estimating a PAG using the FCI algorithm;
  \code{\link{skeleton}} for estimating a skeleton
  using the PC algorithm; \code{\link{pc}} for estimating a CPDAG using
  the PC algorithm; \code{\link{gaussCItest}},
  \code{\link{disCItest}}, \code{\link{binCItest}} and
  \code{\link{dsepTest}} as examples for \code{indepTest}.
}
\author{
  Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch}).
}
\examples{
##################################################
## Example without latent variables
##################################################

set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)
amat <- wgtMatrix(myDAG)
amat[amat!=0] <- 1
amat <- amat + t(amat)
amat[amat!=0] <- 1

myCPDAG <- dag2cpdag(myDAG)
amat2 <- t(wgtMatrix(myCPDAG))

mycor <- cov2cor(trueCov(myDAG))

## find skeleton and CPDAG using the RFCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
indepTest <- gaussCItest
res <- rfci(suffStat, indepTest, p, alpha = 0.99, verbose=TRUE)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

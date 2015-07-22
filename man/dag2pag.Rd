\name{dag2pag}
\title{Convert a DAG with latent variables into a PAG}
\alias{dag2pag}
%%\alias{skeleton.dag2pag}
\encoding{UTF8}
\description{
  Convert a DAG with latent variables into its corresponding (unique)
  Partial Ancestral Graph (PAG).
}
\usage{
dag2pag(suffStat, indepTest, graph, L, alpha, rules = rep(TRUE,10),
        verbose = FALSE)
}
\arguments{
  \item{suffStat}{the sufficient statistics, a \code{\link{list}}
    containing all necessary elements for the conditional independence
    decisions in the function \code{indepTest}.}
  \item{indepTest}{a \code{\link{function}} for testing conditional
    independence.  The function is internally called as
    \code{indepTest(x,y,S,suffStat)}, and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and
    \code{y} are variables, and \code{S} is a (possibly empty) vector of
    variables (all variables are denoted by their column numbers
    in the adjacency matrix).  \code{suffStat} is a list containing
    all relevant elements for the conditional independence
    decisions.  The return value of \code{indepTest()} is the p-value of
    the test for conditional independence.}
  \item{graph}{a DAG with \code{p} nodes, a \pkg{graph} object.  The
    graph must be topological sorted (for example produced using
    \code{\link{randomDAG}}).}
  \item{L}{array containing the labels of the nodes in the \code{graph}
    corresponding to the latent variables.}
  \item{alpha}{significance level in \eqn{(0,1)} for the individual
    conditional independence tests.}
  \item{rules}{logical vector of length 10 indicating which rules
    should be used when directing edges.  The order of the rules is
    taken from Zhang (2009).}
  \item{verbose}{logical; if \code{TRUE}, detailed output is provided.}
}
\details{
  This function converts a DAG (graph object) with latent variables into
  its corresponding (unique) PAG, an \code{\linkS4class{fciAlgo}} class
  object, using the ancestor information and conditional independence
  tests entailed in the true DAG.  The output of this function is
  exactly the same as the one using
\preformatted{fci(suffStat, gaussCItest, p, alpha, rules = rep(TRUE, 10))
}
  using the true correlation matrix in \code{gaussCItest()} with a large
  \dQuote{virtual sample size} and a large alpha, but it is much faster,
  see the example.
}
\value{An object of \code{\link{class}} \code{\linkS4class{fciAlgo}},
  containing the estimated graph (in the form of an adjacency matrix
  with various possible edge marks), the conditioning sets that lead to
  edge removals (sepset) and several other parameters.
}
\references{
  Richardson, T. and Spirtes, P. (2002).
  Ancestral graph Markov models.
  \emph{Ann. Statist.} \bold{30}, 962--1030; Theorem 4.2., page 983.
}
\author{Diego Colombo and
  Markus Kalisch \email{kalisch@stat.math.ethz.ch}.
}
\seealso{
  \code{\link{fci}}, \code{\link{pc}}
}
\examples{
## create the graph
set.seed(78)
g <- randomDAG(10, prob = 0.25)
graph::nodes(g) # "1" "2" ... "10" % FIXME: should be kept in result!

## define nodes 2 and 6 to be latent variables
L <- c(2,6)

## compute the true covariance matrix of g
cov.mat <- trueCov(g)
## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)

## Find PAG
## as dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
system.time(
true.pag <- dag2pag(suffStat = list(C = true.corr, n = 10^9),
                    indepTest = gaussCItest,
                    graph=g, L=L, alpha = 0.9999) )

### ---- Find PAG using fci-function --------------------------

## From trueCov(g), delete rows and columns belonging to latent variable L
true.cov1 <- cov.mat[-L,-L]
## transform covariance matrix into a correlation matrix
true.corr1 <- cov2cor(true.cov1)

## Find PAG with FCI algorithm
## as dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
system.time(
true.pag1 <- fci(suffStat = list(C = true.corr1, n = 10^9),
                 indepTest = gaussCItest,
                 p = ncol(true.corr1), alpha = 0.9999) )

## confirm that the outputs are equal
stopifnot(true.pag@amat == true.pag1@amat)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

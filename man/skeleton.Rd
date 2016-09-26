\name{skeleton}
\alias{skeleton}
\title{Estimate (Initial) Skeleton of a DAG using the PC / PC-Stable Algorithm}
\description{
  Estimate the skeleton of a DAG without latent and selection
  variables using the PC Algorithm or
  estimate an initial skeleton of a DAG with arbitrarily many latent and
  selection variables using the FCI and the RFCI algorithms.

  If used in the PC algorithm, it estimates the order-independent
  \dQuote{PC-stable} (\code{"stable"}) or original PC (\code{"original"})
  \dQuote{skeleton} of a directed acyclic graph (DAG) from observational
  data.

  When used in the FCI and RFCI algorithms, this function estimates
  only an initial order-independent (or PC original) \dQuote{skeleton}.
  Because of the presence of latent and selection variables, to find the
  final skeleton those algorithms need to perform additional tests later on and
  consequently some edges can be further deleted.
}
\usage{
skeleton(suffStat, indepTest, alpha, labels, p,
         method = c("stable", "original", "stable.fast"), m.max = Inf,
         fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
         numCores = 1, verbose = FALSE)
}
\arguments{
  \item{suffStat}{Sufficient statistics: List containing all necessary
    elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Predefined \code{\link{function}} for testing
    conditional independence.  The function is internally called as
    \code{indepTest(x,y,S,suffStat)} and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and
    \code{y} are variables, and \code{S} is a (possibly empty) vector of
    variables (all variables are denoted by their column numbers
    in the adjacency matrix). \code{suffStat} is a list containing
    all relevant elements for the conditional independence decisions.
    The return value of \code{indepTest} is the p-value of the test for
    conditional independence.}
  \item{alpha}{significance level (number in \eqn{(0,1)} for the
    individual conditional independence tests.}
  \item{labels}{(optional) character vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.}
  \item{method}{Character string specifying method; the default,
    \code{"stable"} provides an \emph{order-independent} skeleton, see
    \sQuote{Details} below.}
  \item{m.max}{Maximal size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{fixedGaps}{logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i,j]} is true, the edge \eqn{i--j}{i-j} is removed before
    starting the algorithm.  Therefore, this edge is guaranteed to be
    \emph{absent} in the resulting graph.}
  \item{fixedEdges}{a logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i,j]} is true, the edge \eqn{i--j}{i-j} is never considered
    for removal.  Therefore, this edge is guaranteed to be \emph{present} in
    the resulting graph.}
  \item{NAdelete}{logical needed for the case \code{indepTest(*)}
    returns \code{NA}.  If it is true, the corresponding edge is
    deleted, otherwise not.}
  \item{numCores}{number of processor cores to use for parallel computation.
    Only available for \code{method = "stable.fast"}.}
  \item{verbose}{if \code{TRUE}, detailed output is provided.}
}
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the skeleton of
  the underlying DAG, the conditioning sets (\code{sepset}) that led to
  edge removals and several other parameters.
}
\details{
  Under the assumption that the distribution of the observed variables
  is faithful to a DAG and that there are \bold{no} latent and selection
  variables, this function estimates the skeleton of the DAG.  The
  skeleton of a DAG is the undirected graph resulting from removing all
  arrowheads from the DAG.  Edges in the skeleton of a DAG have the
  following interpretation:
  \cr
  There is an edge between \eqn{i} and \eqn{j}, \eqn{i} -- \eqn{j},
  if and only if variables \eqn{i} and \eqn{j} are conditionally
  dependent given \eqn{S} for all possible subsets \eqn{S} of the
  remaining nodes.

  On the other hand, the distribution of the observed variables
  is faithful to a DAG with \bold{arbitrarily many} latent and selection
  variables, \code{skeleton()} estimates the initial skeleton of the
  DAG.  Edges in this initial skeleton of a DAG have the
  following interpretation:
  \cr
  There is an edge \eqn{i} -- \eqn{j} if and only if variables \eqn{i} and
  \eqn{j} are conditionally dependent given \eqn{S} for all possible
  subsets \eqn{S} of the neighbours of \eqn{i} and the neighbours of
  \eqn{j}.

  The data are not required to follow a specific distribution,
  but one should make sure that the conditional indepedence test used in
  \code{indepTest} is appropriate for the data.  Pre-programmed versions
  of \code{indepTest} are available for Gaussian data
  (\code{\link{gaussCItest}}), discrete data (\code{\link{disCItest}}),
  and binary data (see \code{\link{binCItest}}).  Users may also specify
  their own \code{indepTest} function.

  The PC algorithm (Spirtes, Glymour and Scheines, 2000)
  (\code{method = "original"}) is known to be order-dependent, in the
  sense that the output may depend on the order in which the variables
  are given.  Therefore, Colombo and Maathuis (2014) proposed a simple
  modification, called \dQuote{PC-stable}, which yields
  order-independent adjacencies in the skeleton, provided by \code{pc()}
  with the new default \code{method = "stable"}.  This stable variant
  of the algorithm is also available with the \code{method = "stable.fast"}:
  it runs the algorithm of Colombo and Maathuis (2014) faster than
  \code{method = "stable"} in general, but should be regarded as an
  experimental option at the moment.

  The algorithm starts with a complete undirected graph.  In each
  step, it visits all pairs \eqn{(i,j)} of adjacent nodes in the
  current graph, and determines based on conditional independence tests
  whether the edge \eqn{i--j}{i-j} should be removed.  In particular, for each step
  \eqn{m} (\eqn{m=0,1,\dots}) of the size of the conditioning sets, the
  algorithm at first determines the neighbours \eqn{a(i)} of each node
  \eqn{i} in the graph.  Then, the algorithm visits all pairs \eqn{(i,j)}
  of adjacent nodes in the current graph, and the edge \eqn{i--j}{i-j} is
  kept if and only if the null hypothesis
  \cr
  \emph{\eqn{i} and \eqn{j} are conditionally independent given S}
  \cr
  rejected at significance level \code{alpha} for all subsets \eqn{S} of size
  \eqn{m} of \eqn{a(i)} and of \eqn{a(j)} (as judged by the function
  \code{indepTest}).  For the \code{"stable"} method, the neighborhoods
  \eqn{a(i)} are kept fixed within each value of \eqn{m}, and this
  makes the algorithm order-independent.  Method \code{"original"},
  the original PC algorithm would update the neighbour list after each
  edge change.

  The algorithm stops when \eqn{m} is larger than the largest
  neighbourhood size of all nodes, or when \eqn{m} has reached the limit
  \code{m.max} which may be set by the user.

  Since the FCI (Spirtes, Glymour and Scheines, 2000) and RFCI (Colombo
  et al., 2012) algorithms are built up from the PC algorithm, they are also
  order-dependent in the skeleton.  To resolve their order-dependence
  issues in the skeleton is more involved, see Colombo and Maathuis
  (2014).  However now, with \code{method = "stable"}, this function
  estimates an initial order-independent skeleton in these algorithms
  (for additional details on how to make the final skeleton of FCI fully
  order-independent see \code{\link{fci}} and Colombo and Maathuis (2014)).

  The information in \code{fixedGaps} and \code{fixedEdges} is used as follows.
  The gaps given in \code{fixedGaps} are introduced in the very beginning of
  the algorithm by removing the corresponding edges from the complete
  undirected graph.  Pairs \eqn{(i,j)} in \code{fixedEdges} are skipped
  in all steps of the algorithm, so that these edges remain in the graph.

  %% MM: why say this here?
  %% FIXME that it even does not *return* the correct variable names !
  Note: Throughout, the algorithm works with the column positions of
  the variables in the adjacency matrix, and not with the names of
  the variables.
}
\references{
  D. Colombo and M.H. Maathuis (2014).Order-independent constraint-based
  causal structure learning. \emph{Journal of Machine Learning Research}
  \bold{15} 3741-3782. 

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294-321.

  M. Kalisch and P. Buehlmann (2007).
  \emph{Estimating high-dimensional directed acyclic graphs with the
  PC-algorithm}, JMLR \bold{8} 613-636.

  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition, MIT Press.
}

\seealso{
  \code{\link{pc}} for generating a partially directed graph
  using the PC algorithm; \code{\link{fci}} for generating a partial
  ancestral graph using the FCI algorithm; \code{\link{rfci}} for
  generating a partial ancestral graph using the RFCI algorithm;
  \code{\link{udag2pdag}} for converting the skeleton to a CPDAG.

  Further, \code{\link{gaussCItest}}, \code{\link{disCItest}},
  \code{\link{binCItest}} and \code{\link{dsepTest}} as examples for
  \code{indepTest}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), Martin Maechler,
  Alain Hauser, and Diego Colombo.
}
\examples{
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x) # labels aka node names

## estimate Skeleton
skel.fit <- skeleton(suffStat = list(C = cor(gmG8$x), n = n),
                     indepTest = gaussCItest, ## (partial correlations)
                     alpha = 0.01, labels = V, verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow=c(1,2))
  plot(skel.fit, main = "Estimated Skeleton")
  plot(gmG8$g, main = "True DAG")
}

##################################################
## Using d-separation oracle
##################################################

## define sufficient statistics (d-separation oracle)
Ora.stat <- list(g = gmG8$g, jp = RBGL::johnson.all.pairs.sp(gmG8$g))
## estimate Skeleton
fit.Ora <- skeleton(suffStat=Ora.stat, indepTest = dsepTest, labels = V,
                    alpha=0.01) # <- irrelevant as dsepTest returns either 0 or 1

if (require(Rgraphviz)) {
  ## show estimated Skeleton
  plot(fit.Ora, main = "Estimated Skeleton (d-sep oracle)")
  plot(gmG8$g, main = "True DAG")
}
##################################################
## Using discrete data
##################################################
## Load data
data(gmD)
V <- colnames(gmD$x) # labels aka node names

## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)

## estimate Skeleton
skel.fit <- skeleton(suffStat,
                     indepTest = disCItest, ## (G^2 statistics independence test)
                     alpha = 0.01, labels = V, verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow = c(1,2))
  plot(skel.fit, main = "Estimated Skeleton")
  plot(gmD$g, main = "True DAG")
}

##################################################
## Using binary data
##################################################
## Load binary data
data(gmB)
X <- gmB$x

## estimate Skeleton
skel.fm2 <- skeleton(suffStat = list(dm = X, adaptDF = FALSE),
                     indepTest = binCItest, alpha = 0.01,
                     labels = colnames(X), verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow = c(1,2))
  plot(skel.fm2, main = "Binary Data 'gmB': Estimated Skeleton")
  plot(gmB$g, main = "True DAG")
}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

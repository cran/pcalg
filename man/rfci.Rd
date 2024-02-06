\name{rfci}
\alias{rfci}
\title{Estimate an RFCI-PAG using the RFCI Algorithm}
\description{
  Estimate an RFCI-PAG from observational data, using the RFCI-algorithm.
}
\usage{
rfci(suffStat, indepTest, alpha, labels, p,
     skel.method = c("stable", "original", "stable.fast"),
     fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
     m.max = Inf, rules = rep(TRUE, 10),
     conservative = FALSE, maj.rule = FALSE, 
     numCores = 1, verbose = FALSE)
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
  \item{alpha}{significance level (number in \eqn{(0,1)} for the
    individual conditional independence tests.}
  \item{labels}{(optional) character vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.}
  \item{skel.method}{Character string specifying method; the default,
    \code{"stable"} provides an \emph{order-independent} skeleton, see
    \code{\link{skeleton}}.}
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
  \item{m.max}{Maximum size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{rules}{Logical vector of length 10 indicating which rules
    should be used when directing edges. The order of the rules is taken
    from Zhang (2009).}
  \item{conservative}{Logical indicating if the unshielded triples
    should be checked for ambiguity after the skeleton has been found,
    similar to the conservative PC algorithm.}
  \item{maj.rule}{Logical indicating if the unshielded triples should be
    checked for ambiguity after the skeleton has been found using a
    majority rule idea, which is less strict than the conservative.}
  \item{numCores}{Specifies the number of cores to be used for parallel
    estimation of \code{\link{skeleton}}.}
  \item{verbose}{If true, more detailed output is provided.}
}
\value{An object of \code{\link{class}} \code{fciAlgo} (see
 \code{\linkS4class{fciAlgo}}) containing the estimated graph
 (in the form of an adjacency matrix with various possible edge marks),
 the conditioning sets that lead to edge removals (sepset) and several other
 parameters.
}

\details{
  This function is rather similar to \link{fci}. However, it does not
  compute any Possible-D-SEP sets and thus does not make tests
  conditioning on subsets of Possible-D-SEP. This makes RFCI much faster
  than FCI. The orientation rules for v-structures and rule 4 were
  modified in order to produce an RFCI-PAG which, in the oracle version,
  is guaranteed to have the correct ancestral relationships.

  The first part of the RFCI algorithm is analogous to the PC and FCI
  algorithm. It starts with a complete undirected graph and estimates an
  initial skeleton using the function \code{\link{skeleton}}, which
  produces an initial order-independent skeleton, see
  \code{\link{skeleton}} for more details. All edges
  of this skeleton are of the form o-o. Due to the presence of hidden
  variables, it is no longer sufficient to consider only subsets of the
  neighborhoods of nodes \code{x} and \code{y} to decide whether the
  edge \code{x o-o y} should be removed. The FCI algorithm performs
  independence tests conditioning on subsets of Possible-D-SEP to remove
  those edges. Since this procedure is computationally infeasible, the
  RFCI algorithm uses a different approach to remove some of those
  superfluous edges before orienting the v-structures and the
  discriminating paths in orientation rule 4.

  Before orienting the v-structures, we perform the following additional
  conditional independence tests. For each unshielded triple a-b-c in
  the initial skeleton, we check if both a and b and b and c are
  conditionally dependent given the separating of a and c
  (sepset(a,c)). These conditional dependencies may not have been
  checked while estimating the initial skeleton, since sepset(a,c) does
  not need to be a subset of the neighbors of a nor of the neighbors of
  c. If both conditional dependencies hold and b is not in the
  sepset(a,c), the triple is oriented as a v-structure a->b<-c. On the
  other hand, if an additional conditional independence relationship may
  be detected, say a is independent from b given the sepset(a,c), the
  edge between a and c is removed from the graph and the set responsible
  for that is saved in sepset(a,b). The removal of an edge can destroy
  or create new unshielded triples in the graph. To solve this problem
  we work with lists (for details see Colombo et al., 2012).

  Before orienting discriminating paths, we perform the following additional
  conditional independence tests. For each triple a <-* b o- *c with a
  -> c,  the algorithm searches for a discriminating path p = <d,
  . . . , a,b,c> for b of minimal length, and checks that the vertices
  in every consecutive pair (f1,f2) on p are conditionally dependent
  given all subsets of 
  \eqn{\textrm{sepset}(d,c) \setminus \{f1, f2\}}{sepset(d,c) \ {f1, f2}}
  . If we do not find any
  conditional independence relationship, the path is oriented as in rule
  (R4). If one or more conditional independence relationships are found,
  the corresponding edges are removed, their minimal separating sets are
  stored.

  Conservative RFCI can be computed if the argument of \code{conservative} is
  \code{TRUE}. After the final skeleton is computed and the
  additional local tests on all unshielded triples, as described above,
  have been done, all potential v-structures a-b-c are checked
  in the following way. We test whether a and c are independent
  conditioning on any subset of the neighbors of a or any subset of the
  neighbors of c. When a subset makes a and c conditionally independent,
  we call it a separating set. If b is in no such separating set or in
  all such separating sets, no further action is taken and the normal
  version of the RFCI algorithm is continued. If, however, b is in only
  some separating sets, the triple a-b-c is marked 'ambiguous'. If a is
  independent of c given some S in the skeleton (i.e., the edge a-c
  dropped out), but a and c remain dependent given all subsets of
  neighbors of either a or c, we will call all triples a-b-c
  'unambiguous'. This is because in the RFCI algorithm, the true
  separating set might be outside the neighborhood of either a or c. An
  ambiguous triple is not oriented as a v-structure. Furthermore, no further
  orientation rule that needs to know whether a-b-c is a v-structure or
  not is applied. Instead of using the conservative version, which is
  quite strict towards the v-structures, Colombo and Maathuis (2014)
  introduced a less strict version for the v-structures called majority
  rule. This adaptation can be called using \code{maj.rule = TRUE}. In
  this case, the triple a-b-c is marked as 'ambiguous' if and only if b
  is in exactly 50 percent of such separating sets or no separating set
  was found. If b is in less than 50 percent of the separating sets it
  is set as a v-structure, and if in more than 50 percent it is set as a
  non v-structure (for more details see Colombo and Maathuis,
  2014).

  The implementation uses the stabilized skeleton
  \code{\link{skeleton}}, which produces an initial order-independent
  skeleton. The final skeleton and edge orientations can still be
  order-dependent, see Colombo and Maathuis (2014).
}

\references{
  D. Colombo and M.H. Maathuis (2014).Order-independent constraint-based
  causal structure learning. \emph{Journal of Machine Learning Research}
  \bold{15} 3741-3782. 

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson (2012).
  Learning high-dimensional directed acyclic graphs with latent and
  selection variables. \emph{Ann. Statist.} \bold{40}, 294-321.
}
\seealso{\code{\link{fci}} and \code{\link{fciPlus}} for estimating a
  PAG using the FCI algorithm; 
  \code{\link{skeleton}} for estimating an initial skeleton
  using the RFCI algorithm; \code{\link{pc}} for estimating a CPDAG using
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

## find skeleton and PAG using the RFCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
indepTest <- gaussCItest
res <- rfci(suffStat, indepTest, alpha = 0.9999, p=p, verbose=TRUE)

% The following: 1st part is *identical* to >>>>>>>> ./udag2apag.Rd
##################################################%  --------------
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## create the DAG :
V <- LETTERS[1:5]
edL <- setNames(vector("list", length = 5), V)
edL[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL[[2]] <- list(edges=3,weights=c(1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
## and leave  edL[[ 5 ]] empty
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
if (require(Rgraphviz))
  plot(g)

## define the latent variable
L <- 1

## compute the true covariance matrix of g
cov.mat <- trueCov(g)
## delete rows and columns belonging to latent variable L
true.cov <- cov.mat[-L,-L]
## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(true.cov)

## find PAG with RFCI algorithm
## as dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha :
rfci.pag <- rfci(suffStat = list(C = true.corr, n = 10^9),
                 indepTest = gaussCItest, alpha = 0.9999, labels = V[-L],
                 verbose=TRUE)

## define PAG given in Zhang (2008), Fig. 6, p.1882
corr.pag <- rbind(c(0,1,1,0),
                  c(1,0,0,2),
                  c(1,0,0,2),
                  c(0,3,3,0))
## check that estimated and correct PAG are in agreement:
stopifnot(corr.pag == rfci.pag@amat)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

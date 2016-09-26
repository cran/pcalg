\name{fci}
\alias{fci}
\title{Estimate a PAG by the FCI Algorithm}
\description{
  Estimate a Partial Ancestral Graph (PAG) from observational data, using
  the FCI (Fast Causal Inference) algorithm.
}

\usage{
fci(suffStat, indepTest, alpha, labels, p,
    skel.method = c("stable", "original", "stable.fast"),
    type = c("normal", "anytime", "adaptive"),
    fixedGaps = NULL, fixedEdges = NULL,
    NAdelete = TRUE, m.max = Inf, pdsep.max = Inf,
    rules = rep(TRUE, 10), doPdsep = TRUE, biCC = FALSE,
    conservative = FALSE, maj.rule = FALSE,
    numCores = 1, verbose = FALSE)
}

\arguments{
  \item{suffStat}{sufficient statistics: A named \code{\link{list}}
    containing all necessary elements for the conditional independence
    decisions in the function \code{indepTest}.}
  \item{indepTest}{a \code{\link{function}} for testing conditional
    independence.  The function is internally called as
    \code{indepTest(x,y, S, suffStat)}, and
    tests conditional independence of \code{x} and \code{y} given
    \code{S}.  Here, \code{x} and \code{y} are variables, and \code{S} is
    a (possibly empty) vector of variables (all variables are denoted
    by their column numbers
    in the adjacency matrix). \code{suffStat} is a list with all
    relevant information, see above.  The return value of
    \code{indepTest()} is the p-value of the test for conditional
    independence.}
  \item{alpha}{numeric significance level (in \eqn{(0, 1)}) for the
    individual conditional independence tests.}
  \item{labels}{(optional) \code{\link{character}} vector of variable (or
    \dQuote{node}) names.  Typically preferred to specifying \code{p}.}
  \item{p}{(optional) number of variables (or nodes).  May be specified
    if \code{labels} are not, in which case \code{labels} is set to
    \code{1:p}.}% FIXME: can get *both* from  suffStat$C  in the typical case
  \item{skel.method}{character string specifying method; the default,
    \code{"stable"} provides an \emph{order-independent} skeleton, see
    \code{\link{skeleton}}.}
  \item{type}{character string specifying the version of the FCI
    algorithm to be used.  By default, it is \code{"normal"}, and so the
    normal FCI algorithm is called.  If set to \code{"anytime"}, the
    \sQuote{Anytime FCI} is called and \code{m.max} needs to be specified.
    If set to \code{"adaptive"}, the \sQuote{Adaptive Anytime FCI} is
    called and \code{m.max} is not used.  For more information, see Details.}
  \item{fixedGaps}{\code{\link{logical}} matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    removed before starting the algorithm. Therefore, this edge is
    guaranteed to be absent in the resulting graph.}
  \item{fixedEdges}{logical matrix of dimension p*p. If entry
    \code{[i,j]} or \code{[j,i]} (or both) are TRUE, the edge i-j is
    never considered for removal. Therefore, this edge is
    guaranteed to be present in the resulting graph.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{m.max}{Maximum size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{pdsep.max}{Maximum size of Possible-D-SEP for which subsets are
    considered as conditioning sets in the conditional independence
    tests. If the nodes \code{x} and \code{y} are adjacent in the graph
    and the size of Possible-D-SEP(\code{x})\ {\code{x},\code{y}} is
    bigger than pdsep.max, the edge is simply left in the graph.
    Note that if pdsep.max is less than Inf, the final PAG
    may be a supergraph of the one computed with pdsep.max = Inf,
    because fewer tests may have been performed in the former.}
  \item{rules}{Logical vector of length 10 indicating which rules
    should be used when directing edges. The order of the rules is taken
    from Zhang (2009).}
  \item{doPdsep}{If \code{TRUE}, Possible-D-SEP is computed for all nodes,
    and all subsets of Possible-D-SEP are considered as conditioning
    sets in the conditional independence tests, if not defined otherwise
    in \code{pdsep.max}. If \code{FALSE}, Possible-D-SEP
    is not computed, so that the algorithm simplifies to the
    Modified PC algorithm of Spirtes, Glymour and Scheines (2000, p.84).}
  \item{biCC}{If \code{TRUE}, only nodes on paths between nodes \code{x}
    and \code{y} are considered to be in Possible-D-SEP(\code{x}) when
    testing independence between \code{x} and \code{y}. Uses
    biconnected components, \code{\link[RBGL]{biConnComp}} from \pkg{RBGL}.}
  \item{conservative}{Logical indicating if the unshielded triples
    should be checked for ambiguity the second time when v-structures
    are determined.  For more information, see details.}
  \item{maj.rule}{Logical indicating if the unshielded triples should
    be checked for ambiguity the second time when v-structures are
    determined using a majority rule idea, which is less strict than the
    standard conservative.  For more information, see details.}
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
  This function is a generalization of the PC algorithm (see
  \code{\link{pc}}), in the sense that it allows arbitrarily many
  latent and selection variables. Under the assumption that the data
  are faithful to a DAG that includes all latent and selection variables,
  the FCI algorithm (Fast Causal Inference algorithm) (Spirtes, Glymour
  and Scheines, 2000) estimates the Markov equivalence class of MAGs
  that describe the conditional independence relationships between the
  observed variables.

  We estimate an equivalence class of \bold{maximal ancestral graph}s (\bold{MAG}s)
  instead of DAGs, since DAGs are not closed under marginalization and
  conditioning (Richardson and Spirtes, 2002).

  An equivalence class of a MAG can be uniquely represented by a \bold{partial
  ancestral graph} (\bold{PAG}).  A PAG contains the following types of edges:
  o-o, o-, o->, ->, <->, -. The bidirected edges come from hidden variables,
  and the undirected edges come from selection variables. The edges have
  the following interpretation: (i) there is an edge between \code{x}
  and \code{y} if and only if variables {x} and {y} are conditionally dependent
  given S for all sets S consisting of all selection variables
  and a subset of the observed variables; (ii) a tail on an edge means
  that this tail is present in all MAGs in the Markov equivalence class;
  (iii) an arrowhead on an edge means that this arrowhead is present in
  all MAGs in the Markov equivalence class; (iv) a o-edgemark means that
  there is a at least one MAG in the Markov equivalence class where the
  edgemark is a tail, and at least one where the edgemark is an
  arrowhead. Information on the interpretation of edges in a MAG can be
  found in the references given below.

  The first part of the FCI algorithm is analogous to the PC algorithm. It
  starts with a complete undirected graph and estimates an initial skeleton
  using \code{\link{skeleton}(*, method="stable")} which produces an
  initial order-independent skeleton, see \code{\link{skeleton}} for
  more details.  All edges of this skeleton are of
  the form o-o.  Due to the presence of hidden variables, it is no longer
  sufficient to consider only subsets of the neighborhoods of nodes \code{x}
  and \code{y} to decide whether the edge \code{x o-o y} should be removed.
  Therefore, the initial skeleton may contain some superfluous edges.
  These edges are removed in the next step of the algorithm which
  requires some orientations.  Therefore, the v-structures
  are determined using the conservative method (see discussion on
  \code{conservative} below).

  After the v-structures have been oriented, Possible-D-SEP sets for each
  node in the graph are computed at once. To decide whether edge
  \code{x o-o y} should be removed, one performs conditional indepedence
  tests of x and y given all possible subsets of Possible-D-SEP(x) and
  of Possible-D-SEP(y). The edge is removed if a conditional
  independence is found. This produces a fully order-independent final
  skeleton as explained in Colombo and Maathuis (2014). Subsequently,
  the v-structures are newly determined on the final skeleton (using
  information in sepset). Finally, as many as possible undetermined edge
  marks (o) are determined using (a subset of) the 10 orientation rules
  given by Zhang (2009).

  The \dQuote{Anytime FCI} algorithm was introduced by Spirtes (2001).  It
  can be viewed as a modification of the FCI algorithm that only performs
  conditional independence tests up to and including order m.max when
  finding the initial skeleton, using the function \code{skeleton}, and
  the final skeleton, using the function \code{pdsep}.  Thus, Anytime FCI
  performs fewer conditional independence tests than FCI.  To use the
  Anytime algorithm, one sets \code{type = "anytime"} and needs to
  specify \code{m.max}, the maximum size of the conditioning sets.

  The \dQuote{Adaptive Anytime FCI} algorithm was introduced by Colombo
  et. al (2012).  The first part of the algorithm is identical to the normal
  FCI described above.  But in the second part when the final skeleton is
  estimated using the function \code{\link{pdsep}}, the Adaptive Anytime
  FCI algorithm only performs conditional independence tests up to and
  including order \code{m.max}, where m.max is the maximum size of the
  conditioning sets that were considered to determine the initial
  skeleton using the function \code{skeleton}.  Thus, m.max is chosen
  adaptively and does not have to be specified by the user.

  \emph{Conservative} versions of FCI, Anytime FCI, and Adaptive Anytime FCI
  are computed if \code{conservative = TRUE} is specified.  After the
  final skeleton is computed, all potential
  v-structures a-b-c are checked in the following way. We test whether a
  and c are independent conditioning on any subset of the neighbors of a
  or any subset of the neighbors of c. When a subset makes a and c
  conditionally independent, we call it a separating set. If b is in no
  such separating set or in all such separating sets, no further action
  is taken and the normal version of the FCI, Anytime FCI, or Adaptive
  Anytime FCI algorithm is continued. If, however, b is in only some
  separating sets, the triple a-b-c is marked \sQuote{ambiguous}. If a is
  independent of c given some S in the skeleton (i.e., the edge a-c
  dropped out), but a and c remain dependent given all subsets of
  neighbors of either a or c, we will call all triples a-b-c
  \sQuote{unambiguous}. This is because in the FCI algorithm, the true separating set
  might be outside the neighborhood of either a or c. An ambiguous
  triple is not oriented as a v-structure. Furthermore, no further
  orientation rule that needs to know whether a-b-c is a v-structure or
  not is applied. Instead of using the conservative version, which is
  quite strict towards the v-structures, Colombo and Maathuis (2014)
  introduced a less strict version for the v-structures called majority
  rule. This adaptation can be called using \code{maj.rule = TRUE}. In
  this case, the triple a-b-c is marked as \sQuote{ambiguous} if and only if b
  is in exactly 50 percent of such separating sets or no separating set
  was found. If b is in less than 50 percent of the separating sets it
  is set as a v-structure, and if in more than 50 percent it is set as a
  non v-structure (for more details see Colombo and Maathuis,
  2014). Colombo and Maathuis (2014) showed that with both these
  modifications, the final skeleton and the decisions about the
  v-structures of the FCI algorithm are fully order-independent.

  Note that the order-dependence issues on the 10 orientation rules are
  still present, see Colombo and Maathuis (2014) for more details.
}
\references{
  D. Colombo and M.H. Maathuis (2014). Order-independent
  constraint-based causal structure learning. Journal of Machine
  Learning Research 15 3741-3782. 

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294-321.

  M.Kalisch, M. Maechler, D. Colombo, M. H. Maathuis,
  P. Buehlmann (2012). Causal Inference Using Graphical Models with
  the R Package pcalg. \emph{Journal of Statistical Software}
  \bold{47(11)} 1--26, \url{http://www.jstatsoft.org/v47/i11/}.

  P. Spirtes (2001). An anytime algorithm for causal inference. \emph{In
  Proc. of the Eighth International Workshop on Artifiial Intelligence
  and Statistics} 213-221. Morgan Kaufmann, San Francisco.

  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition, MIT Press,
  Cambridge (MA).

  P. Spirtes, C. Meek, T.S. Richardson (1999). In: \emph{Computation,
  Causation and Discovery}. An algorithm for causal
  inference in the presence of latent variables and selection bias.
  Pages 211-252. MIT Press.

  T.S. Richardson and P. Spirtes (2002). Ancestral graph Markov models.
  \emph{Annals of Statistics} \bold{30} 962-1030.

  J. Zhang (2008). On the completeness of orientation rules for
  causal discovery in the presence of latent confounders and selection bias.
   \emph{Artificial Intelligence} \bold{172} 1873-1896.
}
\seealso{\code{\link{fciPlus}} for a more efficient variation of FCI;
  \code{\link{skeleton}} for estimating a skeleton 
  using the PC algorithm; \code{\link{pc}} for estimating a CPDAG using
  the PC algorithm; \code{\link{pdsep}} for computing
  Possible-D-SEP for each node and testing and adapting the graph
  accordingly; \code{\link{qreach}} for a fast way of finding Possible-D-SEP
  for a given node.

  \code{\link{gaussCItest}}, \code{\link{disCItest}},
  \code{\link{binCItest}} and \code{\link{dsepTest}} as examples for
  \code{indepTest}.
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

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
res <- fci(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p, doPdsep = FALSE)
%FIXME showEdgeList(res)

##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## create the graph g
p <- 4
L <- 1 # '1' is latent
V <- c("Ghost", "Max","Urs","Anna","Eva")
edL <- setNames(vector("list", length=length(V)), V)
edL[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL[[2]] <- list(edges=3,weights=c(1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## delete rows and columns belonging to latent variable L
true.cov <- cov.mat[-L,-L]
## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(true.cov)

## The same, for the following three examples
indepTest <- gaussCItest
suffStat <- list(C = true.corr, n = 10^9)

## find PAG with FCI algorithm.
## As dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:

normal.pag <- fci(suffStat, indepTest, alpha = 0.9999, labels = V[-L],
                  verbose=TRUE)
%FIXME showEdgeList(normal.pag)

## find PAG with Anytime FCI algorithm with m.max = 1
## This means that only conditioning sets of size 0 and 1 are considered.
## As dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
anytime.pag <- fci(suffStat, indepTest, alpha = 0.9999, labels = V[-L],
                   type = "anytime", m.max = 1,
                   verbose=TRUE)

## find PAG with Adaptive Anytime FCI algorithm.
## This means that only conditining sets up to size K are considered
## in estimating the final skeleton, where K is the maximal size of a
## conditioning set found while estimating the initial skeleton.
## As dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
adaptive.pag <- fci(suffStat, indepTest, alpha = 0.9999, labels = V[-L],
                    type = "adaptive",
                    verbose=TRUE)

## define PAG given in Zhang (2008), Fig. 6, p.1882
corr.pag <- rbind(c(0,1,1,0),
                  c(1,0,0,2),
                  c(1,0,0,2),
                  c(0,3,3,0))

## check if estimated and correct PAG are in agreement
all(corr.pag == normal.pag  @ amat) # TRUE
all(corr.pag == anytime.pag @ amat) # FALSE
all(corr.pag == adaptive.pag@ amat) # TRUE
ij <- rbind(cbind(1:4,1:4), c(2,3), c(3,2))
all(corr.pag[ij] == anytime.pag @ amat[ij]) # TRUE

\dontshow{stopifnot(
  corr.pag == normal.pag  @ amat
  ,
  corr.pag[ij] == normal.pag@amat[ij]
  ,
  corr.pag == adaptive.pag@ amat
  )}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

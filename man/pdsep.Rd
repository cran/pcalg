\name{pdsep}
\alias{pdsep}
\title{Estimate Final Skeleton in the FCI algorithm}
\description{
  Estimate the final skeleton in the FCI algorithm (Spirtes et al, 2000),
  as described in Steps 2 and 3 of Algorithm 3.1 in Colombo et
  al. (2012).  The input of this function consists of an initial
  skeleton that was estimated by the PC algorithm (Step 1 of Algorithm
  3.1 in Colombo et al. (2012)).

  Given the initial skeleton, all unshielded triples are considered and
  oriented as colliders when appropriate.  Then, for all nodes x in the
  resulting partially directed graph G, Possible-D-SEP(x,G) is computed,
  using the function \code{qreach}.  Finally, for any edge y-z that is
  present in G, conditional independence between Y and Z is tested given
  all subsets of Possible-D-SEP(y,G) and all subsets of
  Possible-D-SEP(z,G).  These tests are done at level alpha, using
  \code{indepTest}.  If the pair of nodes is judged to be independent
  given some set S, then S is recorded in sepset(y,z) and sepset(z,y)
  and the edge y-z is deleted.  Otherwise, the edge remains and there is
  no change to sepset.
}
\usage{
pdsep(skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
      pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
      biCC = FALSE, verbose = FALSE)
}
\arguments{
  \item{skel}{Graph object returned by \code{\link{skeleton}}.}
  \item{suffStat}{Sufficient statistic: A list containing all necessary
    elements for making conditional independence decisions using
    function \code{indepTest}.}
  \item{indepTest}{Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)} for
    testing conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are node numbers of the
    adjacency matrix, \code{S} is a (possibly empty) vector of node
    numbers of the adjacency matrix and \code{suffStat} is a list
    containing all relevant elements for making conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables.}
  \item{sepset}{List of length \code{p}; each element of the list
    contains another list of length \code{p}.  The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out.  This object is thought to be
    obtained from a \code{pcAlgo}-object or \code{fciAlgo}-object.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{pMax}{Matrix with the maximal p-values of conditional
    independence tests in a previous call of \code{\link{skeleton}},
    \code{\link{pc}} or \code{\link{fci}} which produced \code{G}. This
    object is thought to be obtained from a \code{pcAlgo}-object or
    \code{fciAlgo}-object.}
  \item{m.max}{Maximum size of the conditioning sets that are considered in the
    conditional independence tests.}
  \item{pdsep.max}{Maximum size of Possible-D-SEP for which subsets are
    considered as conditioning sets in the conditional independence
    tests. If the nodes \code{x} and \code{y} are adjacent in the graph
    and the size of Possible-D-SEP(x,G)\ {\code{x},\code{y}}, is
    bigger than pdsep.max, the edge is simply left in the graph.
    Note that if pdsep.max is less than Inf, the final PAG
    is typically a supergraph of the one computed with pdsep.max = Inf,
    because fewer tests may have been performed in the former.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{unfVect}{Vector containing numbers that encode the unfaithful
    triple (as returned by \code{\link{pc.cons.intern}}). This is
    needed in the conservative FCI.}
  \item{biCC}{Logical; if \code{TRUE}, only nodes on paths between nodes
    a and c are considered to be in sepset(a,c).  This uses biconnected
    components, see \code{\link[RBGL]{biConnComp}} from \pkg{RBGL}.}
  \item{verbose}{Logical indicating that detailed output is to be provided.}
}
\value{A list with the following elements:
  \item{G}{Updated adjacency matrix representing the final skeleton}
  \item{sepset}{Updated sepsets}
  \item{pMax}{Updated matrix containing maximal p-values}
  \item{allPdsep}{Possible-D-Sep for each node}
  \item{max.ord}{Maximal order of conditioning sets during independence
    tests}
  \item{n.edgetests}{Number of conditional edgetests performed, grouped
    by the size of the conditioning set.}
}
\details{
  To make the code more efficient, we only perform tests that were not
  performed in the estimation of the initial skeleton.

  Note that the Possible-D-SEP sets are computed once in the beginning. They are  not updated after edge deletions, in order to make sure that the output of the algorithm does not depend on the ordering of the variables (see also Colombo and Maathuis (2013)).
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition. The MIT Press.

  D. Colombo, M.H. Maathuis, M. Kalisch and T.S. Richardson (2012).
  Learning high-dimensional directed acyclic graphs with latent and
  selection variables.
  \emph{Annals of Statistics} \bold{40}, 294--321.

  D. Colombo and M.H. Maathuis (2013).
  Order-independent constraint-based causal structure learning. arXiv:1211.3295v2.
}
\seealso{\code{\link{qreach}} to find Possible-D-SEP(x,G);
  \code{\link{fci}}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
}
\examples{
p <- 10
## generate and draw random DAG:
set.seed(44)
myDAG <- randomDAG(p, prob = 0.2)

## generate 10000 samples of DAG using gaussian distribution
library(RBGL)
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton
indepTest <- gaussCItest
suffStat <- list(C = cor(d.mat), n = n)
alpha <- 0.01
skel <- skeleton(suffStat, indepTest, alpha=alpha, p=p)

## prepare input for pdsep
sepset <- skel@sepset
pMax <- skel@pMax

## call pdsep to find Possible-D-Sep and enhance the skeleton
pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p, sepset, alpha,
                  pMax, verbose = TRUE)
## call pdsep with biconnected components to find Possible-D-Sep and enhance the skeleton
pdsepResBicc <- pdsep(skel@graph, suffStat, indepTest, p, sepset, alpha,
                      pMax, biCC= TRUE, verbose = TRUE)
}

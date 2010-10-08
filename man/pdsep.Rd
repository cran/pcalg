\name{pdsep}
\alias{pdsep}
\title{Compute Possible-D-Sep for Each Node and Adapt Graph Accordingly}
\description{
  Computes the Possible-D-Sep for each node
  and tests conditional independence given all subsets of
  Possible-D-Sep.  The graph is updated accordingly.
}
\usage{
pdsep(skel, suffStat, indepTest, p, sepset,
      pMax, NAdelete = TRUE, verbose = FALSE,
      alpha, unfVect = NULL, biCC = FALSE)
}
\arguments{
  \item{skel}{Graph object returned by \code{\link{skeleton}}.}
  \item{suffStat}{sufficient statistic: A list containing all necessary
    elements for making conditional independence decisions using
    function \code{indepTest}.}
  \item{indepTest}{predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)} for
    testing conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are node numbers of the
    adjacency matrix, \code{S} is a (possibly empty) vector of node
    numbers of the adjacency matrix and \code{suffStat} is a list
    containing all relevant elements for making conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{number of variables}
  \item{sepset}{List of length \code{p}; each element of the list
    contains another list of length \code{p}. The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out. This object is thought to be
    obtained from a \code{pcAlgo}-object or \code{fciAlgo}-object.}
  \item{pMax}{Matrix with the maximal p-values of conditional
    independence tests in a previous call of \code{\link{skeleton}},
    \code{\link{pc}} or \code{\link{fci}} which produced \code{G}. This
    object is thought to be obtained from a \code{pcAlgo}-object or
    \code{fciAlgo}-object.}
  \item{NAdelete}{If indepTest returns \code{NA} and this option is
    \code{TRUE}, the corresponding edge is deleted. If this option is
    \code{FALSE}, the edge is not deleted.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{alpha}{Significance level for the individual conditional
    independence tests}
  \item{unfVect}{vector containing numbers that encode the unfaithful
    triple (as returned by \code{\link{pc.cons.intern}}).  This is
    needed in the conservative FCI.}
  \item{biCC}{logical; if true, only nodes on some path between a and c are
    considered to be in sepset(a,c).  Uses biconnected components,
    \code{\link[RBGL]{biConnComp}} from \pkg{RBGL}.}
}
\value{A list with several elements:
  \item{G}{Updated adjecency matrix}
  \item{sepset}{Updated sepset}
  \item{pMax}{Updated matrix containing maximal p-values}
  \item{allPdsep}{Possible d-sep for each node}
  \item{max.ord}{Maximal order of conditioning sets during independence
    tests}
  \item{n.edgetests}{Number of conditional edgetests performed grouped
    by size of conditioning set.}
}

\details{For a given graph G, a node y is in Possible-D-Sep(x) if x!=y
  and there is an undirected path U between x and y in G such that for
  every subpath a,b,c of U either b is a collider on the subpath, or b
  is not a definite noncollider on U and a,b and c form a trianglein G.

  Each pair of nodes x and y which is connected by
  an edge is tested for conditional independence given every subset in
  Possible-D-Sep of x or of y. The conditional independence is tested on
  significance level alpha by using the test given in
  \code{indepTest}. If the pair of nodes is judged to be independent
  given set S, then S is recorded in the sepset of x,y and in the sepset
  of y,x and the edge is deleted. Otherwise, the edge is not deleted and
  no entry to the sepset is made.

  To make the code more efficient, tests that have been done during
  finding the skeleton are not done again.}

\references{
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition. The MIT Press.
}
\seealso{\code{\link{qreach}} for a fast way of finding Possible-D-Sep
  for a given node; \code{\link{fci}} which uses \code{\link{pdsep}}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
p <- 10
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.2)
## generate 10000 samples of DAG using gaussian distribution
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## estimate skeleton
indepTest <- gaussCItest
suffStat <- list(C = cor(d.mat), n = n)
alpha <- 0.01
skel <- skeleton(suffStat, indepTest, p, alpha)

## prepare input for pdsep
sepset <- skel@sepset
pMax <- skel@pMax
n.edgetestsSKEL <- skel@n.edgetests
max.ordSKEL <- skel@max.ord

## call pdsep to find possible d-sep and enhance the skeleton
pdsepRes <- pdsep(skel@graph, suffStat, indepTest, p, sepset, pMax, NAdelete,
verbose = TRUE, alpha)
}

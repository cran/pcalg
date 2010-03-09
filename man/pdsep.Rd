\name{pdsep}
\alias{pdsep}
\title{Compute possible d-sep for each node and test and adapt graph accordingly}

\description{ This function computes the possible d-sep for each
  node. Then, each pair of nodes x and y which is connected by an edge
  is tested for conditional independence given every subset in possible
  d-sep of x or of y. The conditional independence is tested on
  significance level alpha by using the test given in
  \code{indepTest}. If the pair of nodes is judged to be independent
  given set S, then S is recorded in the sepset of x,y and in the sepset
  of y,x and the edge is deleted. Otherwise, the edge is not deleted and
  no entry to the sepset is made.}

\usage{
pdsep(G, suffStat, indepTest, p, sepset, pMax, NAdelete = TRUE, verbose
= FALSE, alpha)
}
\arguments{
  \item{G}{Boolean adjacency matrix (TRUE if there is an edge)}
  \item{suffStat}{Sufficient statistics: A list containing all necessary
    elements for making conditional independence decisions using
    function \code{indepTest}}
  \item{indepTest}{Predefined function for testing conditional independence. The
    function is internally called as \code{indepTest(x,y,S,suffStat)} for
    testing conditional independence of \code{x} and \code{y} given
    \code{S}. Here, \code{x} and \code{y} are node numbers of the
    adjacency matrix, \code{S} is a (possibly empty) vector of node
    numbers of the adjacency matrix and \code{suffStat} is a list
    containing all relevant elements for making conditional independence
    decisions. The return value of \code{indepTest} is the p-value of
    the test for conditional independence.}
  \item{p}{Number of variables}
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

\details{To make the code more efficient, only condtioning sets that
  have not been tested in a previous call of \code{skeleton} are
  tested. For example, subsets of the adjacency of a node without
  additional nodes from possible d-sep are not considered alone, but
  only in combination with nodes from possible d-sep. Within the nodes
  of possible d-sep that are not in the adjacency of a node, all
  possible subsets are tested.}

\references{
  S.L. Lauritzen (2004)
  \emph{Graphical Models}, Oxford University Press, Chapter 3.2.2
}

\seealso{\code{\link{qreach}} for a fast way of finding possible d-sep
    for a given node}

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
G <- (as(skel@graph,"matrix")!=0)
sepset <- skel@sepset
pMax <- skel@pMax
n.edgetestsSKEL <- skel@n.edgetests
max.ordSKEL <- skel@max.ord

## call pdsep to find possible d-sep and enhance the skeleton
pdsepRes <- pdsep(G, suffStat, indepTest, p, sepset, pMax, NAdelete,
verbose = TRUE, alpha)
}

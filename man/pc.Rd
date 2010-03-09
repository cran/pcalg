\name{pc}
\alias{pc}
\title{Estimate the equivalence class of a DAG (CPDAG) 
  using the PC algorithm}
\description{
  Estimate the equivalence class of a directed acyclic
  graph (DAG) from observational data, using the PC-algorithm. 
}
\usage{
pc(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = "rand")
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
  \item{u2pd}{Method for dealing with conflicting information when
    trying to orient edges (see details below).}
}
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the equivalence 
  class of the underlying DAG. 
}
\details{
  Under the assumption that the distribution of the observed variables 
  is faithful to a DAG, this function estimates the equivalence class of 
  the DAG. We do not estimate the DAG itself, because this is typically 
  impossible (even with an infinite amount of data), since different 
  DAGs can describe the same conditional independence relationships. 
  Since all DAGs in an equivalence class describe the same conditional
  independence relationships, they are equally valid ways to
  describe the conditional dependence structure that was given as
  input.

  All DAGs in an equivalence class have the same skeleton (i.e.,
  the same adjacency information) and the same v-structures (see
  definition below). However, the direction of some edges may be
  undetermined, in the sense that they point one way in one DAG in the
  equivalence class, while they point the other way in another DAG in
  the equivalence class.

  An equivalence class can be uniquely represented by a completed 
  partially directed acyclic graph (CPDAG). A CPDAG
  contains undirected and directed edges. The edges have the following 
  interpretation: (i) there is a (directed or undirected) edge between i 
  and j if and only if variables i and j are conditionally dependent 
  given S for all possible subsets S of the remaining nodes; (ii) a directed
  edge i->j means that this directed edge is present in all DAGs in the 
  equivalence class; (iii) an undirected edge i-j means that there is at least 
  one DAG in the equivalence class with edge i->j and there is at least one DAG
  in the equivalence class with edge i<-j. 
  
  The CPDAG is estimated using the PC algorithm (named after its inventors
  \bold{P}eter Spirtes and \bold{C}lark Glymour). The skeleton is 
  estimated by the function \code{\link{skeleton}} (see the help file 
  of this function for details). Subsequently, as many edges as possible
  are oriented. This is done in two steps.

  First, the algorithm considers all triples 
  \code{(a,b,c)}, where \code{a} and
  \code{b} are adjacent, \code{b} and \code{c} are adjacent, but
  \code{a} and \code{c} are not adjacent. For all such triples, 
  we direct both edges towards b (a->b<-c) if and only if \code{b} was 
  not part of the conditioning set that made the edge between \code{a} 
  and \code{c} drop out. These conditioning sets were saved in
  \code{sepset}. The structure a->b<-c is called a v-structure.

  After determining all v-structures, there may still
  be undirected edges. It may be possible to direct some of these edges, since
  one can deduce that one of the two possible directions of the edge is 
  invalid because it introduces 
  a new v-structure or a directed cycle. Such edges are found by 
  repeatedly applying rules R1-R3 as given in Algorithm 2 of Kalisch 
  and B\"uhlmann (2007). The algorithm stops if none of the rules is 
  applicable to the graph. 

  Sampling errors (or hidden variables) can lead to conflicting
  information about edge directions. For example, one may find that
  a-b-c and b-c-d should both be directed as v-structures.  This gives
  conflicting information about the edge b-c, since it should be
  directed as b<-c in v-structure a->b<-c, while it should be directed
  as b->c in v-structure b->c<-d. In such cases, we simply overwrite the
  directions of the conflicting edge. In the example above this means
  that we obtain a->b->c<-d if a-b-c was visited first, and a->b<-c<-d
  if b-c-d was visited first.

  Sampling errors or hidden variables can also lead to invalid CPDAGs, 
  meaning that there does not exist a DAG
  that has the same skeleton and v-structures as the 
  graph found by the algorithm. An example of this is an undirected cycle 
  consisting of the 
  edges a-b-c-d and d-a. In this case it is impossible to direct the edges 
  without creating a cycle or a new v-structure. The option \code{u2pd} 
  specifies what should be done in such a situation. If the option is set to 
  \code{"relaxed"}, the algorithm simply outputs the invalid CPDAG. If
  the option is set to \code{"rand"},
  all direction information is discarded and a random DAG is generated on the
  skeleton. If the option is set to \code{"retry"}, up to 100 combinations 
  of possible directions of the ambiguous edges are tried, and the first 
  combination that results in a valid CPDAG is chosen. If no valid combination
  is found, an arbitrary DAG is generated on the skeleton 
  as in the option "rand".
  
  Notes: (1) Throughout, the algorithm works with the column positions of 
  the variables in the adjacency matrix, and not with the names of 
  the variables. (2) When plotting the object, undirected and bidirected 
  edges are equivalent.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.

  M. Kalisch and P. B\"uhlmann (2007)
  \emph{Estimating high-dimensional
    directed acyclic graphs with the PC-algorithm},
    JMLR \bold{8} 613-636.
}
\seealso{\code{\link{skeleton}} for estimating a skeleton of a DAG; 
  \code{\link{udag2pdag}} for converting the
  skeleton to a CPDAG; \code{\link{gaussCItest}},
  \code{\link{disCItest}}, \code{\link{binCItest}} and
  \code{\link{dsepTest}} as examples for \code{indepTest}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler.
}
\examples{
## generate a random DAG:
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)

## plot the DAG
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

##################################################
## Using Gaussian Data
##################################################
## generate 10000 samples from the DAG using Gaussian errors
n <- 10000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## define independence test (partial correlations)
indepTest <- gaussCItest 
## define sufficient statistics
suffStat <- list(C = cor(d.mat), n = n)
## estimate CPDAG
alpha <- 0.01
pc.fit <- pc(suffStat, indepTest, p, alpha)
## show estimated CPDAG
plot(pc.fit, main = "Estimated CPDAG")
x11(); plot(myDAG)

##################################################
## Using d-separation oracle
##################################################
## define independence test
indepTest <- dsepTest
## define sufficient statistics (d-separation oracle)
suffStat <- list(g = myDAG, jp = johnson.all.pairs.sp(myDAG))
## estimate CPDAG
alpha <- 0.01 ## value is irrelevant as dsepTest returns either 0 or 1
fit <- pc(suffStat, indepTest, p, alpha)
## show estimated CPDAG
plot(fit, main = "Estimated CPDAG")

##################################################
## Using discrete data
##################################################
data(discreteData)
p <- ncol(dat)
## define independence test (G^2 statistics)
indepTest <- disCItest
## define sufficient statistics
suffStat <- list(dm = dat, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate CPDAG
alpha <- 0.01 
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)
## show estimated CPDAG
plot(pc.fit, main = "Estimated CPDAG")
rm(dat)

##################################################
## Using binary data
##################################################
data(binaryData)
p <- ncol(dat)
## define independence test
indepTest <- binCItest
## define sufficient statistics
suffStat <- list(dm = dat, adaptDF = FALSE)
## estimate CPDAG
alpha <- 0.01 
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)
## show estimated CPDAG
plot(pc.fit, main = "Estimated CPDAG")
rm(dat)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

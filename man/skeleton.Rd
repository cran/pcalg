\name{skeleton}
\alias{skeleton}
\title{Estimate the skeleton of a DAG using the PC Algorithm}
\description{
  Estimate the \dQuote{skeleton} of a directed 
  acyclic graph (DAG) from observational data, 
  using the PC-algorithm.  
}
\usage{
skeleton(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
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
\value{An object of \code{\link{class}} \code{"pcAlgo"} (see
  \code{\linkS4class{pcAlgo}}) containing an estimate of the skeleton of
  the underlying DAG,
  the conditioning sets that led to edge removals (sepset) and several
  other parameters.
}
\details{
  Under the assumption that the distribution of the observed variables 
  is faithful to a DAG, this function estimates the skeleton 
  of the DAG. The skeleton of a DAG is the undirected graph 
  resulting from removing all arrowheads from the DAG. Edges in the skeleton of
  a DAG have the following interpretation: there is an edge between i and j 
  if and only if variables i and j are conditionally dependent given S for 
  all possible subsets S of the remaining nodes. 

  The data are not required to follow a specific distribution, 
  but one should make 
  sure that the conditional indepedence test used in \code{indepTest} 
  is appropriate for the data. Pre-programmed versions of \code{indepTest} 
  are available for Gaussian data (\code{\link{gaussCItest}}), discrete data 
  (\code{\link{disCItest}}), and binary data (see \code{\link{binCItest}}). 
  Users can also specify their own \code{indepTest} function.   

  The PC algorithm (Spirtes, Glymour and Scheines, 2000) starts with a 
  complete undirected graph. In each step, it visits 
  all pairs (\code{i},\code{j}) of adjacent nodes in the current graph, and 
  determines based on conditional independence tests whether the edge i-j 
  should be removed. In particular, in step \code{m} (\code{m}=0,1,...)  
  the algorithm visits all 
  pairs (\code{i}, \code{j}) of adjacent nodes in the current
  graph, and the edge between
  \code{i} and \code{j} is kept if and only if the null hypothesis 
  "\code{i} and \code{j} are
  conditionally independent given S" is
  rejected at significance level \code{alpha} for all subsets \code{S} of size
  \code{m} of the neighbours of \code{i} and of the neighbors of \code{j} in the
  current graph (as judged by the function \code{indepTest}).
  The algorithm stops when m is larger than the
  largest neighbourhood size of all nodes, or when m has reached 
  the limit \code{m.max} defined by the user.

  The information in \code{fixedGaps} and \code{fixedEdges} is used as follows. 
  The gaps given in \code{fixedGaps} are introduced in the very beginning of 
  the algorithm by removing the corresponding edges from the complete 
  undirected graph. Pairs  (\code{i},\code{j}) in \code{fixedEdges} are 
  skipped in all steps of the algorithm, so that these edges remain in the 
  graph. 

  Note: Throughout, the algorithm works with the column positions of 
  the variables in the adjacency matrix, and not with the names of 
  the variables.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, MIT Press.

  M. Kalisch and P. B\"uhlmann (2007)
  \emph{Estimating high-dimensional
    directed acyclic graphs with the PC-algorithm},
    JMLR \bold{8} 613-636.
}

\seealso{\code{\link{pc}} for generating a partially directed graph
  using the PC algorithm; \code{\link{udag2pdag}} for converting the
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
set.seed(46)
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
## estimate skeleton
alpha <- 0.01
skel.fit <- skeleton(suffStat, indepTest, p, alpha)
## show estimated skeleton
plot(skel.fit, main = "Estimated Skeleton")
x11(); plot(myDAG)

##################################################
## Using d-separation oracle
##################################################
## define independence test
indepTest <- dsepTest
## define sufficient statistics (d-separation oracle)
suffStat <- list(g = myDAG, jp = johnson.all.pairs.sp(myDAG))
## estimate skeleton
alpha <- 0.01 ## value is irrelevant as dsepTest returns either 0 or 1
skel.fit <- skeleton(suffStat, indepTest, p, alpha)
## show estimated skeleton
plot(skel.fit, main = "Estimated Skeleton")

##################################################
## Using discrete data
##################################################
data(discreteData)
p <- ncol(dat)
## define independence test (G^2 statistics)
indepTest <- disCItest
## define sufficient statistics
suffStat <- list(dm = dat, nlev = c(3,2,3,4,2), adaptDF = FALSE)
## estimate skeleton 
alpha <- 0.01 
skel.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = TRUE)
## show estimated skeleton
plot(skel.fit, main = "Estimated Skeleton")
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
## estimate skeleton
alpha <- 0.01 
skel.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = TRUE)
## show estimated skeleton
plot(skel.fit, main = "Estimated Skeleton")
rm(dat)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

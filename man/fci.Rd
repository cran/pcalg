\name{fci}
\alias{fci}
\title{Estimate the equivalence class of a MAG (PAG) using the FCI Algorithm}
\description{
  Estimate the equivalence class of a maximal ancestral graph (MAG) from
  observational data, using the FCI-algorithm.
}

\usage{
fci(suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL,
    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10),
    doPdsep = TRUE, conservative = c(FALSE, FALSE), biCC = FALSE,
    cons.rules = FALSE, labels = NA)
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
  \item{rules}{Logical vector of length 10 indicating which rules
    should be used when directing edges. The order of the rules is taken
    from Zhang (2009).}
  \item{doPdsep}{If TRUE, Possible-D-SEP is computed for all nodes,
    and all subsets of Possible-D-SEP are considered as conditioning
    sets in the conditional independence tests. If FALSE, Possible-D-SEP
    is not computed, so that the algorithm simplifies to the
    Modified PC algorithm of Spirtes, Glymour and Scheines (2000, page
    ...).}
  \item{conservative}{If the first argument is TRUE, the skeleton is
  computed in a conservative way. If the second argument is TRUE, the
  triples are checked for faithfulness again after possible deletion of
  edges when finding Possible-D-SEP. For more information, see Details.}
  \item{biCC}{If TRUE, only nodes on paths between a and c are
  considered to be in sepset(a,c). Uses biconnected components.}
  \item{cons.rules}{If TRUE, an orientation rule that needs information
  on definite non-colliders is only applied, if the corresponding
  subgraph relevant for the rule does not involve an unfaithful triple.}
\item{labels}{Vector of length p with node names. If \code{NA},
  \code{as.character{1:p}} is used instead.}
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
  the FCI algorithm (Fast Causal Inference algorithm)
  estimates the equivalence class of MAGs that describe
  the conditional independence relationships between the observed variables.

  We estimate an equivalence class of MAGs instead of DAGs, since DAGs are
  not closed under marginalization and conditioning
  (Richardson and Spirtes, 2002).

  An equivalence class of a MAG can be uniquely represented by a partial
  ancestral graph (PAG). A PAG contains the following types of edges:
  o-o, o-, o->, ->, <->, -. The bidirected edges come from hidden variables,
  and the undirected edges come from selection variables. The edges have
  the following interpretation: (i) there is an edge between \code{x}
  and \code{y} if and only if variables {x} and {y} are conditionally dependent
  given S for all sets S consisting of all selection variables
  and a subset of the observed variables; (ii) a tail on an edge means
  that this tail is present in all MAGs in the equivalence class;
  (iii) an errowhead on an  edge means that this arrowhead is present in
  all MAGs in the equivalence
  class; (iv) a o-edgemark means that there is a at least one MAG in the
  equivalence class where the edgemark is a tail, and at least one where
  the edgemark is an arrowhead. Information on the interpretation of edges
  in a MAG can be found in the references given below.

  The first part of the FCI algorithm is analogous to the PC algorithm. It
  starts with a complete undirected graph and estimates an initial skeleton
  using the function \code{\link{skeleton}}. All edges of this skeleton are of
  the form o-o. Due to the presence of hidden variables, it is no longer
  sufficient to consider only subsets of the neighborhoods of nodes \code{x}
  and \code{y} to decide whether the edge \code{x-y} should be removed.
  Therefore, the initial skeleton may contain some superfluous edges.
  These edges are removed in the next step of the algorithm. To decide whether
  edge \code{x o-o y} should be removed, one computes Possible-D-SEP(x) and
  Possible-D-SEP(y) and performs
  conditional indepedence tests of x and y given all possible
  subsets of Possible-D-SEP(x) and
  of Possible-D-SEP(y) (see helpfile of \code{\link{pdsep}}).
  Subsequently, the v-structures are determined (using information in sepset).
  Finally, as many as possible
  undetermined edge marks (o) are determined using (a subset of) the 10
  orientation rules given by Zhang (2009).

  The conservative FCI consists of two parts. In the first part (done if
  the first argument of \code{conservative} is \code{TRUE}), we call
  \code{pc.cons.internal} with option \code{version.unf = c(1,2)} after
  computing the skeleton. This is a slight variation of the conservative
  PC (which used \code{version.unf = c(2,2)}): If a is independent of c
  given some S in the skeleton (i.e., the edge a-c dropped out), but a
  and c remain dependent given all subsets of neighbors of either a or
  c, we will call all triples a-b-c 'faithful'. This is because in the
  FCI, the true separating set might be outside the neighborhood of
  either a or c. In the second part (done if the second argument of
  \code{conservative} is \code{TRUE}), we call \code{pc.cons.internal}
  with option \code{version.unf = c(1,2)} again after Possilbe-D-Sep was
  found and the graph potentially lost some edges. Therefore, new
  triples might have occured. If this second part is done, the resulting
  information on sepset and faithful triples overwrites the previous and
  will be used for the subsequent orientation rules. }

\references{
  T.S. Richardson and P. Spirtes (2002). Ancestral graph Markov models.
  \emph{Annals of Statistics} \bold{30} 962-1030.

  Markus Kalisch, Martin Maechler, Diego Colombo, Marloes H. Maathuis,
  Peter Buehlmann (2012). Causal Inference Using Graphical Models with
  the R Package pcalg. \emph{Journal of Statistical Software}
  \bold{47(11)} 1--26, \url{http://www.jstatsoft.org/v47/i11/}. 
  
  P. Spirtes, C. Glymour and R. Scheines (2000).
  \emph{Causation, Prediction, and Search}, 2nd edition, MIT Press,
  Cambridge (MA).

  P. Spirtes, C. Meek, T.S. Richardson (1999). In: \emph{Computation,
  Causation and Discovery}. An algorithm for causal
  inference in the presence of latent variables and selection bias.
  Pages 211-252. MIT Press.

  J. Zhang (2008). On the completeness of orientation rules for
  causal discovery in the presence of latent confounders and selection bias.
   \emph{Artificial Intelligence} \bold{172} 1873-1896.
}
\seealso{\code{\link{skeleton}} for estimating a skeleton
  using the PC algorithm; \code{\link{pc}} for estimating a CPDAG using
  the PC algorithm; \code{\link{pdsep}} for computing
  Possible-D-Sep for each node and testing and adapting the graph
  accordingly; \code{\link{qreach}} for a fast way of finding possible d-sep
    for a given node; \code{\link{gaussCItest}},
  \code{\link{disCItest}}, \code{\link{binCItest}} and
  \code{\link{dsepTest}} as examples for \code{indepTest}.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo.
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

## find skeleton and CPDAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
indepTest <- gaussCItest
res <- fci(suffStat, indepTest, p, alpha = 0.99, verbose=TRUE)


##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## create the graph
p <- 4
amat1 <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
colnames(amat1) <- rownames(amat1) <- as.character(1:5)
L1 <- 1
V1 <- as.character(1:5)
edL1 <- vector("list",length=5)
names(edL1) <- V1
edL1[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL1[[2]] <- list(edges=3,weights=c(1))
edL1[[3]] <- list(edges=5,weights=c(1))
edL1[[4]] <- list(edges=5,weights=c(1))
g1 <- new("graphNEL", nodes=V1, edgeL=edL1,edgemode="directed")

## compute the true covariance matrix of g1
cov.mat1 <- trueCov(g1)

## delete rows and columns belonging to latent variable L1
true.cov1 <- cov.mat1[-L1,-L1]

## transform covariance matrix into a correlation matrix
true.corr1 <- cov2cor(true.cov1)

## find PAG with FCI algorithm
## as dependence "oracle", we use the true correlation matrix in the
## function gaussCItest with a large "virtual sample size" and a large
## alpha
suffStat1 <- list(C = true.corr1, n = 10^9)
indepTest1 <- gaussCItest
true.pag1 <- fci(suffStat1, indepTest1, p, alpha = 0.99, verbose=TRUE)

## define PAG given in paper
corr.pag1 <- matrix(0,4,4)
corr.pag1[1,] <- c(0,1,1,0)
corr.pag1[2,] <- c(1,0,0,2)
corr.pag1[3,] <- c(1,0,0,2)
corr.pag1[4,] <- c(0,3,3,0)

## check if estimated and correct PAG are in agreement
all(corr.pag1==true.pag1@amat)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

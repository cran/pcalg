\name{udag2apag}
\alias{udag2apag}
\title{Last step of RFCI algorithm: Transform partially oriented graph into RFCI-PAG}
\description{
  This function performs the last step of the RFCI algorithm: It
  transforms a partially oriented graph in which
  the v-structures have been oriented into an RFCI Partial Ancestral
  Graph (PAG) (see Colombo et al (2012)).

  While orienting the edges, this function performs some additional
  conditional independence tests in orientation rule 4 to ensure
  correctness of the ancestral relationships.  As a result of these
  additional tests, some additional edges can be deleted.  The result is
  the final adjacency matrix indicating also the edge marks and the
  updated sepsets.
}
\usage{
udag2apag(apag, suffStat, indepTest, alpha, sepset,
          rules = rep(TRUE, 10), unfVect = NULL, verbose = FALSE)
}
\arguments{
  \item{apag}{Adjacency matrix of the pcAlgo-object of size p*p where
    the v-structures have already been oriented (coding 0,1,2,3 for no
  edge, circle, arrowhead, tail; e.g., \code{amat[a,b] = 2} and
  \code{amat[b,a] = 3} implies a -> b).}
  \item{suffStat}{Sufficient statistics: A \code{\link{list}} containing
    all necessary elements for the conditional independence decisions in the
    function \code{indepTest}.}
  \item{indepTest}{Pre-defined function for testing conditional
    independence.  The function is internally called as
    \code{indepTest(x,y,S,suffStat)}, and tests conditional independence
    of \code{x} and \code{y} given \code{S}.  Here, \code{x} and \code{y}
    are variables, and \code{S} is a (possibly empty) set of variables
    (all variables are coded by their column numbers in the adjacency
    matrix).  \code{suffStat} is a list containing all relevant elements
    for the conditional independence decisions.  The return value of
    \code{indepTest} is the p-value of the test for conditional
    independence.}
  \item{alpha}{Significance level for the individual conditional
    independence tests.}
  \item{sepset}{List of length p; each element of the list
    contains another list of length p.  The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out.  This object is thought to be
    obtained from a \code{pcAlgo}-object.}
  \item{rules}{Logical vector of length 10 with \code{TRUE} or \code{FALSE} for
    each rule, where \code{TRUE} in position i means that rule i (Ri) will be
    applied.  By default, all rules are active.}
  \item{unfVect}{Vector containing numbers that encode the ambiguous
    triples (as returned by \code{\link{pc.cons.intern}}).  This is
    needed in the conservative and in the majority rule versions of RFCI.}
  \item{verbose}{Logical indicating if detailed output is to be given.}
}
\details{
  The partially oriented graph in which the v-structures have been
  oriented is transformed into an RFCI-PAG using adapted rules of Zhang
  (2008).  This function is similar to \code{\link{udag2pag}} used to
  orient the skeleton into a PAG in the FCI algorithm.  However, it is
  slightly more complicated because we perform additional conditional
  independence tests when applying rule 4, to ensure correctness of the
  ancestral relationships.  As a result, some additional edges can be
  deleted, see Colombo et al. (2012).  Because of these addiitonal
  tests, we need to give \code{suffStat}, \code{indepTest}, and
  \code{alpha} as inputs. Since edges can be deleted, the input
  adjacency matrix \code{apag} and the input separating sets
  \code{sepset} can change in this algorithm.

  If \code{unfVect = NULL} (no ambiguous triples), the orientation rules
  are applied to each eligible structure until no more edges can be
  oriented.  On the other hand, hand, if one uses conservative or
  majority rule FCI and ambiguous triples have been found in
  \code{\link{pc.cons.intern}}, \code{unfVect} contains the numbers of
  all ambiguous triples in the graph.  In this case, the orientation
  rules take this information into account.  For example, if a *-> b o-*
  c and <a,b,c> is an unambigous unshielded triple and not a
  v-structure, then we obtain b -* c (otherwise we would create an
  additional v-structure).  On the other hand, if a *-> b o-* c but
  <a,b,c> is an ambiguous unshielded triple, then the circle mark at b
  is not oriented.

  Note that the algorithm works with columns' position of the adjacency
  matrix and not with the names of the variables.

  Note that this function does not resolve possible order-dependence in
  the application of the orientation rules, see Colombo and Maathuis (2013).
}
\value{
  \item{apag}{Final adjacency matrix (coding 0,1,2,3 for no edge,
  circle, arrowhead, tail; e.g., \code{amat[a,b] = 2} and
  \code{amat[b,a] = 3} implies a -> b)}

  \item{sepset}{Updated list of separating sets}
}
\references{
  D. Colombo and M.H. Maathuis (2013).
  Order-independent constraint-based causal structure learning.
  arXiv:1211.3295v2.

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson (2012).
  Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294--321.

  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172}, 1873--1896.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{rfci}}, \code{\link{udag2pag}},
  \code{\link{dag2pag}}, \code{\link{udag2pdag}},
  \code{\link{udag2pdagSpecial}}, \code{\link{udag2pdagRelaxed}}
}
\examples{%% first part is  *identical* to >>>>>>>> ./rfci.Rd and ./udag2pag.Rd
##################################################%  --------       -----------
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## create the DAG :
amat <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
V <- LETTERS[1:5]
colnames(amat) <- rownames(amat) <- V
edL <- setNames(vector("list",length=5), V) %% MM{FIXME}: nicer as in ./jointIda.Rd
edL[[1]] <- list(edges= c(2,4),weights=c(1,1))
edL[[2]] <- list(edges= 3,     weights=c(1))
edL[[3]] <- list(edges= 5,     weights=c(1))
edL[[4]] <- list(edges= 5,     weights=c(1))
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

%% Here, ./rfci.Rd diverts and calls rfci() ...
n <- 100000
alpha <- 0.01
p <- ncol(true.corr)

if (require("MASS")) {
  ## generate 100000 samples of DAG using standard normal error distribution
  set.seed(289)
  d.mat <- mvrnorm(n, mu = rep(0, p), Sigma = true.cov)

  ## estimate the skeleton of given data
  suffStat <- list(C = cor(d.mat), n = n)
  indepTest <- gaussCItest
  resD <- skeleton(suffStat, indepTest, alpha = alpha, labels=colnames(true.corr))

  ## estimate all ordered unshielded triples
  amat.resD <- as(resD@graph, "matrix")
  print(u.t <- find.unsh.triple(amat.resD)) # four of them

  ## check and orient v-structures
  vstrucs <- rfci.vStruc(suffStat, indepTest, alpha=alpha,
			 sepset = resD@sepset, g.amat = amat.resD,
			 unshTripl= u.t$unshTripl, unshVect = u.t$unshVect,
			 verbose = TRUE)

  ## Estimate the final skeleton and extend it into a PAG
  ## (using all 10 rules, as per default):
  resP <- udag2apag(vstrucs$amat, suffStat, indepTest=indepTest, alpha=alpha,
		    sepset=vstrucs$sepset, verbose = TRUE)
  print(Amat <- resP$graph)
} # only if "MASS" is there

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

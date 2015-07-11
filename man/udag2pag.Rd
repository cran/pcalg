\name{udag2pag}
\alias{udag2pag}
\title{Last steps of FCI algorithm: Transform Final Skeleton into FCI-PAG}
\description{
  This function perform the last steps of the FCI algorithm, as it
  transforms an un-oriented final skeleton into a Partial Ancestral
  Graph (PAG).  The final skeleton must have been estimated with
  \code{\link{pdsep}()}.  The result is an adjacency matrix
  indicating also the edge marks.
}
\usage{
udag2pag(pag, sepset, rules = rep(TRUE, 10), unfVect = NULL,
  verbose = FALSE, orientCollider = TRUE)
}
\arguments{
  \item{pag}{Adjacency matrix of the final skeleton of size p*p (coding
    0,1,2,3 for no edge, circle, arrowhead, tail; e.g., \code{amat[a,b]
      = 2} and \code{amat[b,a] = 3} implies a -> b).}
  \item{sepset}{List of length p; each element of the list
    contains another list of length p. The element
    \code{sepset[[x]][[y]]} contains the separation set that made the edge
    between \code{x} and \code{y} drop out. This object is thought to be
    obtained from a \code{pcAlgo}-object.}
  \item{rules}{Array of length 10 containing \code{TRUE} or \code{FALSE}
    for each rule. \code{TRUE} in position i means that rule i (Ri) will
    be applied.  By default, all rules are used.}
  \item{unfVect}{Vector containing numbers that encode ambiguous unshielded
    triples (as returned by \code{\link{pc.cons.intern}}).  This is
    needed in the conservative and majority rule versions of FCI.}
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
  \item{orientCollider}{if \code{TRUE}, collider are oriented.}
}
\details{
  The skeleton is transformed into an FCI-PAG using rules by Zhang (2008).

  If \code{unfVect = NULL} (i.e., one uses standard FCI or one uses
  conservative/majority rule FCI but there are no ambiguous triples),
  then the orientation rules are applied to each eligible structure
  until no more edges can be oriented.  On the other hand, if one uses
  conservative or majority rule FCI and ambiguous triples have been
  found in \code{\link{pc.cons.intern}}, \code{unfVect} contains the
  numbers of all ambiguous triples in the graph.  In this case, the
  orientation rules take this information into account.  For example, if
  a *-> b o-* c and <a,b,c> is an unambigous unshielded triple and not a
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
  Adjacency matrix M with edge marks. The edge marks are
  coded in the following way: M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i]
  != 0: i *-o j; M[i,j]=2, M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.
}
\references{
  D. Colombo and M.H. Maathuis (2013). Order-independent constraint-based
  causal structure learning. arXiv:1211.3295v2

  D. Colombo, M. H. Maathuis, M. Kalisch, T. S. Richardson
  (2012). Learning high-dimensional directed acyclic graphs with latent
  and selection variables. \emph{Ann. Statist.} \bold{40}, 294--321.

  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172}, 1873--1896.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{fci}}, \code{\link{udag2apag}}, \code{\link{dag2pag}};
  further, \code{\link{udag2pdag}} (incl. \code{udag2pdagSpecial} and
  \code{udag2pdagRelaxed}).
}
\examples{
##################################################
## Example with hidden variables
## Zhang (2008), Fig. 6, p.1882
##################################################

## draw a DAG with latent variables
## this example is taken from Zhang (2008), Fig. 6, p.1882 (see references)
amat <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
V <- as.character(1:5)
colnames(amat) <- rownames(amat) <- V
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges= c(2,4),weights=c(1,1))
edL[[2]] <- list(edges= 3,     weights=c(1))
edL[[3]] <- list(edges= 5,     weights=c(1))
edL[[4]] <- list(edges= 5,     weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL,edgemode="directed")

if(require("Rgraphviz"))  plot(g) else print(g)

## define the latent variable
L <- 1

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]

## transform it in a correlation matrix
true.corr <- cov2cor(true.cov)

if (require("MASS")) {
  ## generate 100000 samples of DAG using standard normal error distribution
  n <- 100000
  alpha <- 0.01
  set.seed(314)
  d.mat <- mvrnorm(n, mu = rep(0,dim(true.corr)[1]), Sigma = true.cov)

  ## estimate the skeleton of given data
  suffStat <- list(C = cor(d.mat), n = n)
  indepTest <- gaussCItest
  resD <- skeleton(suffStat, indepTest, p=dim(true.corr)[2], alpha = alpha)

  ## estimate v-structures conservatively
  tmp <- pc.cons.intern(resD, suffStat, indepTest, alpha, version.unf = c(1, 1))
  ## tripleList <- tmp$unfTripl
  resD <- tmp$sk

  ## estimate the final skeleton of given data using Possible-D-Sep
  pdsepRes <- pdsep(resD@graph, suffStat, indepTest, p=dim(true.corr)[2],
		    resD@sepset, alpha = alpha, m.max = Inf,
		    pMax = resD@pMax)

  ## extend the skeleton into a PAG using all 10 rules
  resP <- udag2pag(pag = pdsepRes$G, pdsepRes$sepset, rules = rep(TRUE,10),
		   verbose = TRUE)
  colnames(resP) <- rownames(resP) <- as.character(2:5)
  print(resP)

} # only if "MASS" is there

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



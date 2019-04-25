\name{backdoor}
\alias{backdoor}
\title{Find Set Satisfying the Generalized Backdoor Criterion (GBC)}
\description{
  This function first checks if the total causal effect of
  one variable (\code{x}) onto another variable (\code{y}) is
  identifiable via the GBC, and if this is
  the case it explicitly gives a set of variables that satisfies the
  GBC with respect to \code{x} and \code{y}
  in the given graph.
}
\usage{
backdoor(amat, x, y, type = "pag", max.chordal = 10, verbose=FALSE)
}
\arguments{
  \item{amat}{adjacency matrix of type \code{\link{amat.cpdag}} or
    \code{\link{amat.pag}}.}
  \item{x,y}{(integer) position of variable \eqn{X} and \eqn{Y},
    respectively, in the adjacency matrix.}
  \item{type}{string specifying the type of graph of the adjacency matrix
    \code{amat}.  It can be a DAG (type="dag"), a CPDAG (type="cpdag");
    then the type of the adjacency matrix is assumed to be
    \link{amat.cpdag}.  It can also be a MAG (type="mag"), or a PAG
    (type="pag"); then the type of the adjacency matrix is assumed to be
  \link{amat.pag}.}
  \item{max.chordal}{only if \code{type = "mag"}, is used in
    \code{\link{pag2magAM}} to determine paths too large to be checked
    for chordality.}
  \item{verbose}{logical; if true, some output is produced during
    computation.}
}
\details{
  This function is a generalization of Pearl's backdoor criterion, see
  Pearl (1993), defined for directed acyclic graphs (DAGs), for single
  interventions and single outcome variable to more general types of
  graphs (CPDAGs, MAGs, and PAGs) that describe Markov equivalence
  classes of DAGs with and without latent variables but without
  selection variables. For more details see Maathuis and Colombo (2015).

  The motivation to find a set W that satisfies the GBC with respect to
  \code{x} and \code{y}
  in the given graph relies on the result of the generalized backdoor adjustment:

  \emph{If a set of variables W satisfies the GBC relative to \code{x}
    and \code{y} in the given graph, then
    the causal effect of \code{x} on \code{y} is identifiable and is given
    by} \deqn{%
    P(Y|do(X = x)) = \sum_W P(Y|X,W) \cdot P(W).}{%
    P(Y|do(X = x)) =  sum_W P(Y|X,W) * P(W).}

  This result allows to write post-intervention densities (the one
  written using Pearl's do-calculus) using only observational densities
  estimated from the data.

  If the input graph is a DAG (\code{type="dag"}), this function reduces
  to Pearl's backdoor criterion for single interventions and single
  outcome variable, and the parents of \code{x} in the DAG satisfy the
  backdoor criterion unless \code{y} is a parent of \code{x}.

  If the input graph is a CPDAG C (\code{type="cpdag"}), a MAG M
  (\code{type="mag"}), or a PAG P (\code{type="pag"}) (with both M and P
  not allowing selection variables), this function first checks if the
  total causal effect of \code{x} on \code{y} is identifiable via the
  GBC (see Maathuis and Colombo, 2015). If
  the effect is not identifiable in this way, the output is
  NA. Otherwise, an explicit set W that satisfies the GBC with respect
  to \code{x} and \code{y} in the given graph is found.

  At this moment this function is not able to work with an RFCI-PAG.

  It is important to note that there can be pair of nodes \code{x} and
  \code{y} for which there is no set W that satisfies the GBC, but the
  total causal effect might be identifiable via some other technique.

  For the coding of the adjacency matrix see \link{amatType}.
}
\value{
  Either NA if the total causal effect is not identifiable via the
  GBC, or a set if the effect is identifiable
  via the GBC. Note that if the set W is
  equal to the empty set, the output is NULL.
}
\references{
  M.H. Maathuis and D. Colombo (2015). A generalized backdoor
  criterion. Annals of Statistics 43 1060-1088.

  J. Pearl (1993). Comment: Graphical models, causality and intervention.
  \emph{Statistical Science} \bold{8}, 266--269.
}
\author{Diego Colombo and Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{gac}} for the Generalized Adjustment Criterion
  (GAC), which is a generalization of GBC; \code{\link{pc}} for
  estimating a CPDAG, \code{\link{dag2pag}}
  and \code{\link{fci}} for estimating a PAG, and
  \code{\link{pag2magAM}} for estimating a MAG.
}
\examples{%% note:  Tests in  ../tests/test_backdoor.R
#####################################################################
##DAG
#####################################################################
## Simulate the true DAG
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG

## Extract the adjacency matrix of the true DAG
true.amat <- (amat <- as(myDAG, "matrix")) != 0 # TRUE/FALSE <==> 1/0
print.table(1*true.amat, zero.=".") # "visualization"

## Compute set satisfying the GBC:
backdoor(true.amat, 5, 7, type="dag")
\dontshow{stopifnot(backdoor(true.amat, 5, 7, type="dag") == 3)}
#####################################################################
##CPDAG
#####################################################################
##################################################
## Example not identifiable
## Maathuis and Colombo (2015), Fig. 3a, p.1072
##################################################
## create the graph
p <- 5
. <- 0
amat <- rbind(c(.,.,1,1,1),
              c(.,.,1,1,1),
              c(.,.,.,1,.),
              c(.,.,.,.,1),
              c(.,.,.,.,.))
colnames(amat) <- rownames(amat) <- as.character(1:5)
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,5),weights=c(1,1,1))
edL[[3]] <- list(edges=4,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## estimate the true CPDAG
myCPDAG <- dag2cpdag(g)
## Extract the adjacency matrix of the true CPDAG
true.amat <- (as(myCPDAG, "matrix") != 0) # 1/0 <==> TRUE/FALSE

## The effect is not identifiable, in fact:
backdoor(true.amat, 3, 5, type="cpdag")
\dontshow{stopifnot(identical(NA, backdoor(true.amat, 3, 5, type="cpdag")))}

##################################################
## Example identifiable
## Maathuis and Colombo (2015), Fig. 3b, p.1072
##################################################

## create the graph
p <- 6
amat <- rbind(c(0,0,1,1,0,1), c(0,0,1,1,0,1), c(0,0,0,0,1,0),
              c(0,0,0,0,1,1), c(0,0,0,0,0,0), c(0,0,0,0,0,0))
colnames(amat) <- rownames(amat) <- as.character(1:6)
V <- as.character(1:6)
edL <- vector("list",length=6)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=c(5,6),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

## estimate the true CPDAG
myCPDAG <- dag2cpdag(g)
## Extract the adjacency matrix of the true CPDAG
true.amat <- as(myCPDAG, "matrix") != 0

## The effect is identifiable and the set satisfying GBC is:
backdoor(true.amat, 6, 3, type="cpdag")
\dontshow{stopifnot(backdoor(true.amat, 6, 3, type="cpdag") == 1:2)}

##################################################################
##PAG
##################################################################
##################################################
## Example identifiable
## Maathuis and Colombo (2015), Fig. 5a, p.1075
##################################################

## create the graph
p <- 7
amat <- t(matrix(c(0,0,1,1,0,0,0, 0,0,1,1,0,0,0, 0,0,0,1,0,1,0,
                   0,0,0,0,0,0,1, 0,0,0,0,0,1,1, 0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0),  7, 7))
colnames(amat) <- rownames(amat) <- as.character(1:7)
V <- as.character(1:7)
edL <- vector("list",length=7)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4),weights=c(1,1))
edL[[2]] <- list(edges=c(3,4),weights=c(1,1))
edL[[3]] <- list(edges=c(4,6),weights=c(1,1))
edL[[4]] <- list(edges=7,weights=c(1))
edL[[5]] <- list(edges=c(6,7),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
L <- 5

## compute the true covariance matrix of g
cov.mat <- trueCov(g)

## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)
suffStat <- list(C=true.corr, n=10^9)
indepTest <- gaussCItest

## estimate the true PAG
true.pag <- dag2pag(suffStat, indepTest, g, L, alpha = 0.9999)@amat

## The effect is identifiable  and the backdoor set is:
backdoor(true.pag, 3, 5, type="pag")
\dontshow{stopifnot(backdoor(true.pag, 3, 5, type="pag") == 1:2)}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

\name{udag2pag}
\alias{udag2pag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extend a pcAlgo-object containing a skeleton to a PAG}
\description{
This function extends a pcAlgo-object containing a skeleton and
corresponding conditional independence information to a Partial
Ancestral Graph (PAG). The pcAlgo-object must have been estimated with
the option \code{psepset=TRUE}. The result is an adjacency matrix
indicating also the edge marks.
}
\usage{
udag2pag(gInput, rules, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gInput}{pcAlgo-object containing skeleton and
    cond. ind. information; must have been estimated with
    \code{psepset=TRUE}}
  \item{rules}{array of length 10 containing TRUE or FALSE for each
    rule. TRUE in position i means that rule i (Ri) will be used. Per
    default, all rules are set to TRUE.}
  \item{verbose}{0: No output; 1: Details. Default is 0.}
}
\details{
  The skeleton is extended to a PAG using rules by Zhang (see
  References). Note that the algorithm works with columns' position of the
  adjacency matrix and not with the names of the variables. 
}
\value{
  The output is an adjacency matrix M with edge marks. The edge marks are
  coded in the following way: M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i]
  != 0: i *-o j; M[i,j]=2, M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j. 
}
\references{
  J. Zhang (2008). On the completeness of orientation rules for causal
  discovery in the presence of latent confounders and selection bias.
  \emph{Artificial Intelligence} \bold{172} 1873-1896.
}
\author{Diego Colombo (\email{colombo@stat.math.ethz.ch})}

\seealso{\code{\link{pc}}}

\examples{
## generate and draw random DAG
## this example is taken from Zhang (2008), Fig. 6, p.1882 (see references)
amat <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
colnames(amat) <- rownames(amat) <- as.character(1:5)
V <- as.character(1:5)
edL <- vector("list",length=5)
names(edL) <- V
edL[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL[[2]] <- list(edges=3,weights=c(1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=5,weights=c(1))
g <- new("graphNEL", nodes=V, edgeL=edL,edgemode="directed")
if (require(Rgraphviz)) {
plot(g)
}

## define the latent variable
L <- 1

## generate 100000 samples of DAG using standard normal error distribution
n <- 100000
d.mat <- rmvDAG(n, g, errDist = "normal")

## delete rows and columns corresponding to the latent variable
d.mat <- d.mat[-L,-L]

## estimate the skeleton of given data using psepset=TRUE
resD <- pcAlgo(d.mat, alpha=0.05, psepset=TRUE)
resD

## extend the pcalgo-object into a PAG using all 10 rules
rules <- rep(TRUE,10)
resP <- udag2pag(resD, rules=rules, verbose=1)
colnames(resP) <- rownames(resP) <- nodes(g)[-L]

if (require(Rgraphviz)) {
## plot the original DAG and the PAG
op <- par(mfrow=c(1,2))
plot(g, main="original (random) DAG")
plotAG(resP)
par(op)
}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



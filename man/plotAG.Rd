\name{plotAG}
\alias{plotAG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot partial ancestral graphs (PAG)}
\description{
This function plots partial ancestral graphs, namely graphs with the
following six kinds of edges: -, ->, <->, o-, o-o, o->. 
}
\usage{
plotAG(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{adjacency matrix M with edge marks. The edge marks are
  coded in the following way: M[i,j]=M[j,i]=0: no edge; M[i,j]=1, M[j,i]
  != 0: i *-o j; M[i,j]=2, M[j,i]!=0: i*->j; M[i,j]=3, M[j,i]!=0: i*-j.}
}
\details{
  Used to plot the output of the function udag2pag. Note that Rgraphviz
  places the edge marks sometimes not quite at the end of the edge. This
  is just a layout problem and should not hinder a correct interpretation.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}

\seealso{\code{\link{pcAlgo}}}

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
plot(g)

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

## plot the original DAG and the PAG
op <- par(mfrow=c(1,2))
plot(g, main="original (random) DAG")
plotAG(resP)
par(op)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



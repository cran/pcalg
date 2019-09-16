\name{optAdjSet}
\alias{optAdjSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the optimal adjustment set}
\description{
\code{optAdjSet} computes the optimal valid adjustment set relative to the variables (\code{X},\code{Y}) in the given graph.
}
\usage{
optAdjSet(graphEst,x.pos,y.pos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{graphEst}{graphNel object or adjacency matrix of type amat.cpdag.}
  \item{x.pos, x}{Positions of variables \code{X} in the covariance matrix.}
  \item{y.pos, y}{Positions of variables \code{Y} in the covariance matrix.}
}
\details{
Suppose we have data from a linear SEM compatible with a known causal graph \code{G} and our aim is to estimate the total joint effect of \code{X} on \code{Y}. Here the total joint effect of \code{X} \eqn{= (X_1,X_2)} on \code{Y} is defined via Pearl's do-calculus as the vector \eqn{(E[Y|do(X_1=x_1+1,X_2=x_2)]-E[Y|do(X_1=x_1,X_2=x_2)], E[Y|do(X_1=x_1,X_2=x_2+1)]-E[Y|do(X_1=x_1,X_2=x_2)])}, with a similar definition for more than two variables. These values are equal to the partial derivatives (evaluated at \eqn{x_1,x_2}) of \eqn{E[Y|do(X=x_1',X_2=x_2')]} with respect to \eqn{x_1}' and \eqn{x_2}'. Moreover, under the linearity assumption, these partial derivatives do not depend on the values at which they are evaluated.

It is possible to estimate the total joint effect of \code{X} on \code{Y} with a simple linear regression of the form \code{lm(Y ~ X + Z)}, if and only if the covariate set \code{Z} is a valid adjustment set (see Perkovic et al. (2018)). Often, however, there are multiple such valid adjustment sets, providing total effect estimates with varying accuracies. Suppose that there exists a valid adjustment set relative to (\code{X},\code{Y}) in causal graph \code{G}, and each node in Y is a descendant of X, then there exists a valid adjustment which provides the total effect estimate with the optimal asymptotic variance, which we will refer to as \code{O(X,Y,G)} (Henckel et al., 2019). This function returns this optimal valid adjustment set \code{O(X,Y,G)}. 

The restriction that each node in \code{Y} be a descendant of the node set \code{X} is not notable, as the total effect of the node set \code{X} on a non-descendant is always 0. If provided with a node set \code{Y} that does not fulfill this condition this function computes a pruned node set \code{Y2} by removing all nodes from \code{Y} that are not descendants of \code{X} and returns \code{O(X,Y2,G)} instead. The user will be alerted to this and given the pruned set \code{Y2}.


}
\value{
A vector with the positions of the nodes of the optimal set O(\code{X},\code{Y},\code{G}).
}
\references{
  E. Perković, J. Textor, M. Kalisch and M.H. Maathuis (2018). Complete graphical characterization and construction of adjustment   sets in Markov equivalence classes of ancestral graphs. \emph{Journal of Machine Learning Research.} \bold{18}(220) 1--62,
    
  L. Henckel, E. Perkovic and M.H. Maathuis (2019). Graphical criteria for efficient total 
  effect estimation via adjustment in causal linear models.
  \emph{Working Paper.}
    
}
\author{
Leonard Henckel
}

\examples{
## Simulate a true DAG, its CPDAG and an intermediate max. PDAG
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
p <- 10
## true DAG
myDAG <- randomDAG(p, prob = 0.3) 
## true CPDAG
myCPDAG <- dag2cpdag(myDAG) 
## true PDAG with added background knowledge 5 -> 6
myPDAG <- addBgKnowledge(myCPDAG,5,6) 
par(mfrow = c(1,3))
plot(myDAG)
plot(myPDAG)
plot(myCPDAG) ## plot of the graphs

## if the CPDAG C is amenable relative to (X,Y),
## the optimal set will be the same for all DAGs 
## and any max. PDAGs obtained by adding background knowledge to C 
(optAdjSet(myDAG,3,10))
(optAdjSet(myPDAG,3,10))
(optAdjSet(myCPDAG,3,10))


## the optimal adjustment set can also be compute for sets X and Y
(optAdjSet(myDAG,c(3,4),c(9,10)))
(optAdjSet(myPDAG,c(3,4),c(9,10)))
(optAdjSet(myCPDAG,c(3,4),c(9,10)))

## The only restriction is that it requires all nodes in Y to be
## descendants of X.
## However, if a node in Y is non-descendant of X the lowest variance
## partial total effect estimate is simply 0.
## Hence, we can proceed with a pruned Y. This function does this automatically!
optAdjSet(myDAG,1,c(3,9))

## Note that for sets X there may be no valid adjustment set even
## if the PDAG is is amenable relative to (X,Y).
\dontrun{optAdjSet(myPDAG,c(4,9),7)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}% __ONLY ONE__ keyword per line

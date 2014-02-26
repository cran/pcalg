\name{pag2mag}
\alias{pag2magAM}
\title{Transform a PAG into a MAG in the Corresponding Markov Equivalence Class}
\description{
  Transform a Partial Ancestral Graph (PAG) into a valid
  Maximal Ancestral Graph (MAG) that belongs to the Markov equivalence class
  represented by the given PAG, with no additional edges into node \code{x}.
}
\usage{
pag2magAM(amat.pag, x, max.chordal = 10, verbose = FALSE)
}
\arguments{
  \item{amat.pag}{Adjacency matrix (coding 0,1,2,3 for no edge, circle,
    arrowhead, tail; e.g., \code{amat[a,b] = 2} and \code{amat[b,a] = 3}
    implies a -> b)} 
  \item{x}{Node in the PAG into which no additional edges are oriented.}
  \item{max.chordal}{Positive integer: graph paths larger than
    \code{max.chordal} are considered to be too large to be checked for
    chordality.}
  \item{verbose}{Logical; if true, some output is produced during
    computation.}  
}
\details{
  This function converts a PAG (adjacency matrix) to a valid MAG
  (adjacency matrix) that belongs to the Markov equivalence class
  represented by the given PAG.  Note that we assume that there are no
  selection variables, meaning that the edges in the PAG can be of the
  following types: ->, <->, o->, and o-o.  In a first step, it uses the
  Arrowhead Augmentation of Zhang (2006), i.e., any o-> edge is oriented
  into ->.  Afterwards, it orients each chordal component into a valid DAG
  without orienting any additional edges into x.

  This function is used in the Generalized Backdoor Criterion
  \code{\link{backdoor}} with \code{type="pag"}, see Maathuis and Colombo
  (2013) for details.
}
\value{
  The output is an adjacency matrix M of a valid MAG with edge
  marks.  The edge marks are coded in the following way: M[i,j]=M[j,i]=0:
  no edge; M[i,j]=2, M[j,i]=3: i -> j; M[i,j]=2, M[j,i]=2: i <-> j.
}
\references{
  M.H. Maathuis and D. Colombo (2013). A generalized backdoor
  criterion. arXiv preprint arXiv:1307.5636.

  Zhang, J. (2006). Causal Inference and Reasoning in Causally
  Insufficient Systems. Ph. D. thesis, Carnegie Mellon University.
}
\author{
  Diego Colombo, Markus Kalisch and Martin Maechler.
}
\seealso{\code{\link{fci}}, \code{\link{dag2pag}}, \code{\link{backdoor}}}
\examples{
## create the graph
set.seed(78)
p <- 12
g <- randomDAG(p, prob = 0.4)
## Compute the true covariance and then correlation matrix of g:
true.corr <- cov2cor(trueCov(g))

## define nodes 2 and 6 to be latent variables
L <- c(2,6)

## Find PAG
## As dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
true.pag <- dag2pag(suffStat = list(C= true.corr, n= 10^9),
                    indepTest= gaussCItest, graph=g, L=L, alpha= 0.9999)

## find a valid MAG such that no additional edges are directed into
(amat.mag <- pag2magAM(true.pag@amat, 4)) # -> the adj.matrix of the MAG
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

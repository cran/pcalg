\name{idaFast}
\alias{idaFast}
\title{Estimate multiset of possible total causal effects for 
  several target variables at the same time.}
\description{
  This function estimates the multiset of possible total causal effects of 
  one variable (\code{x}) on a vector of target variables (\code{y}) from
  observational data.
  This method is more efficient than looping over \code{\link{ida}}.  
  Only method="local" (see help file of \code{\link{ida}}) is available.
}
\usage{
idaFast(x.pos,y.pos.set,mcov,graphEst)
}
\arguments{
\item{x.pos}{Position of variable \code{x} in the covariance matrix}
\item{y.pos.set}{Vector containing the position of the target variables
  \code{y} in the covariance matrix}
\item{mcov}{Covariance matrix that was used to estimate \code{graphEst}}
\item{graphEst}{Estimated CPDAG from the function
  \code{\link{pc}}. If the output of \code{\link{pc}} is \code{pc.fit},
  then the estimated CPDAG can be obtained by \code{pc.fit@graph}.}
}
\details{
  This function performs
  \code{ida(x.pos, y.pos, mcov, graphEst, method="local", y.notparent=FALSE, 
  verbose=FALSE)} for
  all values of \code{y.pos} in \code{y.pos.set} at the same time, 
  in an efficient way. 
  See the help file of \code{\link{ida}} for more details. Note that the
  option \code{y.notparent = TRUE} is not implemented, since it is not
  clear how to do that efficiently without orienting all edges away from
  \code{y.pos.set} at the same time, which seems not to be
  desirable. Suggestions are welcome.
}
\value{
  Matrix with \code{length(y.pos.set)} rows. Row i contains the multiset
  of estimated possible total causal effects of \code{x} on
  \code{y.pos.set[i]}. Note that all multisets in the matrix have the
  same length, since the parents of \code{x} are the same for all elements
  of \code{y.pos.set}.
}
\references{
  M.H. Maathuis, M. Kalisch, P. Buehlmann (2009). Estimating
    high-dimensional intervention effects from observational data. 
    \emph{Annals of Statistics}. \bold{37}. 3122 - 3164.

    M.H. Maathuis, D. Colombo, M. Kalisch,
  P. Buehlmann (2010). 
  Predicting causal effects in large-scale systems from observational data.
  \emph{Nature Methods, to appear}. 
}

\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pc}} for estimating a CPDAG, \code{\link{ida}}
  estimating the multiset of possible total causal effects from
  observational data on only one goal variable but with many more
  options than in \code{\link{idaFast}}.
}

\examples{
## Simulate the true DAG
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

## simulate data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)

## estimate CPDAG (see help on the function "pc")
alpha <- 0.01
indepTest <- gaussCItest 
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest, p, alpha)

if (require(Rgraphviz)) {
plot(myDAG)
}

(eff.est1 <- ida(2,5,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
(eff.est2 <- ida(2,6,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
(eff.est3 <- ida(2,7,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
## These three computations can be combinded in an efficient way by using idaFast
(eff.estF <- idaFast(2,c(5,6,7),cov(dat),pc.fit@graph))

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



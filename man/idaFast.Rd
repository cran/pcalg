\name{idaFast}
\alias{idaFast}
\title{Estimate multiset of possible total causal effects for 
  several target nodes at the same time.}
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
\item{x.pos}{Column index of \code{x} in the covariance matrix}
\item{y.pos.set}{Vector of column indices of the target variables \code{y} in 
the covariance matrix}
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
  See the help file of \code{\link{ida}} for more details.
}
\value{
  Matrix with \code{length(y.pos.set)} rows. Row i contains the multiset
  of estimated possible total causal effects of \code{x} on
  \code{y.pos.set[i]}. Note that all multisets in the matrix have the
  same length, since the parents of \code{x} are the same for elements
  of \code{y.pos.set}.
}
\references{
  M.H. Maathuis, M. Kalisch, P. B\"uhlmann (2009), \emph{Estimating
    high-dimensional intervention effects from observational data};
  Annals of Statistics, 2009.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pc}},\code{\link{idaFast}}}
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

plot(myDAG)

(eff.est1 <- ida(2,5,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
(eff.est2 <- ida(2,6,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
(eff.est3 <- ida(2,7,cov(dat),pc.fit@graph,method="local",verbose=FALSE))
## These three computations can be combinded in an efficient way by using idaFast
(eff.estF <- idaFast(2,c(5,6,7),cov(dat),pc.fit@graph))

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



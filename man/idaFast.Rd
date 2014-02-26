\name{idaFast}
\alias{idaFast}
\title{Multiset of Possible Total Causal Effects for Several Target Var.s}
\description{
  This function estimates the multiset of possible total causal effects of
  one variable (\code{x}) on a \emph{several} (i.e., a vector of) target
  variables (\code{y}) from observational data.

  \code{idaFast()} is more efficient than looping over
  \code{\link{ida}}.  Only \code{method="local"} (see \code{\link{ida}})
  is available.
}
\usage{
idaFast(x.pos, y.pos.set, mcov, graphEst)
}
\arguments{
\item{x.pos}{integer position of variable \code{x} in the covariance matrix.}
\item{y.pos.set}{integer vector of positions of the target variables
  \code{y} in the covariance matrix.}
\item{mcov}{covariance matrix that was used to estimate \code{graphEst}}
\item{graphEst}{estimated CPDAG from the function
  \code{\link{pc}}.  If the output of \code{\link{pc}} is \code{pc.fit},
  then the estimated CPDAG can be obtained by \code{pc.fit@graph}.}
}
\details{
  This function performs
  \code{\link{ida}(x.pos, y.pos, mcov, graphEst, method="local",
    y.notparent=FALSE, verbose=FALSE)} for all values of \code{y.pos} in
  \code{y.pos.set} simultaneously, in an efficient way.
  See (the help about) \code{\link{ida}} for more details.  Note that the
  option \code{y.notparent = TRUE} is not implemented, since it is not
  clear how to do that efficiently without orienting all edges away from
  \code{y.pos.set} at the same time, which seems not to be
  desirable.  Suggestions are welcome.
}
\value{
  Matrix with \code{length(y.pos.set)} rows.  Row \eqn{i} contains the multiset
  of estimated possible total causal effects of \code{x} on
  \code{y.pos.set[i]}.  Note that all multisets in the matrix have the
  same length, since the parents of \code{x} are the same for all elements
  of \code{y.pos.set}.
}
\references{
  see the list in \code{\link{ida}}.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pc}} for estimating a CPDAG, and
  \code{\link{ida}} for estimating the multiset of possible total causal
  effects from observational data on only one target variable but with many more
  options (than here in \code{idaFast}).
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
cov.d <- cov(dat)

## estimate CPDAG (see help on the function "pc")
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)

if(require(Rgraphviz)) {
  op <- par(mfrow=c(1,3))
  plot(myDAG,        main="true DAG")
  plot(myCPDAG,      main="true CPDAG")
  plot(pc.fit@graph, main="pc()-estimated CPDAG")
  par(op)
}

(eff.est1 <- ida(2,5, cov.d, pc.fit@graph))## method = "local" is default
(eff.est2 <- ida(2,6, cov.d, pc.fit@graph))
(eff.est3 <- ida(2,7, cov.d, pc.fit@graph))
## These three computations can be combinded in an efficient way
## by using idaFast :
(eff.estF <- idaFast(2, c(5,6,7), cov.d, pc.fit@graph))
\dontshow{
stopifnot(all.equal(eff.estF, tol = 1e-15,
                    rbind("5"= eff.est1, "6"= eff.est2, "7"= eff.est3)) ,
          all.equal(idaFast(2, c(7,5), covTrue, myCPDAG),
                    rbind("7"= ida(2, 7, covTrue, myCPDAG),
                          "5"= ida(2, 5, covTrue, myCPDAG)), tol = 1e-15)
 )
}%dont
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

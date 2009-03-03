\name{beta.special.pcObj}
\alias{beta.special.pcObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute set of intervention effects in a fast way}
\description{
  This function computes a set of intervention effects of one variable
  onto another variable in a fast way. The pc-object has to be precomputed.
}
\usage{
beta.special.pcObj(x.pos,y.pos,pcObj,mcov=NA,amat=NA,amatSkel=NA,t.amat=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x.pos}{Column of x in dat}
\item{y.pos}{Column of y in dat}
\item{pcObj}{Precomputed pc-object}
\item{mcov}{covariance that was used in the pc-object fit}
\item{amat, amatSkel, t.amat}{matrices that can be precomputed, if
  needed (see code for details on how to precompute)}
}
\details{
Estimate the intervention effect of x on y; the pcObj has to be
precomputed. This method is intended to be a fast version of

beta.special(dat=NA,x.pos,y.pos,verbose=0,a=NA,myDAG=NA,myplot=FALSE,
perfect=FALSE,method="local",collTest=TRUE,pcObj=pcObj,all.dags=NA,
scaled.data=FALSE,u2pd="relaxed")

Thus, this is a faster version for the local method given a
precomputed PC-Algo Object (relaxed udag2pdag, so CPDAG might not
be a real CPDAG; this does not matter, since we try not to extend).

If this function is used repeatedly on the same pc-object (e.g. for
different values of x and y), the matrices amat, amatSkel and t.amat can
be precomputed in order to be more efficient (see the first entry in the
code for furhter details).
}
\value{
  estimates of intervention effects
}
\references{
  M.H. Maathuis, M. Kalisch, P. B\"uhlmann (2009), \emph{Estimating
    high-dimensional intervention effects from observational data};
  Annals of Statistics, 2009.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pcAlgo}},\code{\link{dag2cpdag}},\code{\link{beta.special.pcObj}}}
\examples{
set.seed(125)
p <- 10
n <- 10000

myDAG <- randomDAG(p, prob = 0.4)

cov.t <- trueCov(myDAG)
pcObj <- pcAlgo.Perfect(cov2cor(cov.t),directed=TRUE)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
cov.e <- cov(d.mat)
pcObj.est <- pcAlgo(d.mat,alpha=0.01,directed=TRUE)

x.pos <- 1
y.pos <- 5
true.eff <- causalEffect(myDAG,y.pos,x.pos)
est.eff.p <- beta.special.pcObj(x.pos=x.pos,y.pos=y.pos,pcObj=pcObj,mcov=cov.t)
est.eff.e <- beta.special.pcObj(x.pos=x.pos,y.pos=y.pos,pcObj=pcObj.est,mcov=cov.e)

true.eff
est.eff.p
est.eff.e

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



\name{beta.special}
\alias{beta.special}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute set of intervention effects}
\description{
  This function computes a set of intervention effects of one variable
  onto another variable. 
}
\usage{
beta.special(dat=NA,x.pos,y.pos,verbose=0,a=0.01,myDAG=NA,myplot=FALSE,perfect=FALSE,method="local",collTest=TRUE,pcObj=NA,all.dags=NA,u2pd="rand")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{data matrix}
\item{x.pos}{Column of x in dat}
\item{y.pos}{Column of y in dat}
\item{verbose}{0=no comments, 2=detail on estimates}
\item{a}{significance level of tests for finding CPDAG}
\item{myDAG}{needed if true correlation matrix shall be computed}
\item{myplot}{plot estimated graph}
\item{perfect}{True cor matrix is calculated from myDAG}
\item{method}{"local" - local (all combinations of parents in regr.);
  "global" - all DAGs}
\item{collTest}{True - Exclude orientations of undirected edges that
  introduce a new collider}
\item{pcObj}{Fit of PC Algorithm (CPDAG); if this is available, no
  new fit is done}
\item{all.dags}{All DAGs in the format of function allDags; if this is
  available, no new function call allDags is done}
\item{u2pd}{Function for converting udag to pdag;"rand":
  udag2pdag;"relaxed": udag2pdagRelaxed;"retry": udag2pdagSpecial} 
}
\details{
  After estimating a CPDAG, intervention effects of x on y are
  computed. For details, see references.
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
##############################
## Global
##############################
set.seed(125)
p <- 10
n <- 100000

myDAG <- randomDAG(p, prob = 0.4)
plot(dag2cpdag(myDAG))
cmat <- cov2cor(trueCov(myDAG))
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

x.pos <- 1
y.pos <- 5
## true value
true.eff <- causalEffect(myDAG,y.pos,x.pos)
## value using oracle (true value should be element of this)
est.eff.p <- beta.special(x.pos=x.pos,y.pos=y.pos,myDAG=myDAG,perfect=TRUE,method="global")
## value using data (true value should come close to some element of this)
est.eff.e <- beta.special(dat=d.mat,x.pos=x.pos,y.pos=y.pos,perfect=FALSE,method="global",myplot=TRUE)

true.eff
est.eff.p
est.eff.e

## Setting perfect&global always contains exact value of true value
## OK

##############################
## Local
##############################
set.seed(125)
p <- 10
n <- 100000

myDAG <- randomDAG(p, prob = 0.4)
cmat <- cov2cor(trueCov(myDAG))
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

x.pos <- 1
y.pos <- 5
## true value
true.eff <- causalEffect(myDAG,y.pos,x.pos)
## value using oracle (true value should be element of this)
est.eff.p <- beta.special(x.pos=x.pos,y.pos=y.pos,myDAG=myDAG,perfect=TRUE,method="local")
## value using data (true value should come close to some element of this)
est.eff.e <- beta.special(dat=d.mat,x.pos=x.pos,y.pos=y.pos,perfect=FALSE,method="local")
plot(pcAlgo(d.mat,alpha=0.01,directed=TRUE))

true.eff
est.eff.p
est.eff.e

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



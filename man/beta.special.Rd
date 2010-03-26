\name{beta.special}
\alias{beta.special}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute set of intervention effects}
\description{
  This function is DEPRECATED! Use \code{\link{ida}} instead.
}

\usage{ beta.special(dat=NA, x.pos, y.pos, verbose=0, a=0.01, myDAG=NA,
myplot=FALSE, perfect=FALSE, method="local", collTest=TRUE, pcObj=NA,
all.dags=NA, u2pd="rand") }

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
\value{
  estimates of intervention effects
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



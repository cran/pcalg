\name{beta.special}
\alias{beta.special}
\title{Compute set of intervention effects}
\description{
  This function is DEPRECATED! Use \code{\link{ida}} instead.
}
\usage{
beta.special(dat=NA, x.pos, y.pos, verbose=0, a=0.01, myDAG=NA,
             myplot=FALSE, perfect=FALSE, method="local", collTest=TRUE,
             pcObj=NA, all.dags=NA, u2pd="rand")
}
\arguments{
\item{dat}{Data matrix}
\item{x.pos, y.pos}{integer column positions of \eqn{x} and \eqn{y} in \code{dat}.}
\item{verbose}{0=no comments, 2=detail on estimates}
\item{a}{Significance level of tests for finding CPDAG}
\item{myDAG}{Needed if true correlation matrix shall be computed}
\item{myplot}{Plot estimated graph}
\item{perfect}{True cor matrix is calculated from myDAG}
\item{method}{"local" - local (all combinations of parents in regr.);
  "global" - all DAGs}
\item{collTest}{True - Exclude orientations of undirected edges that
  introduce a new collider}
\item{pcObj}{Fit of PC Algorithm (CPDAG); if this is available, no
  new fit is done}
\item{all.dags}{All DAGs in the format of function allDags; if this is
  available, no new function call allDags is done}
\item{u2pd}{function for converting a  UDAG to a PDAG;
  "rand": \code{\link{udag2pdag}};
  "relaxed": \code{\link{udag2pdagRelaxed}};
  "retry": \code{\link{udag2pdagSpecial}}.}
}
\value{
  estimates of intervention effects
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{\code{\link{pcAlgo}}, \code{\link{dag2cpdag}};
  \code{\link{beta.special.pcObj}} for a \emph{fast} version of
  \code{beta.special()}, using a precomputed pc-object.
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}



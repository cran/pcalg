\name{pcSelect}
\alias{pcSelect}
\title{PC-Select: Estimate subgraph around a response variable}
\description{
  This function is intended for feature selection: If you have a
  response variable \eqn{y} and a data matrix \eqn{dm}, which columns
  are "strongly influential" on \eqn{y}. The type of influence is the
  same as in the PC-Algorithm, i.e., \eqn{y} and \eqn{x} (a column of
  \eqn{dm}) are associated if they are correlated even when conditioning
  on any subset of the remaining columns in \eqn{dm}. Therefore, only
  very strong relations will be found and the result is typically a
  subset of other feature selection techniques. Note that there are also
  robust correlation methods available which render this method robust.
}
\usage{
pcSelect(y,dm, alpha, corMethod = "standard", verbose = 0, directed=FALSE)
}
\arguments{
  \item{y}{Response Vector (length(y)=nrow(dm))}
  \item{dm}{Data matrix (rows: samples, cols: nodes)}
  \item{alpha}{Significance level of individual partial correlation tests}
  \item{corMethod}{"standard" or "Qn" for standard or robust correlation
  estimation}
  \item{verbose}{0-no output, 1-small output, 2-details (using 1 and 2
  makes the function very much slower)}
  \item{directed}{Boolean; should the output graph be directed?}
}
\value{
  \item{G}{A boolean vector indicating which column of \code{dm} is
  associated with \code{y}}
  \item{zMin}{The minimal z-values when testing partial correlations
  between \code{y} and each column of \code{dm}. The larger the number,
  the more consistent is the edge with the data.}
}
\details{
  This function basically applies \code{\link{pcAlgo}} on the data
  matrix obtained by joining \code{y} and \code{dm}. Since the output is
  not concerned with the edges found within the columns of \code{dm},
  the algorithm is adapted accordingly. Therefore, the runtime and the
  ability to deal with large datasets is typically increased quite a lot.
}
\references{
  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.
}
\seealso{\code{\link{pcAlgo}} which is the more general version of this function.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler.
}
\examples{
p <- 10
## generate and draw random DAG :
set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
plot(myDAG, main = "randomDAG(10, prob = 0.2)")

## generate 1000 samples of DAG using standard normal error distribution
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## let's pretend that the 10th column is the response and the first 9
## columns are explanatory variable. Which of the first 9 variables
## "cause" the tenth variable?
y <- d.mat[,10]
dm <- d.mat[,-10]
pcSelect(d.mat[,10],d.mat[,-10],alpha=0.05)
## You see, that variable 4,5,6 are considered as important
## By inspecting zMin you can also see, that the influence of variable 6
## is very evident from the data (zMin is 21.32, so quite large - as a
## rule
## of thumb for judging what is large, you could use quantiles of the
## Standard Normal Distribution)

## The result should be the same when using pcAlgo
resU <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard",directed=TRUE)
resU
plot(resU,zvalue.lwd=TRUE)
## as can be seen, the pcAlgo function also finds 4,5,6 as the important
## variables
## Again, variable 6 seems to be very evident from the data
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

\name{pcSelect}
\alias{pcSelect}
\title{PC-Select: Estimate subgraph around a response variable}

\description{ The goal is feature selection: If you
  have a response variable \eqn{y} and a data matrix \eqn{dm}, we want
  to know which variables are \dQuote{strongly influential} on \eqn{y}. The
  type of influence is the same as in the PC-Algorithm, i.e., \eqn{y}
  and \eqn{x} (a column of \eqn{dm}) are associated if they are
  correlated even when conditioning on any subset of the remaining
  columns in \eqn{dm}. Therefore, only very strong relations will be
  found and the result is typically a subset of other feature selection
  techniques. Note that there are also robust correlation methods
  available which render this method robust.
}
\usage{
pcSelect(y, dm, alpha, corMethod = "standard",
         verbose = FALSE, directed = FALSE)
}
\arguments{
  \item{y}{response vector.}
  \item{dm}{data matrix (rows: samples/observations, columns: variables);
    \code{nrow(dm) == length(y)}.}
  \item{alpha}{significance level of individual partial correlation tests.}
  \item{corMethod}{a string determining the method for correlation
    estimation via \code{\link{mcor}()}; specifically any of the
    \code{mcor(*, method = "..")} can be used, e.g., \code{"Qn"} for
    one kind of robust correlation estimate.}
  \item{verbose}{\code{\link{logical}} or in \eqn{\{0,1,2\}};
    \describe{
      \item{FALSE, 0:}{No output,}
      \item{TRUE, 1:}{Little output,}
      \item{2:}{Detailed output.}
    }
    Note that such diagnostic output may make the function considerably slower.
  }
  \item{directed}{logical; should the output graph be directed?}
}
\value{
  \item{G}{A \code{\link{logical}} vector indicating which column of
    \code{dm} is associated with \code{y}.}
  \item{zMin}{The minimal z-values when testing partial correlations
    between \code{y} and each column of \code{dm}.  The larger the
    number, the more consistent is the edge with the data.}
}
\details{
  This function basically applies \code{\link{pc}} on the data
  matrix obtained by joining \code{y} and \code{dm}.  Since the output is
  not concerned with the edges found within the columns of \code{dm},
  the algorithm is adapted accordingly.  Therefore, the runtime and the
  ability to deal with large datasets is typically increased
  substantially.
}

\references{
  Buehlmann, P., Kalisch, M. and Maathuis, M.H. (2010).
  Variable selection for high-dimensional linear models:
  partially faithful distributions and the PC-simple algorithm.
  \emph{Biometrika} \bold{97}, 261--278.
}
\seealso{
  \code{\link{pc}} which is the more general version of this function;
  \code{\link{pcSelect.presel}} which applies \code{pcSelect()} twice.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler.
}
\examples{
p <- 10
## generate and draw random DAG :
suppressWarnings(RNGversion("3.5.0"))
set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
if (require(Rgraphviz)) {
  plot(myDAG, main = "randomDAG(10, prob = 0.2)")
}
## generate 1000 samples of DAG using standard normal error distribution
n <- 1000
d.mat <- rmvDAG(n, myDAG, errDist = "normal")

## let's pretend that the 10th column is the response and the first 9
## columns are explanatory variable. Which of the first 9 variables
## "cause" the tenth variable?
y <- d.mat[,10]
dm <- d.mat[,-10]
(pcS <- pcSelect(d.mat[,10], d.mat[,-10], alpha=0.05))
## You see, that variable 4,5,6 are considered as important
## By inspecting zMin,
with(pcS, zMin[G])
## you can also see that the influence of variable 6
## is most evident from the data (its zMin is 18.64, so quite large - as
## a rule of thumb for judging what is large, you could use quantiles
## of the Standard Normal Distribution)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

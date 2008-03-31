\name{mcor}
\alias{mcor}
\title{Compute Correlation Matrix}
\description{
  Compute a correlation matrix, possibly by robust methods, applicable also
  for the case of a large number of variables.
}
\usage{
mcor(dm, method= c("standard", "Qn","QnStable", "ogkScaleTau2", "ogkQn","shrink"))
}
\arguments{
  \item{dm}{numeric matrix of data; rows are samples, columns are variables.}
  \item{method}{"standard" (default), "Qn", "QnStable", "ogkQn" and "shrink"
    envokes standard, 
    elementwise robust (based on \eqn{Q_n} scale estimator, see
    \code{\link[robustbase]{Qn}}), robust (Qn using OGK, see
    \code{\link[robustbase]{covOGK}}) or shrinked ()
    correlation estimate respectively.} 
}
\value{
  A correlation matrix estimated according to the specified method.
}
\details{
  The "standard" method envokes a standard correlation estimator. "Qn"
  envokes a robust, elementwise correlation estimator based on the Qn scale
  estimte. "QnStable" also uses the Qn scale estimator, but uses an
  improved way of transforming that into the correlation
  estimator. "ogkQn" envokes a correlation estimator based on Qn using
  OGK. "shrink" is only useful when used wiht pcSelect. An optimal
  shrinkage parameter is used. Only correlation between response and
  covariates is shrinked. 
}
\seealso{
  \code{\link[robustbase]{Qn}}) and \code{\link[robustbase]{covOGK}}
  from package \pkg{robustbase}.
  \code{\link{pcorOrder}} for computing partial correlations and
  \code{\link{condIndFisherZ}} for testing zero partial correlation.
}

\references{
  See those in the help pages for \code{Qn} and \code{covOGK} from package
  \pkg{robustbase}.
%   Ricardo A. Maronna and Ruben H. Zamar (2002). "Robust Estimates of
%   Location and Dispersion for High-Dimensional Datasets"; Technometrics,
%   Nov. 2002, Vol. 44, No. 4.
%
%   C. Croux and P.J. Rousseeuw (1992). "Time-efficient algorithms for two
%   highly robust estimators of scale"; Comput. Stat. 1, 411-428.
}

\author{
  Markus Kalisch \email{kalisch@stat.math.ethz.ch} and Martin Maechler
}
\examples{
## produce uncorrelated normal random variables
set.seed(42)
x <- rnorm(100)
y <- 2*x + rnorm(100)
## compute correlation of var1 and var2
mcor(cbind(x,y), method="standard")

## repeat but this time with heavy-tailed noise
yNoise <- 2*x + rcauchy(100)
mcor(cbind(x,yNoise), method="standard") ## shows almost no correlation
mcor(cbind(x,yNoise), method="Qn")       ## shows a lot correlation
mcor(cbind(x,yNoise), method="QnStable")       ## shows a lot correlation
mcor(cbind(x,yNoise), method="ogkQn")    ## dito
}
\keyword{multivariate}
\keyword{robust}

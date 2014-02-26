\name{rmvnorm.ivent}
\alias{rmvnorm.ivent}
\encoding{UTF-8}
\concept{Gaussian causal model}
\title{Simulate from a Gaussian Causal Model}
\description{
  Produces one or more samples from the observational or an interventional
  distribution associated to a Gaussian causal model.
}
\usage{
rmvnorm.ivent(n, object, target = integer(0), target.value = numeric(0))
}
\arguments{
  \item{n}{Number of samples required.}
  \item{object}{An instance of \code{\linkS4class{GaussParDAG}}}
  \item{target}{Intervention target: vector of intervened nodes.  If the vector
    is empty, samples from the observational distribution are generated.
    Otherwise, samples from an interventional distribution are simulated.}
  \item{target.value}{Values of the intervened variables.  If 
    \code{target.value} is a vector of the same length as \code{target}, the
    indicated intervention levels are used for all \code{n} samples.  If
    \code{target.value} is a matrix of dimension \code{n} by 
    \code{length(target)}, the \eqn{i}-th sample is simulated using the
    \eqn{i}-th row of the matrix as intervention levels.}
}
\value{
  If \code{n = 1} a vector of length \code{p} is returned, where \code{p} 
  denotes the number of nodes of \code{object}. Otherwise an \code{n} by
  \code{p} matrix is returned with one sample per row.
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\examples{
set.seed(307)
myDAG <- r.gauss.pardag(5, 0.5)
var(rmvnorm.ivent(n = 1000, myDAG))
myDAG$cov.mat()
var(rmvnorm.ivent(n = 1000, myDAG, target = 1, target.value = 1))
myDAG$cov.mat(target = 1, ivent.var = 0)
}
\keyword{datagen}
\keyword{models}

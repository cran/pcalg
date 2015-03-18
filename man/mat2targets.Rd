\name{mat2targets}
\alias{mat2targets}
\encoding{UTF-8}
\concept{intervention}
\concept{interventional data}
\title{Construct a list of intervention targets and a target index vector}
\description{
  This function constructs a list of intervention targets and a corresponding
  vector of target indices from a matrix specifying intervened vertices.  The
  input matrix has the same dimensions as the usual data matrix; the output
  can be used to create scoring objects (see \code{\linkS4class{Score}}) and
  to run causal inference methods based on interventional data such as 
  \code{\link{gies}} or \code{\link{simy}}.
}
\usage{
mat2targets(A)
}
\arguments{
  \item{A}{Logical matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n}
  is the sample size of a data set with jointly interventional and 
  observational data, and \eqn{p} is the number of variables.  \code{A[i, j]}
  is \code{TRUE} iff variable \code{j} is intervened in data point \code{i}.}
}
\value{
  \code{mat2targets} returns a list with two components:
  \item{targets}{A list of unique intervention targets.}
  \item{target.index}{A vector of intervention target indices. The intervention
    target of data point \code{i} is encoded as 
    \code{targets[[target.index[i]]]}.}
}
\author{
  Alain Hauser (\email{alain.hauser@bfh.ch})
}
\seealso{
  \code{\linkS4class{Score}}, \code{\link{gies}}, \code{\link{simy}}
}
\examples{
## Specify interventions using a matrix
p <- 5
n <- 10
A <- matrix(FALSE, nrow = n, ncol = p)
for (i in 1:n) A[i, (i-1) \%\% p + 1] <- TRUE

## Generate list of intervention targets and corresponding indices
target.list <- mat2targets(A)

for (i in 1:length(target.list$target.index))
  sprintf("Intervention target of \%d-th data point: \%d", 
    i, target.list$targets[[target.list$target.index[i]]]) 
}
\keyword{manip}
\keyword{list}

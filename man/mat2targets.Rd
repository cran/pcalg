\name{mat2targets}
\alias{mat2targets}
\alias{targets2mat}
\encoding{UTF-8}
\concept{intervention}
\concept{interventional data}
\title{Conversion between an intervention matrix and a list of intervention 
  targets}
\description{
  In a data set with \eqn{n} measurements of \eqn{p} variables, intervened
  variables can be specified in two ways:
  \itemize{
    \item with a \code{\link{logical}} intervention matrix of dimension
    \eqn{n \times p}{n × p}, where the entry \code{[i, j]} indicates whether 
    variable \eqn{j} has been intervened in measurement \eqn{i}; or
    
    \item with a list of (unique) intervention targets and a 
    \eqn{p}-dimensional vector indicating the indices of the intervention
    targets of the \eqn{p} measurements.
  }
  
  The function \code{mat2targets} converts the first representation to the
  second one, the function \code{targets2mat} does the reverse conversion.  The
  second representation can be used to create scoring objects (see 
  \code{\linkS4class{Score}}) and to run causal inference methods based on 
  interventional data such as \code{\link{gies}} or \code{\link{simy}}.
}
\usage{
mat2targets(A)
targets2mat(p, targets, target.index)
}
\arguments{
  \item{A}{Logical matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n}
  is the sample size of a data set with jointly interventional and 
  observational data, and \eqn{p} is the number of variables.  \code{A[i, j]}
  is \code{TRUE} iff variable \code{j} is intervened in data point \code{i}.}
  \item{p}{Number of variables}
  \item{targets}{List of unique intervention targets}
  \item{target.index}{Vector of intervention target indices.  The intervention
    target of data point \code{i} is encoded as 
    \code{targets[[target.index[i]]]}.}
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

## Convert back to matrix representation
all(A == targets2mat(p, target.list$targets, target.list$target.index))
}
\keyword{manip}
\keyword{list}

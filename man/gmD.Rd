\name{gmD}
\alias{gmD}
\docType{data}
\title{Graphical Model Discrete 5-Dim Example Data}
\description{
  This data set contains a matrix containing information on five
  discrete variables (levels are coded as numbers) and the corresonding
  DAG model.
}
\usage{data(gmD)}
\format{
  A \code{\link{list}} of two components
  \describe{
    \item{x:}{a \code{\link{data.frame}} with 5 columns \code{X1}
      .. \code{X5} each coding a discrete variable (aka
      \code{\link{factor}}) with interagesInt [1:10000, 1:5] 2 2 1 1 1 2 2 0 2 0 ...}
    \item{g:}{Formal class 'graphNEL' [package "graph"] with 6 slots\cr
      .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...\cr
      .. ..@ edgeL     :List of 5\cr
      ........
    }
  }
  where \code{x} is the data matrix and \code{g} is the DAG from which
  the data were generated.
}
\details{
The data was generated using Tetrad in the following way. A random DAG
on five nodes was generated; discrete variables were assigned to each
node (with 3, 2, 3, 4 and 2 levels); then conditional probability tables
corresponding
to the structure of the generated DAG were constructed. Finally, 10000
samples were drawn using the conditional probability tables.
}
% \references{
% %%  possibly secondary usages
% }
\examples{
data(gmD)
str(gmD, max=1)
if(require("Rgraphviz"))
  plot(gmD$ g, main = "gmD $ g --- the DAG of the gmD (10'000 x 5 discrete data)")
## >>>  1 --> 3 <-- 2 --> 4 --> 5
str(gmD$x)
## The number of unique values of each variable:
sapply(gmD$x, function(v) nlevels(as.factor(v)))
## X1 X2 X3 X4 X5
##  3  2  3  4  2
lapply(gmD$x, table) ## the (marginal) empirical distributions
## $X1
##    0    1    2
## 1933 3059 5008
##
## $X2
##    0    1
## 8008 1992
##
## $X3
## .....
}
\keyword{datasets}

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
  The format is a list of two components
  \describe{
    \item{x:}{Int [1:10000, 1:5] 2 2 1 1 1 2 2 0 2 0 ...}
    \item{g:}{Formal class 'graphNEL' [package "graph"] with 6 slots\cr
      .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...\cr
      .. ..@ edgeL     :List of 5\cr
      ........      
    }
  }
}
\details{
The data was generated using Tetrad in the following way. A random DAG
on five nodes was generated; discrete variables were assigned to each
node (with 3, 2, 3, 4 and 2 levels); then conditional probability tables
corresponding 
to the structure of the generated DAG were constructed. Finally, 10000
samples were drawn using the conditional probability tables.
}
% \source{
% %%  ~~ reference to a publication or URL from which the data were obtained ~~
% }
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(gmD)
## maybe str(gmD) ; plot(gmD) ...
}
\keyword{datasets}

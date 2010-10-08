\name{gmB}
\alias{gmB}
\docType{data}
\title{Graphical Model 5-Dim Binary Example Data}
\description{
This data set contains a matrix containing information on five binary
variables (coded as 0/1) and the corresonding DAG model.
}
\usage{data(gmB)}
\format{
  The format is a list of two components
  \describe{
    \item{x:}{int [1:5000, 1:5] 0 1 1 0 0 1 1 0 1 1 ...}
    \item{g:}{Formal class 'graphNEL' [package "graph"] with 6 slots\cr
      .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...\cr
      .. ..@ edgeL     :List of 5\cr
      ........      
    }
  }
}

\details{
The data was generated using Tetrad in the following way. A random DAG
on five nodes was generated; binary variables were assigned to each
node; then conditional probability tables corresponding
to the structure of the generated DAG were constructed. Finally, 5000
samples were drawn using the conditional probability tables.
}
% \source{
% %%  ~~ reference to a publication or URL from which the data were obtained ~~
% }
% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(gmB)
## maybe str(gmB) ; plot(gmB) ...
}
\keyword{datasets}

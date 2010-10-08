\name{gmL}
\alias{gmL}
\docType{data}
\title{Latent Variable Graphical Model Data Example}
\description{
This data set contains a matrix containing information on four gaussian
variables and the corresonding DAG model containing four observed and
one latent variable.
}
\usage{data(gmL)}
\format{
  The format is a list of 2 components
  \describe{
    \item{x:}{ $ x: num [1:10000, 1:4] 0.924 -0.189 1.016 0.363 0.497 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "2" "3" "4" "5"}
    \item{g:}{ $ g:Formal class 'graphNEL' [package "graph"] with 6 slots
  .. ..@ nodes     : chr [1:5] "1" "2" "3" "4" ...
  .. ..@ edgeL     :List of 5
  ........
      }
  }
}
\details{
The data was generated as indicated below. First, a random DAG model was
generated with five nodes; then 10000 samples were drawn from this
model; finally, variable one was declared to be latent and
the corresponding column was deleted from the simulated data set.
}
\source{
  % \code{
  \preformatted{
    ## Used to generate "gmL"
    set.seed(47)
    p <- 5
    n <- 10000
    gGtrue <- randomDAG(p, prob = 0.3) ## true DAG
    myX <- rmvDAG(n, gGtrue)
    colnames(myX) <- as.character(1:5)
    gmL <- list(x = myX[,-1], g = gGtrue)
  }
  %}% code
}
% \references{
 %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(gmL)
str(gmL, max=3)

## the graph:
gmL$g
graph::nodes(gmL$g) ; str(graph::edges(gmL$g))
if(require("Rgraphviz"))
  plot(gmL$g, main = "gmL $ g -- latent variable example data")

pairs(gmL $x) # the data
}
\keyword{datasets}

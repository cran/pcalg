\name{dsep}
\alias{dsep}
\title{Test for d-separation in a DAG}
\description{
  This function tests for d-separation of nodes in a DAG.
}
\usage{
dsep(a, b, S, g, john.pairs = NA)
}
\arguments{
  \item{a}{Label (sic!) of node A}
  \item{b}{Label (sic!) of node B}
  \item{S}{Labels (sic!) of set of nodes on which is conditioned}
  \item{g}{Element "g" containing
    the Directed Acyclic Graph (object of \code{\link{class}}
    "graph", see \code{\link[graph]{graph-class}} from the
    package \pkg{graph})}
  \item{john.pairs}{the shortest path distance matrix for all pairs of
    nodes as computed by \code{\link[RBGL]{johnson.all.pairs.sp}} from
    package \pkg{RBGL}.}
}
\value{TRUE if a and b are d-separated by S in G, otherwise FALSE.
}

\details{This function checks separation in the moralized graph as
  explained in Lauritzen (2004).}

\references{
  S.L. Lauritzen (2004),
  Graphical Models, \emph{Oxford University Press}, Chapter 3.2.2
}

\seealso{\code{\link{dsepTest}} for a wrapper of this function that can
  easily be included into \code{\link{skeleton}}, \code{\link{pc}} or
      \code{\link{fci}}
    }

\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
## generate random DAG
p <- 8
set.seed(45)
myDAG <- randomDAG(p, prob = 0.3)
if (require(Rgraphviz)) {
plot(myDAG)
}

## Examples for d-separation
dsep("1","7",NULL,myDAG)
dsep("4","5",NULL,myDAG)
dsep("4","5","2",myDAG)
dsep("4","5",c("2","3"),myDAG)

## Examples for d-connection
dsep("1","3",NULL,myDAG)
dsep("1","6","3",myDAG)
dsep("4","5","8",myDAG)
}

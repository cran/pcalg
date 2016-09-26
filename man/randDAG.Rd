\name{randDAG}
\title{Random DAG Generation}
\alias{randDAG}
\concept{DAG}
\concept{random directed acyclic graph}
\description{
  Generating random directed acyclic graphs (DAGs) with fixed expected
  number of neighbours.  Several different methods are provided, each
  intentionally biased towards certain properties.  The methods are based
  on the analogue \code{*.game} functions in the \pkg{igraph} package.
}
\usage{
randDAG(n, d, method ="er", par1=NULL, par2=NULL,
        DAG = TRUE, weighted = TRUE, wFUN = list(runif, min=0.1, max=1))
}
\arguments{
  \item{n}{integer, at least \code{2}, indicating the number of nodes in
    the DAG.} 
  \item{d}{a positive number, corresponding to the expected number of
    neighbours per node, more precisely the expected sum of the in- and
    out-degree.} 
  \item{method}{a string, specifying the method used for generating the
    random graph.  See details below.}
  \item{par1, par2}{optional additional arguments, dependent on the
    method.  See details.}
  \item{DAG}{logical, if \code{TRUE}, labelled graph is directed to a labelled acyclic graph.}
  \item{weighted}{logical indicating if edge weights are computed according to \code{wFUN}.}
  \item{wFUN}{a \code{\link{function}} for computing the edge weights in
    the DAG.  It takes as first argument a number of edges \code{m} for
    which it returns a vector of length \code{m} containing the
    weights.  Alternatively, \code{wFUN} can be a \code{\link{list}}
    consisting of the function in the first entry and of further
    arguments of the function in the additional entries.  The default
    (only if \code{weighted} is true) is a uniform weight in \eqn{[0.1,
      1]}.  See the examples for more.}
}
\details{
  A (weighted) random graph with \code{n} nodes and expected number of
  neighbours \code{d} is constructed.  For \code{DAG=TRUE}, the graph is
  oriented to a DAG.  There are eight different random graph models
  provided, each selectable by the parameters \code{method},
  \code{par1} and \code{par2}, with \code{method}, a string,
  taking one of the following values:
  \describe{
    \item{\code{regular}:}{Graph where every node has exactly \code{d}
      incident edges. \code{par1} and \code{par2} are not used.}

    \item{\code{watts}:}{Watts-Strogatz graph that interpolates between
      the regular (\code{par1->0}) and Erdoes-Renyi graph
      (\code{par1->1}).  The parameter \code{par1} is per default
      \code{0.5} and has to be in \code{(0,1)}. \code{par2} is not used.}

    \item{\code{er}:}{Erdoes-Renyi graph where every edge is present
    independently. \code{par1} and \code{par2} are not used.}

    \item{\code{power}:}{A graph with power-law degree distribution with
    expectation \code{d}.\code{par1} and \code{par2} are not used.} 

    \item{\code{bipartite}:}{Bipartite graph with at least \code{par1*n}
      nodes in group 1 and at most \code{(1-par1)*n} nodes in group 2.
      The argument \code{par1} has to be in \code{[0,1]} and is per
      default \code{0.5}. \code{par2} is not used.}

    \item{\code{barabasi}:}{A graph with power-law degree distribution
      and preferential attachement according to parameter \code{par1}.  It
      must hold that \code{par1 >= 1} and the default is
      \code{par1=1}. \code{par2} is not used.}

    \item{\code{geometric}:}{A geometric random graph in dimension
      \code{par1}, where  \code{par1} can take values from
      \code{{2,3,4,5}} and is per default \code{2}.  If \code{par2="geo"}
      and \code{weighted=TRUE}, then the weights are computed according to
      the Euclidean distance. There are currently no other option for
      \code{par2} implemented.}

    \item{\code{interEr}:}{A graph with \code{par1} islands of
      Erdoes-Renyi graphs, every pair of those connected by a certain
      number of edges proportional to \code{par2} (fraction of
      inter-connectivity).  It is required that
      \eqn{n/s} be integer and \code{par2} in \eqn{(0,1)}.  Defaults are
      \code{par1=2} and \code{par2=0.25}, respectively.}
  }
}
\value{
  A graph object of class \code{\link[graph]{graphNEL}}.
}
\references{
  These methods are mainly based on the analogue functions in the \pkg{igraph} package.
}

\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Manuel Schuerch.}

\note{
  The output is \emph{not} topologically sorted (as opposed to the
  output of \code{\link{randomDAG}}).
}%% MM: --> ../TODO - we should add optional sorting

\seealso{the package \code{\link{igraph}}, notably help pages such as
  \code{\link[igraph]{random.graph.game}} or \code{\link[igraph]{barabasi.game}};

  \code{\link{unifDAG}} for generating uniform random DAGs.

  \code{\link{randomDAG}} a limited and soon deprecated version of \code{\link{randDAG}};
  \code{\link{rmvDAG}} for generating multivariate data according to a DAG.
}
\examples{
set.seed(37)
dag1 <- randDAG(10, 4, "regular")
dag2 <- randDAG(10, 4, "watts")
dag3 <- randDAG(10, 4, "er")
dag4 <- randDAG(10, 4, "power")
dag5 <- randDAG(10, 4, "bipartite")
dag6 <- randDAG(10, 4, "barabasi")
dag7 <- randDAG(10, 4, "geometric")
dag8 <- randDAG(10, 4, "interEr", par2 = 0.5)

## require("Rgraphviz")
par(mfrow=c(4,2))
plot(dag1,main="Regular graph")
plot(dag2,main="Watts-Strogatz graph")
plot(dag3,main="Erdoes-Renyi graph")
plot(dag4,main="Power-law graph")
plot(dag5,main="Bipartite graph")
plot(dag6,main="Barabasi graph")
plot(dag7,main="Geometric random graph")
plot(dag8,main="Interconnected island graph")

set.seed(45)
dag0 <- randDAG(6,3)
dag1 <- randDAG(6,3, weighted=FALSE)
dag2 <- randDAG(6,3, DAG=FALSE)
par(mfrow=c(1,2))
plot(dag1)
plot(dag2)      ## undirected graph
dag0@edgeData   ## note the uniform weights between 0.1 and 1
dag1@edgeData   ## note the constant weights

wFUN <- function(m,lB,uB) { runif(m,lB,uB) }
dag <- randDAG(6,3,wFUN=list(wFUN,1,4))
dag@edgeData   ## note the uniform weights between 1 and 4
}
\keyword{graphs}
\keyword{models}

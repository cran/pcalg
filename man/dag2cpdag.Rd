\name{dag2cpdag}
\alias{dag2cpdag}
\title{Convert a DAG to a CPDAG}
\description{
  Convert a DAG (Directed Acyclic Graph) to a Completed Partially
  Directed Acyclic Graph (CPDAG).
}
\usage{
dag2cpdag(g)
}
\arguments{
  \item{g}{an \R object of class \code{"graph"} (package \pkg{graph}),
    representing a DAG.}
}
\details{
  This function converts a DAG into its corresponding
  (unique) CPDAG as follows.  Because every DAG in the
  Markov equivalence class described by a CPDAG shares the same skeleton
  and the same v-structures, this function takes the skeleton and the
  v-structures of the given DAG \code{g}.  Afterwards it simply uses the
  3 orientation rules of the PC algorithm (see references) to orient as
  many of the remaining undirected edges as possible.
  
  The function is a simple wrapper function for \code{\link{dag2essgraph}}
  which is more powerfull since it also allows the calculation of the 
  Markov equivalence class in the presence of interventional data.

  The output of this function is exactly the same as the one using
  \preformatted{pc(suffStat, indepTest, alpha, labels)}
  using the true correlation matrix in the function \code{gaussCItest}
  with a large virtual sample size and a large alpha, but it is much
  faster.
}
\value{
  A graph object containing the CPDAG.
}
\references{
  C. Meek (1995). Causal inference and causal explanation with
  background knowledge. In \emph{Proceedings of the Eleventh Conference on
  Uncertainty in Artificial Intelligence (UAI-95)}, pp. 403-411. Morgan
  Kaufmann Publishers, Inc.

  P. Spirtes, C. Glymour and R. Scheines (2000)
  \emph{Causation, Prediction, and Search}, 2nd edition, The MIT Press.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and 
  Alain Hauser(\email{alain.hauser@bfh.ch})}
\seealso{\code{\link{dag2essgraph}}, \code{\link{randomDAG}}, \code{\link{pc}}}
\examples{
## A -> B <- C
am1 <- matrix(c(0,1,0, 0,0,0, 0,1,0), 3,3)
colnames(am1) <- rownames(am1) <- LETTERS[1:3]
g1 <- as(t(am1), "graphNEL") ## convert to graph
cpdag1 <- dag2cpdag(g1)

if(requireNamespace("Rgraphviz")) {
    par(mfrow = c(1,2))
    plot(g1)
    plot(cpdag1)
}

## A -> B -> C
am2 <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3,3)
colnames(am2) <- rownames(am2) <- LETTERS[1:3]
g2 <- as(t(am2), "graphNEL") ## convert to graph
cpdag2 <- dag2cpdag(g2)

if(requireNamespace("Rgraphviz")) {
    par(mfrow = c(1,2))
    plot(g2)
    plot(cpdag2)
}
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

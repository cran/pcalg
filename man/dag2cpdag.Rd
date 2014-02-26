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
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Diego Colombo}
\seealso{\code{\link{randomDAG}}, \code{\link{pc}}}
\examples{
p <- 10  ## number of random variables
s <- 0.4 ## sparseness of the graph

## generate a random DAG
set.seed(42)
g <- randomDAG(p, s)
g01 <- 1*(as(g,"matrix") > 0) # 0/1-version of adjacency matrix
print.table(g01, zero.=".")

## compute its unique CPDAG
system.time(
  res <- dag2cpdag(g)
)
## res has some bidirected edges
##  ==> adj.matrix: no longer upper triangular, but {0,1}
print.table(as(res, "matrix"), zero.=".")
dm <- as(res, "matrix") - g01 ## difference: 2 entries in lower tri.
print.table(dm, zero.=".")
stopifnot(all(dm \%in\% 0:1), sum(dm == 1) == 2)

## Find CPDAG with PC algorithm:
## As dependence oracle, we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
system.time(
rpc <- pc(suffStat = list(C = cov2cor(trueCov(g)), n = 10^9),
          indepTest = gaussCItest, alpha = 0.9999, p = p)
)
## confirm that it coincides with dag2cpdag()'s result:
stopifnot(all.equal(as( res,     "matrix"),
                    as(rpc@graph,"matrix"), tol=0))
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

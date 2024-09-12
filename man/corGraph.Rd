\name{corGraph}
\alias{corGraph}
\title{Computing the correlation graph}
\description{
  Computes the correlation graph.  This is the graph in which an edge is
  drawn between node i and node j, if the null hypothesis \dQuote{\emph{Correlation
      between \eqn{X_i} and \eqn{X_j} is zero}} can be rejected at the
  given significance level \eqn{\alpha (alpha)}{alpha}.
}
\usage{
corGraph(dm, alpha=0.05, Cmethod="pearson")
}
\arguments{
  \item{dm}{numeric matrix with rows as samples and columns as variables.}
  \item{alpha}{significance level for correlation test (numeric)}
  \item{Cmethod}{a \code{\link{character}} string indicating which
    correlation coefficient is to be used for the test.  One of
    \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}, can be
    abbreviated.}
}
\value{
  Undirected correlation graph, a \code{\link[graph]{graph-class}} object
  (package \pkg{graph});  \code{\link{getGraph}} for the \dQuote{fitted}
  graph.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler
}
\examples{
## create correlated samples
x1 <- rnorm(100)
x2 <- rnorm(100)
mat <- cbind(x1,x2, x3 = x1+x2)

if (require(Rgraphviz)) {
## ``analyze the data''
(g <- corGraph(mat)) # a 'graphNEL' graph, undirected
plot(g) # ==> (1) and (2) are each linked to (3)

## use different significance level and different method
(g2 <- corGraph(mat, alpha=0.01, Cmethod="kendall"))
plot(g2) ## same edges as 'g'
}
}
\keyword{graphs}
\keyword{models}

\name{ges}
\alias{ges}
\encoding{UTF-8}
\concept{greedy equivalence search}
\concept{essential graph}
\concept{CPDAG}
\title{Estimate the Markov equivalence class of a DAG using GES}
\description{
  Estimate the observational essential graph representing the Markov
  equivalence class of a DAG using the greedy equivalence search (GES)
  algorithm of Chickering (2002).
}
\usage{
ges(p, score, fixedGaps = NULL,
    turning = TRUE, maxDegree = integer(0), verbose = FALSE, ...)
}
\arguments{
  \item{p}{Number of variables.}
  \item{score}{An instance of a class derived from \code{\linkS4class{Score}}
    which only accounts for observational data.}
  \item{fixedGaps}{logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i, j]} is \code{TRUE}, the result is guaranteed to have no edge
    between nodes \eqn{i} and \eqn{j}.}
  \item{turning}{Logical indicating whether the function should try to augment
    the score by turning edges (cf. details).}
  \item{maxDegree}{Parameter used to limit the vertex degree of the estimated
    graph.  Valid arguments:
    \enumerate{
      \item Vector of length 0 (default): vertex degree is not limited.
      \item Real number \eqn{r}, \eqn{0 < r < 1}: degree of vertex \eqn{v} is
        limited to \eqn{r \cdot n_v}, where \eqn{n_v} denotes the number of
        data points where \eqn{v} was not intervened.
      \item Single integer: uniform bound of vertex degree for all vertices
        of the graph.
      \item Integer vector of length \code{p}: vector of individual bounds
        for the vertex degrees.
    }
    }
  \item{verbose}{if \code{TRUE}, detailed output is provided.}
  \item{\dots}{Additional arguments for debugging purposes and fine tuning.}
}
\details{
  Under the assumption that the distribution of the observed variables is
  faithful to a DAG, this function estimates the Markov equivalence class of
  the DAG. It does not estimate the DAG itself, because this is typically
  impossible (even with an infinite amount of data): different DAGs (forming a
  Markov equivalence class) can describe the same conditional independence
  relationships and be statistically indistinguishable from observational data
  alone.

  All DAGs in an equivalence class have the same skeleton (i.e., the same
  adjacency information) and the same v-structures (i.e., the same induced
  subgraphs of the form \eqn{a \longrightarrow b \longleftarrow c}{a → b ← c}).
  However, the direction of some edges may be undetermined, in the sense that
  they point one way in one DAG in the equivalence class, while they point the
  other way in another DAG in the equivalence class.

  An equivalence class can be uniquely represented by a partially directed graph
  called (observational) essential graph or CPDAG (completed partially directed
  acyclic graph).  Its edges have the following interpretation:
  \enumerate{
    \item a directed edge \eqn{a \longrightarrow b}{a → b} stands for an arrow
    that has the same orientation in all representatives of the Markov
    equivalence class;
    \item an undirected edge a -- b stands for an arrow that is oriented in one
    way in some representatives of the equivalence class and in the other way
    in other representatives of the equivalence class.
  }
  Note that when plotting the object, undirected and bidirected edges are
  equivalent.

  GES (greedy equivalence search) is a score-based algorithm that greedily
  maximizes a score function (typically the BIC, passed to the function via the
  argument \code{score}) in the space of (observational) essential graphs in
  three phases, starting from the empty graph:
  \describe{
    \item{Forward phase}{In the forward phase, GES moves through the space of
    essential graphs in steps that correspond to the addition of a single edge
    in the space of DAGs; the phase is aborted as soon as the score cannot be
    augmented any more.}
    \item{Backward phase}{In the backward phase, the algorithm performs moves
    that correspond to the removal of a single edge in the space of DAGs until
    the score cannot be augmented any more.}
    \item{Turning phase}{In the turning phase, the algorithm performs moves
    that correspond to the reversal of a single arrow in the space of DAGs
    until the score cannot be augmented any more.}
  }
  GES cycles through these three phases until no augmentation of the score is
  possible any more.  Note that the turning phase (activated with
  \code{turning = TRUE}, the default behaviour) was not part of the original
  implementation of Chickering (2002), but was introduced by Hauser and
  Bühlmann (2012) and shown to improve the overall estimation performance.

  GES has the same purpose as the PC algorithm (see \code{\link{pc}}).  While
  the PC algorithm is based on conditional independence tests (requiring the
  choice of an independence test and a significance level, see
  \code{\link{pc}}), the GES algorithm is a score-based method (requiring the
  choice of a score function) and does not depend on conditional independence
  tests.  Since GES always operates in the space of essential graphs, it
  returns a valid essential graph (or CPDAG) in any case.
}
\value{
  \code{ges} returns a list with the following two components:
  \item{essgraph}{An object of class \code{\linkS4class{EssGraph}} containing an
    estimate of the equivalence class of the underlying DAG.}
  \item{repr}{An object of a class derived from \code{\linkS4class{ParDAG}}
    containing a (random) representative of the estimated equivalence class.}
}
\references{
  D.M. Chickering (2002).  Optimal structure identification with greedy search.
  \emph{Journal of Machine Learning Research} \bold{3}, 507--554

  A. Hauser and P. Bühlmann (2012).  Characterization and greedy learning of
  interventional Markov equivalence classes of directed acyclic graphs.
  \emph{Journal of Machine Learning Research} \bold{13}, 2409--2464.

  P. Spirtes, C.N. Glymour, and R. Scheines (2000).  \emph{Causation,
  Prediction, and Search}, MIT Press, Cambridge (MA).
}
\author{
  Alain Hauser (\email{alain.hauser@bfh.ch})
}
\seealso{
  \code{\link{pc}}, \code{\linkS4class{Score}}, \code{\linkS4class{EssGraph}}
}
\examples{
## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x)

## Estimate the essential graph
ges.fit <- ges(ncol(gmG8$x), score)

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(ges.fit$essgraph, main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(ges.fit, max=2)
  as(as(ges.fit$essgraph,"graphNEL"),"Matrix")
}
}
\keyword{models}
\keyword{graphs}

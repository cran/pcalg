\name{opt.target}
\alias{opt.target}
\encoding{UTF-8}
\concept{essential graph}
\concept{CPDAG}
\concept{interventions}
\title{Get an optimal intervention target}
\description{
  Given a (observational or interventional) essential graph (or "CPDAG"), find
  the optimal intervention target that maximizes the number of edges that can
  be oriented after the intervention.
}
\usage{opt.target(essgraph, max.size, use.node.names = TRUE)}
\arguments{
  \item{essgraph}{An \code{\linkS4class{EssGraph}} or
    \code{\linkS4class{graphNEL}} object representing a (observational or
    interventional) essential graph (or "CPDAG").}
  \item{max.size}{Maximum size of the intervention target. Only 1 and the
    number of vertices of \code{essgraph} are allowed; the latter means no size
    limit is applied (the default if the parameter is missing).}
  \item{use.node.names}{Indicates if the intervention target should be
    returned as a list of node names (if \code{TRUE}) or indices (if
    \code{FALSE}).}
}
\details{
  This function implements active learning strategies for structure learning
  from interventional data, one that calculates an optimal single-vertex
  intervention target, and one that calculates an optimal intervention target
  of arbitrary size. "Optimal" means the proposed intervention target guarantees
  the highest number of edges that can be oriented after performing the
  intervention, assuming the essential graph provided as input is the true
  essential graph under the currently available interventional data (i.e.,
  neglecting possible estimation errors).

  Implementation corresponds to algorithms "OptSingle" and "OptUnb" published
  in Hauser and Bühlmann (2012).
}
\value{
  A character vector of node names (if \code{use.node.names = TRUE}), or an
  integer vector of node indices (if \code{use.node.names = FALSE}) indicating
  the optimal intervention target.
}
\references{
  A. Hauser and P. Bühlmann (2012). Two optimal strategies for active learning
  of causal models from interventions. \emph{Proceedings of the 6th European
  Workshop on Probabilistic Graphical Models (PGM-2012)}, 123--130
}
\author{
  Alain Hauser (\email{alain.hauser@math.ethz.ch})
}
\seealso{
  \code{\linkS4class{EssGraph}}
}
\examples{
## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x)

## Estimate the essential graph using GES
ges.fit <- ges(score)
essgraph <- ges.fit$essgraph

## Plot the estimated essential graph
if (require(Rgraphviz)) {
  plot(essgraph, main = "Estimated CPDAG")
}
## The CPDAG has 1 unoriented component with 3 edges (Author <-> Bar, Bar <->
## Ctrl, Bar <-> V5)

## Get optimal single-vertex and unbounded intervention target
opt.target(essgraph, max.size = 1)
opt.target(essgraph, max.size = essgraph$node.count())
}

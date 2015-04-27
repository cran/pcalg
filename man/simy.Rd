\name{simy}
\alias{simy}
\concept{dynamic programming}
\concept{interventional Markov equivalence}
\encoding{UTF-8}
\title{Estimate Interventional Markov Equivalence Class of a DAG}
\description{
  Estimate the interventional essential graph representing the Markov
  equivalence class of a DAG using the dynamic programming (DP) approach of
  Silander and Myllymäki (2006).  This algorithm maximizes a decomposable
  scoring criterion in exponential runtime.
}
\usage{
simy(p, targets, score, verbose = FALSE, ...)
}
\arguments{
  \item{p}{Number of variables.}
  \item{targets}{A list of intervention targets (cf. details).  A list of
    vectors, each vector listing the vertices of one intervention target.}
  \item{score}{An instance of a class derived from \code{\linkS4class{Score}}.}
  \item{verbose}{if \code{TRUE}, detailed output is provided.}
  \item{\dots}{Additional arguments for debugging purposes and fine tuning.}
}
\details{
  This function estimates the interventional Markov equivalence class of a DAG
  based on a data sample with interventional data originating from various
  interventions and possibly observational data. The intervention targets used
  for data generation must be specified by the argument \code{targets} as a
  list of (integer) vectors listing the intervened vertices; observational
  data is specified by an empty set, i.e. a vector of the form
  \code{integer(0)}.  As an example, if data contains observational samples
  as well as samples originating from an intervention at vertices 1 and 4,
  the intervention targets must be specified as \code{list(integer(0),
  as.integer(1), as.integer(c(1, 4)))}.

  An interventional Markov equivalence class of DAGs can be uniquely
  represented by a partially directed graph called interventional essential
  graph.  Its edges have the following interpretation:
  \enumerate{
    \item a directed edge \eqn{a \longrightarrow b}{a → b} stands for an arrow that has the same
      orientation in all representatives of the interventional Markov
      equivalence class;
    \item an undirected edge a -- b stands for an arrow that is oriented in one
      way in some representatives of the equivalence class and in the other way
      in other representatives of the equivalence class.
  }
  Note that when plotting the object, undirected and bidirected edges are
  equivalent.

  The DP approach of Silander and Myllymäki (2006) is a score-based algorithm
  that guarantees to find the optimum of any decomposable scoring criterion.
  Its CPU and memory consumption grow exponentially with the number of
  variables in the system, irrespective of the sparseness of the true or
  estimated DAG.  The implementation in the pcalg package is feasible up to
  approximately 20 variables, depending on the user's computer.
}
\value{
  \code{simy} returns a list with the following two components:
  \item{essgraph}{An object of class \code{\linkS4class{EssGraph}} containing an
    estimate of the equivalence class of the underlying DAG.}
  \item{repr}{An object of a class derived from \code{\linkS4class{ParDAG}}
    containing a (random) representative of the estimated equivalence class.}
}
\references{
  T. Silander and P. Myllymäki (2006).  A simple approach for finding the
  globally optimal Bayesian network structure.  \emph{Proceedings of the 22nd
  Conference on Uncertainty in Artificial Intelligence (UAI 2006)}, 445--452
}
\author{
  Alain Hauser (\email{alain.hauser@bfh.ch})
}
\seealso{
  \code{\link{gies}}, \code{\linkS4class{Score}}, \code{\linkS4class{EssGraph}}
}
\examples{
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
data(gmInt)

## Define the score (BIC)
score <- new("GaussL0penIntScore", gmInt$x, gmInt$targets, gmInt$target.index)

## Estimate the essential graph
simy.fit <- simy(ncol(gmInt$x), gmInt$targets, score)
eDAG <- simy.fit$essgraph
as(eDAG, "graph")

## Look at the graph incidence matrix (a "sparseMatrix"):
if(require(Matrix))
  show( as(as(eDAG, "graphNEL"), "Matrix") )


## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(eDAG, main = "Estimated ess. graph")
  plot(gmInt$g, main = "True DAG")
}
}
\keyword{graphs}

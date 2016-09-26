\name{gies}
\alias{gies}
\encoding{UTF-8}
\concept{greedy interventional equivalence search}
\concept{essential graph}
\title{Estimate Interventional Markov Equivalence Class of a DAG by GIES}
\description{
  Estimate the interventional essential graph representing the Markov
  equivalence class of a DAG using the greedy interventional equivalence search
  (GIES) algorithm of Hauser and Bühlmann (2012).
}
\usage{
gies(score, labels = score$getNodes(), targets = score$getTargets(),
     fixedGaps = NULL, adaptive = c("none", "vstructures", "triples"), 
     phase = c("forward", "backward", "turning"), iterate = length(phase) > 1,
     turning = TRUE, maxDegree = integer(0), verbose = FALSE, ...)
}
\arguments{
  \item{score}{An \R object inheriting from \code{\linkS4class{Score}}.}
  \item{labels}{Node labels; by default, they are determined from the scoring
    object.}
  \item{targets}{A list of intervention targets (cf. details).  A list of
    vectors, each vector listing the vertices of one intervention target.}
  \item{fixedGaps}{Logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i, j]} is \code{TRUE}, the result is guaranteed to have no edge
    between nodes \eqn{i} and \eqn{j}.}
  \item{adaptive}{indicating whether constraints should be adapted to
    newly detected v-structures or unshielded triples (cf. details).}
  \item{phase}{Character vector listing the phases that should be used; possible
    values: \code{forward}, \code{backward}, and \code{turning} (cf. details).}
  \item{iterate}{Logical indicating whether the phases listed in the argument
    \code{phase} should be iterated more than once (\code{iterate = TRUE}) or
    not.}
  \item{turning}{Setting \code{turning = TRUE} is equivalent to setting
    \code{phases = c("forward", "backward")} and \code{iterate = FALSE}; the
    use of the argument \code{turning} is deprecated.}
  \item{maxDegree}{Parameter used to limit the vertex degree of the estimated
    graph.  Possible values:
    \enumerate{
      \item Vector of length 0 (default): vertex degree is not limited.
      \item Real number \eqn{r}, \eqn{0 < r < 1}: degree of vertex \eqn{v} is
        limited to \eqn{r \cdot n_v}{r . n_v}, where \eqn{n_v} denotes
	the number of data points where \eqn{v} was not intervened.
      \item Single integer: uniform bound of vertex degree for all vertices
        of the graph.
      \item Integer vector of length \code{p}: vector of individual bounds
        for the vertex degrees.
    }
  }
  \item{verbose}{If \code{TRUE}, detailed output is provided.}
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
    \item a directed edge \eqn{a \longrightarrow b}{a → b} stands for an arrow
      that has the same orientation in all representatives of the
      interventional Markov equivalence class;
    \item an undirected edge \eqn{a} -- \eqn{b} stands for an arrow that is
      oriented in one  way in some representatives of the equivalence class and
      in the other way in other representatives of the equivalence class.
  }
  Note that when plotting the object, undirected and bidirected edges are
  equivalent.

  GIES (greedy interventional equivalence search) is a score-based algorithm
  that greedily maximizes a score function (typically the BIC, passed to the
  function via the argument \code{score}) in the space of interventional
  essential graphs in three phases, starting from the empty graph:
  \describe{
    \item{Forward phase}{In the forward phase, GIES moves through the space of
      interventional essential graphs in steps that correspond to the addition
      of a single edge in the space of DAGs; the phase is aborted as soon as
      the score cannot be augmented any more.}
    \item{Backward phase}{In the backward phase, the algorithm performs moves
      that correspond to the removal of a single edge in the space of DAGs
      until the score cannot be augmented any more.}
    \item{Turning phase}{In the turning phase, the algorithm performs moves
      that correspond to the reversal of a single arrow in the space of DAGs
      until the score cannot be augmented any more.}
  }
  The phases that are actually run are specified with the argument 
  \code{phase}.  GIES cycles through the specified phases until no augmentation 
  of the score is possible any more if \code{iterate = TRUE}.  GIES is an 
  interventional extension of the GES (greedy equivalence search) algorithm of 
  Chickering (2002) which is limited to observational data and hence operates 
  on the space of observational instead of interventional Markov equivalence 
  classes.
  
  Using the argument \code{fixedGaps}, one can make sure that certain edges
  will \emph{not} be present in the resulting essential graph: if the entry
  \code{[i, j]} of the matrix passed to \code{fixedGaps} is \code{TRUE}, there
  will be no edge between nodes \eqn{i} and \eqn{j}.  Using this argument 
  can speed up the execution of GIES and allows the user to account for
  previous knowledge or other constraints.  The argument \code{adaptive} can be
  used to relax the constraints encoded by \code{fixedGaps} as follows:
  \itemize{
    \item When \code{adaptive = "vstructures"} and the algorithm introduces a 
    new v-structure \eqn{a \longrightarrow b \longleftarrow c}{a → b ← c} in the 
    forward phase, then the edge \eqn{a - c} is removed from the list of fixed 
    gaps, meaning that the insertion of an edge between \eqn{a} and \eqn{c} 
    becomes possible even if it was forbidden by the initial matrix passed to 
    \code{fixedGaps}.
    
    \item When \code{adaptive = "triples"} and the algorithm introduces a new
    unshielded triple in the forward phase (i.e., a subgraph of three nodes
    \eqn{a}, \eqn{b} and \eqn{c} where \eqn{a} and \eqn{b} as well as \eqn{b}
    and \eqn{c} are adjacent, but \eqn{a} and \eqn{c} are not), then the edge
    \eqn{a - c} is removed from the list of fixed gaps.
  }
  This modifications of the forward phase of GIES are inspired by the 
  analog modifications in the forward phase of GES, which makes the successive
  application of a skeleton estimation method and GES restricted to an 
  estimated skeleton a \emph{consistent} estimator of the DAG (cf. Nandy,
  Hauser and Maathuis, 2015). 
}
\value{
  \code{gies} returns a list with the following two components:
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
  
  P. Nandy, A. Hauser and M. Maathuis (2015).  Understanding consistency in 
  hybrid causal structure learning.  \emph{arXiv preprint} 1507.02608
}
\author{
  Alain Hauser (\email{alain.hauser@bfh.ch})
}
\seealso{
  \code{\link{ges}}, \code{\linkS4class{Score}}, \code{\linkS4class{EssGraph}}
}
\examples{
## Load predefined data
data(gmInt)

## Define the score (BIC)
score <- new("GaussL0penIntScore", gmInt$x, gmInt$targets, gmInt$target.index)

## Estimate the essential graph
gies.fit <- gies(score)

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(gies.fit$essgraph, main = "Estimated ess. graph")
  plot(gmInt$g, main = "True DAG")
}
}
\keyword{models}
\keyword{graph}

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
gies(p, targets, score, fixedGaps = NULL, 
  turning = TRUE, maxDegree = integer(0), verbose = FALSE, ...)
}
\arguments{
  \item{p}{Number of variables.}
  \item{targets}{A list of intervention targets (cf. details).  A list of
    vectors, each vector listing the vertices of one intervention target.}
  \item{score}{An instance of a class derived from \code{\linkS4class{Score}}.}
  \item{fixedGaps}{logical \emph{symmetric} matrix of dimension p*p.  If entry
    \code{[i, j]} is \code{TRUE}, the result is guaranteed to have no edge
    between nodes \eqn{i} and \eqn{j}.}
  \item{turning}{Logical indicating whether the function should try to augment
    the score by turning edges (cf. details).}
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
  GIES cycles through these three phases until no augmentation of the score is
  possible any more.  GIES is an interventional extension of the GES (greedy
  equivalence search) algorithm of Chickering (2002) which is limited to
  observational data and hence operates on the space of observational instead
  of interventional Markov equivalence classes.
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
gies.fit <- gies(ncol(gmInt$x), gmInt$targets, score) 

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(gies.fit$essgraph, main = "Estimated ess. graph")
  plot(gmInt$g, main = "True DAG")
}
}
\keyword{models}
\keyword{graph}

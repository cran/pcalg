\name{gds}
\alias{gds}
\encoding{UTF-8}
\concept{greedy DAG search}
\concept{essential graph}
\title{Greedy DAG Search to Estimate Markov Equivalence Class of DAG}
\description{
  Estimate the observational or interventional essential graph representing the 
  Markov equivalence class of a DAG by greedily optimizing a score function in 
  the space of DAGs.  In practice, greedy search should always be done in the 
  space of equivalence classes instead of DAGs, giving the functions
  \code{\link{gies}} or \code{\link{ges}} the preference over \code{gds}.
}
\usage{
gds(p, targets, score, verbose = FALSE, ...)
}
\arguments{
  \item{p}{number of variables.}
  \item{targets}{a \code{\link{list}} of intervention targets
    (cf. details).  A list of vectors, each vector listing the vertices
    of one intervention target.}
  \item{score}{an instance of a class derived from \code{\linkS4class{Score}}.}
  \item{verbose}{if \code{TRUE}, detailed output is provided.}
  \item{\dots}{additional arguments for debugging purposes and fine tuning.}
  %% those in    getClass("EssGraph")@refMethods$causal.inf.options
}
\details{
  This function estimates the observational or interventional Markov 
  equivalence class of a DAG
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
    \item an undirected edge a -- b stands for an arrow that is oriented in one 
      way in some representatives of the equivalence class and in the other way 
      in other representatives of the equivalence class.
  }
  Note that when plotting the object, undirected and bidirected edges are 
  equivalent.
  
  Greedy DAG search (GDS) maximizes a score function (typically the BIC, passed
  to the function via the argument \code{score}) of a DAG in three phases, 
  starting from the empty DAG:
  \describe{
    \item{Forward phase}{In the forward phase, GDS adds single arrows to the 
      DAG as long as this augments the score.}
    \item{Backward phase}{In the backward phase, the algorithm removes arrows
      from the DAG as long as this augments the score.}
    \item{Turning phase}{In the turning phase, the algorithm reverts arrows of
      the DAG as long as this augments the score.}
  }
  GIES cycles through these three phases until no augmentation of the score is
  possible any more.  In the end, \code{gds} returns the (interventional or
  observational) essential graph of the last visited DAG.
  
  It is well-known that a greedy search in the space of DAGs instead of 
  essential graphs is more prone to be stuck in local optima of the score
  function and hence expected to yield worse estimation results than GIES
  (function \code{\link{gies}}) or GES (function \code{\link{ges}}) (Chickering,
  2002; Hauser and Bühlmann, 2012).  The 
  function \code{gds} is therefore not of practical use, but can be used
  to compare causal inference algorithms to an elementary and straight-forward
  approach.
}
\value{
  \code{gds} returns a list with the following two components:
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
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\seealso{
  \code{\link{gies}}, \code{\link{ges}}, \code{\linkS4class{Score}}, 
  \code{\linkS4class{EssGraph}}
}
\examples{
## Load predefined data
data(gmInt)

## Define the score (BIC)
score <- new("GaussL0penIntScore", gmInt$x, gmInt$targets, gmInt$target.index) 

## Estimate the essential graph
gds.fit <- gds(ncol(gmInt$x), gmInt$targets, score) 

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(gds.fit$essgraph, main = "Estimated ess. graph")
  plot(gmInt$g, main = "True DAG")
}
}
\keyword{models}
\keyword{graphs}

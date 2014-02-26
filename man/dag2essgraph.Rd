\name{dag2essgraph}
\title{Convert a DAG to an Essential Graph}
\alias{dag2essgraph}
\encoding{UTF-8}
\concept{DAG}
\concept{CPDAG}
\concept{essential graph}
\description{
Convert a DAG to an (interventional or observational) essential graph.
}
\usage{
dag2essgraph(dag, targets = list(integer(0)))
}
\arguments{
  \item{dag}{The DAG whose essential graph has to be calculated, represented as 
    an instance of a class derived from \code{\linkS4class{Score}}.}
  
  \item{targets}{List of intervention targets with respect to which the 
    essential graph has to be calculated.  An observational setting is
    represented by one single empty target (\code{list(integer(0))}).}
}
\details{
  This function converts a DAG to its corresponding (interventional or 
  observational) essential graph, using the algorithm of Hauser and Bühlmann
  (2012).
  
  The essential graph is a partially directed graph that represents the
  (interventional or observational) Markov equivalence class of a DAG.  It has
  the same has the same skeleton as the DAG; a directed edge represents an
  arrow that has a common orientation in all representatives of the 
  (interventional or observational) Markov equivalence class, whereas an
  undirected edge represents an arrow that has different orientations in
  different representatives of the equivalence class.  In the observational 
  case, the essential graph is also known as ``CPDAG'' (Spirtes \emph{et al.},
  2000).
  
  In a purely observational setting (\emph{i.e.}, if \code{targets = 
  list(integer(0))}), the function yields the same graph as \link{dag2cpdag},
  although it uses different classes for parameters and return value.
}
\value{
  An instance of \code{\linkS4class{EssGraph}} representing the essential graph.
}
\references{
  A. Hauser and P. Bühlmann (2012).  Characterization and greedy learning of 
  interventional Markov equivalence classes of directed acyclic graphs.
  \emph{Journal of Machine Learning Research} \bold{13}, 2409--2464.
  
  P. Spirtes, C.N. Glymour, and R. Scheines (2000).  \emph{Causation,
  Prediction, and Search}, MIT Press, Cambridge (MA).
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\seealso{
  \code{\link{dag2cpdag}}, \code{\linkS4class{Score}}, \code{\linkS4class{EssGraph}}
}
\examples{
p <- 10     # Number of random variables
s <- 0.4    # Sparseness of the DAG

## Generate a random DAG
set.seed(42)
dag <- as(randomDAG(p, s), "GaussParDAG")

## Calculate observational essential graph
res.obs <- dag2essgraph(dag)

## Calculate interventional essential graph for intervention targets
## {1} and {3}
res.int <- dag2essgraph(dag, as.list(c(1, 3)))
}
\keyword{models}
\keyword{graphs}

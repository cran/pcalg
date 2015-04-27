\name{EssGraph-class}
\docType{class}
\alias{EssGraph-class}
\alias{plot,EssGraph,ANY-method}
\encoding{UTF-8}
\title{Class \code{"EssGraph"}}
\description{
  This class represents an (observentional or interventional) essential graph.
}
\details{
  An observational or interventional Markov equivalence class of DAGs can be
  uniquely represented by a partially directed graph, the essential
  graph.  Its edges have the following interpretation:
  \enumerate{
    \item a directed edge \eqn{a \longrightarrow b}{a → b} stands for an arrow
      that has the same orientation in all representatives of the
      Markov equivalence class;
    \item an undirected edge \eqn{a - b} stands for an arrow that is oriented in one
      way in some representatives of the equivalence class and in the other way
      in other representatives of the equivalence class.
  }
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\section{Constructor}{
  \code{new("EssGraph", nodes, in.edges, ...)}
  \describe{
    \item{\code{nodes}}{Vector of node names; cf. also field \code{.nodes}.}
    \item{\code{in.edges}}{A list of length \code{p} consisting of index
      vectors indicating the edges pointing into the nodes of the DAG.}
  }
}
\section{Fields}{
  \describe{
    \item{\code{.nodes}:}{Vector of node names; defaults to \code{as.character(1:p)},
      where \code{p} denotes the number of nodes (variables) of the model.}
    \item{\code{.in.edges}:}{A list of length \code{p} consisting of index
      vectors indicating the edges pointing into the nodes of the DAG.}
    \item{\code{targets}}{List of mutually exclusive intervention targets
      with respect to which Markov equivalence is defined.}
    \item{\code{score}:}{Object of class \code{\linkS4class{Score}}; used
      internally for score-based causal inference.}
  }
}
\section{Class-Based Methods}{
  Most class-based methods are only for internal use.  Methods of interest for
  the user are:
  \describe{
    \item{\code{repr()}:}{Yields a representative causal model of the
      equivalence class, an object of a class derived from
      \code{\linkS4class{Score}}.  Since the representative is not only
      characterized by the DAG, but also by appropriate parameters, the field
      \code{score} must be assigned for this method to work.  The DAG is
      drawn at random; note that all representatives are statistically
      indistinguishable under a given set of intervention targets.}
    \item{\code{node.count()}:}{Yields the number of nodes of the essential
      graph.}
    \item{\code{edge.count()}:}{Yields the number of edges of the essential
      graph.  Note that \emph{unoriented} edges count as 2, whereas
      \emph{oriented} edges count as 1 due to the internal representation.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "EssGraph", y = "ANY")}: plots the
      essential graph.  In the plot, undirected and bidirected edges are equivalent.}
  }
}
\author{
  Alain Hauser (\email{alain.hauser@bfh.ch})
}
\seealso{
  \code{\linkS4class{ParDAG}}
}
\examples{
showClass("EssGraph")
}
\keyword{classes}

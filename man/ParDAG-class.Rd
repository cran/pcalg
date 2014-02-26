\name{ParDAG-class}
\Rdversion{1.1}
\docType{class}
\alias{ParDAG-class}
\alias{plot,ParDAG,ANY-method}
\title{Class \code{"ParDAG"}}
\description{
  This virtual base class represents a parametric causal model.
}
\details{
  The class \code{"ParDAG"} serves as a basis for simulating observational
  and/or interventional data from causal models as well as for parameter
  estimation (maximum-likelihood estimation) for a given causal model in the
  presence of a data set with jointly observational and interventional data.
  
  The virtual base class \code{"ParDAG"} provides a "skeleton" for all 
  functions relied to the aforementioned task.  In practical cases, a user may 
  always choose an appropriate class derived from \code{ParDAG} which 
  represents a specific parametric model class.  The base class itself does
  \emph{not} represent such a model class.
}
\section{Extends}{
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
\section{Constructor}{
  \code{new("ParDAG", nodes, in.edges, params)}
  \describe{
    \item{\code{nodes}}{Vector of node names; cf. also field \code{.nodes}.}
    \item{\code{in.edges}}{A list of length \code{p} consisting of index 
      vectors indicating the edges pointing into the nodes of the DAG.}
    \item{\code{params}}{A list of length \code{p} consisting of parameter 
      vectors modeling the conditional distribution of a node given its 
      parents; cf. also field \code{.params}.}
  }
}
\section{Fields}{
  \describe{
    \item{\code{.nodes}:}{Vector of node names; defaults to \code{as.character(1:p)},
      where \code{p} denotes the number of nodes (variables) of the model.}
    \item{\code{.in.edges}:}{A list of length \code{p} consisting of index 
      vectors indicating the edges pointing into the nodes of the DAG.}
    \item{\code{.params}:}{A list of length \code{p} consisting of parameter 
      vectors modeling the conditional distribution of a node given its 
      parents.  The entries of the parameter vectors only get a concrete 
      meaning in derived classes belonging to specific parametric model classes.}
  }
}
\section{Class-Based Methods}{
  \describe{
    \item{\code{node.count()}:}{Yields the number of nodes (variables) of the 
      model.}
    \item{\code{simulate(n, target, int.level)}:}{Generates \eqn{n}
      (observational or interventional) samples from the parametric causal
      model.  The intervention target to be used is specified by the parameter
      \code{target}; if the target is empty (\code{target = integer(0)}),
      observational samples are generated.  \code{int.level} indicates
      the values of the intervened variables; if it is a vector of the same
      length as \code{target}, all samples are drawn from the same intervention
      levels; if it is a matrix with \eqn{n} rows and as many columns as 
      \code{target} has entries, its rows are interpreted as individual
      intervention levels for each sample.}
    \item{\code{edge.count()}:}{Yields the number of edges (arrows) in the DAG.}
    \item{\code{mle.fit(score)}:}{Fits the parameters using an appropriate
      \code{\linkS4class{Score}} object.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "ParDAG", y = "ANY")}: plots the underlying 
    DAG of the causal model.  Parameters are not visualized.}
  }
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\seealso{
  \code{\linkS4class{GaussParDAG}}
}
\keyword{classes}


\name{Score-class}
\docType{class}
\alias{Score-class}
\encoding{UTF-8}
\title{Virtual Class "Score"}
\description{
  This virtual base class represents a score for causal inference; it is used
  in the causal inference functions \code{\link{ges}}, \code{\link{gies}} and
  \code{\link{simy}}. 
}
\section{Fields}{
  The fields of \code{Score} are mainly of interest for users who aim at 
  deriving an own class from this virtual base class, i.e., implementing an own
  score function.
  \describe{
    \item{\code{.nodes}:}{Node labels.  They are passed to causal inference
      methods by default to label the nodes of the resulting graph.}
    \item{\code{decomp}:}{Indicates whether the represented score is
      decomposable (cf. details).  At the moment, only decomposable scores are 
      supported by the implementation of the causal inference algorithms; 
      support for non-decomposable scores is planned.}
    \item{\code{pp.dat}:}{List representing the preprocessed input data;
      this is typically a statistic which is sufficient for the 
      calculation of the score.}
    \item{\code{.pardag.class}:}{Name of the class of the parametric DAG model
      corresponding to the score.  This must name a class derived from 
      \code{ParDAG}.}
    \item{\code{c.fcn}:}{Only used internally; must remain empty for (user
      specified) classes derived from \code{Score}.} 
  }
}
\section{Constructor}{
\preformatted{
new("Score",
  data = matrix(1, 1, 1), 
  targets = list(integer(0)),
  target.index = rep(as.integer(1), nrow(data)), 
  nodes = colnames(data),
  ...)
}
  \describe{
    \item{\code{data}}{Data matrix with \eqn{n} rows and \eqn{p} columns.  Each 
      row corresponds to one realization, either interventional or 
      observational.}
    \item{\code{targets}}{List of mutually exclusive intervention targets
      that have been used for data generation.}
    \item{\code{target.index}}{Vector of length \eqn{n}; the \eqn{i}-th entry 
      specifies the index of the intervention
      target in \code{targets} under which the \eqn{i}-th row of \code{data}
      was measured.}
    \item{\code{nodes}}{Node labels}
    \item{\code{...}}{Additional parameters used by derived (and non-virtual)
      classes.}
  }
}
\section{Methods}{
  Note that since \code{Score} is a virtual class, its methods cannot be called
  directly, but only on derived classes.
  \describe{
    \item{\code{local.score(vertex, parents, ...)}}{For decomposable scores,
      this function calculates the local score of a vertex and its parents.
      Must throw an error in derived classes that do not represent a 
      decomposable score.}
    \item{\code{global.score.int(edges, ...)}}{Calculates the global score of
      a DAG, represented as a list of in-edges: for each vertex in the DAG,
      this list contains a vector of parents.}
    \item{\code{global.score(dag, ...)}}{Calculates the global score of a DAG,
      represented as object of a class derived from \code{ParDAG}.}
    \item{\code{local.fit(vertex, parents, ...)}}{Calculates a local model fit
      of a vertex and its parents, e.g. by MLE. 
      The result is a vector of parameters whose 
      meaning depends on the model class; it matches the convention used in the
      corresponding causal model (cf. \code{.pardag.class}).
      % TODO: adjust this!!
      }
    \item{\code{global.fit(dag, ...)}}{Calculates the global MLE of a DAG,
      represented by an object of the class specified by \code{.pardag.class}.
      The result is a list of vectors, one per vertex, each in the same format
      as the result vector of \code{local.mle}.
      % TODO: adjust this!!
      }
  }
}
\details{
  Score-based structure learning algorithms for causal inference such as 
  Greedy Equivalence Search
  (GES, implemented in the function \code{\link{ges}}), Greedy Interventional
  Equivalence Search (GIES, implemented in the function \code{\link{gies}}) and
  the dynamic programming approach of Silander and Myllymäki (2006)
  (implemented in the function \code{\link{simy}}) try to find the DAG model which
  maximizes a scoring criterion for a given data set.  A widely-used scoring 
  criterion is the Bayesian Information Criterion (BIC).
  
  The virtual class \code{Score} is the base class for providing a scoring
  criterion to the mentioned causal inference algorithms.  It does not 
  implement a concrete scoring criterion, but it defines the functions that 
  must be provided by its descendants (cf. methods).
  
  Knowledge of this class is only required if you aim to implement an own 
  scoring criterion.  At the moment, it is recommended to use the predefined
  scoring criteria for multivariate Gaussian data derived from \code{Score}, 
  \code{\linkS4class{GaussL0penIntScore}} and
  \code{\linkS4class{GaussL0penObsScore}}. 
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
  \code{\link{ges}}, \code{\link{gies}}, \code{\link{simy}}, 
  \code{\linkS4class{GaussL0penIntScore}}, 
  \code{\linkS4class{GaussL0penObsScore}}
}
\keyword{classes}

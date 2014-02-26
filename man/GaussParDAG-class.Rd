\name{GaussParDAG-class}
\title{Class \code{"GaussParDAG"}}
\alias{GaussParDAG-class}
\Rdversion{1.1}
\docType{class}
\encoding{UTF-8}
\concept{Gaussian causal model}
\concept{structural equation model}
\description{
  The \code{"GaussParDAG"} class represents a Gaussian causal model.
}
\details{
  The class \code{"GaussParDAG"} is used to simulate observational
  and/or interventional data from Gaussian causal models as well as for parameter
  estimation (maximum-likelihood estimation) for a given DAG structure in the
  presence of a data set with jointly observational and interventional data.
  
  A Gaussian causal model can be represented as a set of \eqn{p} linear 
  structural equations with Gaussian noise variables.  Those equations are 
  fully specified by indicating the regression parameters, the intercept
  and the variance of the noise or error terms.  More details can be found e.g.
  in Kalisch and Bühlmann (2007) or Hauser and Bühlmann (2012).
}
\section{Extends}{
  Class \code{"\linkS4class{ParDAG}"}, directly.

  All reference classes extend and inherit methods from 
  \code{"\linkS4class{envRefClass}"}.
}
\section{Constructor}{
  \code{new("GaussParDAG", nodes, in.edges, params)}
  \describe{
    \item{\code{nodes}}{Vector of node names; cf. also field \code{.nodes}.}
    \item{\code{in.edges}}{A list of length \code{p} consisting of index 
      vectors indicating the edges pointing into the nodes of the DAG.}
    \item{\code{params}}{A list of length \code{p} consisting of parameter 
      vectors modeling the conditional distribution of a node given its 
      parents; cf. also field \code{.params} for the meaning of the 
      parameters.}
  }
}
\section{Fields}{
  \describe{
    \item{\code{.nodes}:}{Vector of node names; defaults to \code{as.character(1:p)},
      where \code{p} denotes the number of nodes (variables) of the model.}
    \item{\code{.in.edges}:}{A list of length \code{p} consisting of index 
      vectors indicating the edges pointing into the nodes of the DAG.  The 
      \eqn{i}-th entry lists the indices of the parents of the \eqn{i}-th node.}
    \item{\code{.params}:}{A list of length \code{p} consisting of parameter 
      vectors modeling the conditional distribution of a node given its 
      parents.  The \eqn{i}-th entry models the conditional (normal) 
      distribution of the \eqn{i}-th variable in the model given its parents.
      It is a vector of length \eqn{k + 2}, where \eqn{k} is the number of 
      parents of node \eqn{i}; the first entry encodes the error variance of
      node \eqn{i}, the second entry the intercept, and the remaining entries
      the regression coefficients (see above).  In most cases, it is easier
      to access the parameters via the wrapper functions \code{err.var}, 
      \code{intercept} and \code{weight.mat}.}
  }
}
\section{Class-Based Methods}{
  \describe{
    \item{\code{set.err.var(value)}:}{Sets the error variances.  The argument
      must be a vector of length \eqn{p}, where \eqn{p} denotes the number
      of nodes in the model.}
    \item{\code{err.var()}:}{Yields the vector of error variances.}
    \item{\code{intercept()}:}{Yields the vector of intercepts.}
    \item{\code{set.intercept(value)}:}{Sets the intercepts.  The argument
      must be a vector of length \eqn{p}, where \eqn{p} denotes the number
      of nodes in the model.}
    \item{\code{weight.mat(target)}:}{Yields the (observational or 
      interventional) weight matrix of the model. The weight matrix is an 
      \eqn{p \times p}{p x p} matrix whose \eqn{i}-th columns contains the 
      regression coefficients of the \eqn{i}-th structural equation, if node
      \eqn{i} is not intervened (i.e., if \code{i} is not contained in the 
      vector \code{target}), and is empty otherwise.}
    \item{\code{cov.mat(target, ivent.var)}:}{Yields the covariance matrix
      of the observational or an interventional distribution of the causal
      model.  If \code{target} has length 0, the covariance matrix of the 
      observational distribution is returned; otherwise \code{target} is a 
      vector of the intervened nodes, and \code{ivent.var} is a vector of the
      same length indicating the variances of the intervention levels.
      Deterministic interventions with fix intervention levels would correspond
      to vanishing intervention variances; with non-zero intervention variances,
      stochastic interventions are considered in which intervention values are
      realizations of Gaussian variables (Korb et al., 2004).}
  }

  The following methods are inherited (from the corresponding class):
  node.count ("ParDAG"), edge.count ("ParDAG"), simulate ("ParDAG")
}
\references{
  A. Hauser and P. Bühlmann (2012).  Characterization and greedy learning of 
  interventional Markov equivalence classes of directed acyclic graphs.
  \emph{Journal of Machine Learning Research} \bold{13}, 2409--2464.
  
  M. Kalisch and P. Buehlmann (2007).  Estimating high-dimensional directed 
  acyclic graphs with the PC-algorithm. \emph{Journal of Machine Learning 
  Research} \bold{8}, 613--636.
  
  K.B. Korb, L.R. Hope, A.E. Nicholson, and K. Axnick (2004).  Varieties of 
  causal intervention.  \emph{Proc. of the Pacific Rim International Conference 
  on Artificial Intelligence (PRICAI 2004)}, 322--331
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\seealso{
  \code{\linkS4class{ParDAG}}
}
\examples{
set.seed(307)
myDAG <- r.gauss.pardag(p = 5, prob = 0.4)
myDAG$weight.mat()
myDAG$err.var()
myDAG$intercept()
myDAG$set.intercept(runif(5, 3, 4))
myDAG$intercept()
if (require(Rgraphviz)) plot(myDAG)
}
\keyword{classes}
\keyword{models}

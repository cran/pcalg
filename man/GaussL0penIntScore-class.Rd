\name{GaussL0penIntScore-class}
\docType{class}
\alias{GaussL0penIntScore-class}
\alias{local.score,GaussL0penIntScore-method}
\alias{global.score,GaussL0penIntScore-method}
\alias{local.mle,GaussL0penIntScore-method}
\alias{global.mle,GaussL0penIntScore-method}
\encoding{UTF-8}
\title{Class \code{"GaussL0penIntScore"}}
\description{
  This class represents a score for causal inference from jointly interventional
  and observational Gaussian data; it is used in the causal inference functions 
  \code{\link{gies}} and \code{\link{simy}}. 
}
\details{
  The class implements an \eqn{\ell_0}{ℓ0}-penalized Gaussian maximum 
  likelihood estimator.  The penalization is a constant (specified by
  the argument \code{lambda} in the constructor) times the number of 
  parameters of the DAG model.  By default, the constant \eqn{\lambda}{λ} is 
  chosen as \eqn{\log(n)/2}{log(n)/2}, which corresponds to the BIC score.
}
\section{Extends}{
  Class \code{"\linkS4class{Score}"}, directly.
  
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
%% TODO:
%% - replace this by \section{Constructor}; see the documentation of graphNEL (package graph)
%% - include methods and slots, see documentation skeleton produced by promptClass("GaussL0penIntScore")
%% - include aliases for methods and slots, see documentation of graphNEL
\section{Fields}{
  The class \code{GaussL0penIntScore} has the same fields as \code{\linkS4class{Score}}.   
  They need not be accessed by the user.
}
\section{Constructor}{
\preformatted{
new("GaussL0penIntScore",
  data = matrix(1, 1, 1),
  targets = list(integer(0)), 
  target.index = rep(as.integer(1), nrow(data)), 
  lambda = 0.5*log(nrow(data)), 
  intercept = FALSE, 
  use.cpp = TRUE, 
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
    \item{\code{lambda}}{Penalization constant (cf. details)}
    \item{\code{intercept}}{Indicates whether an intercept is allowed in the linear
      structural equations, or, equivalently, if a mean different from zero is
      allowed for the observational distribution.}
    \item{\code{use.cpp}}{Indicates whether the calculation of the score should be done
      by the C++ library of the package, which speeds up calculation.  This
      parameter should only be set to \code{FALSE} in the case of problems.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{local.score(vertex, parents, ...)}}{Calculates the local score 
      of a vertex and its parents.  Since this score has no obvious 
      interpretation, it is rather for internal use.}
    \item{\code{global.score.int(edges, ...)}}{Calculates the global score of
      a DAG, represented as a list of in-edges: for each vertex in the DAG,
      this list contains a vector of parents.}
    \item{\code{global.score(dag, ...)}}{Calculates the global score of a DAG,
      represented as an object of a class derived from
      \code{\linkS4class{ParDAG}}.} 
    \item{\code{local.mle(vertex, parents, ...)}}{Calculates the local MLE of 
      a vertex and its parents.  The result is a vector of parameters encoded
      as follows:
      \itemize{
        \item First element: variance of the Gaussian error term
        \item Second element: intercept
        \item Following elements: regression coefficients; one per parent vertex
      }
      }
    \item{\code{global.mle(dag, ...)}}{Calculates the global MLE of a DAG,
      represented by an object of a class derived from
      \code{\linkS4class{ParDAG}}. 
      The result is a list of vectors, one per vertex, each in the same format
      as the result vector of \code{local.mle}.}
  }
}
\author{
  Alain Hauser (\email{alain.hauser@biology.unibe.ch})
}
\examples{
##################################################
## Using Gaussian Data
##################################################
## Load predefined data
data(gmInt)

## Define the score object
score <- new("GaussL0penIntScore", gmInt$x, gmInt$targets, gmInt$target.index) 

## Score of the true underlying DAG
score$global.score(as(gmInt$g, "GaussParDAG"))

## Score of the DAG that has only one edge from 1 to 2
A <- matrix(0, ncol(gmInt$x), ncol(gmInt$x))
A[1, 2] <- 1
score$global.score(as(A, "GaussParDAG"))
## (Note: this is lower than the score of the true DAG.)
}
\seealso{
  \code{\link{gies}}, \code{\link{simy}}, 
  \code{\linkS4class{GaussL0penObsScore}}, 
  \code{\linkS4class{Score}}
}
\keyword{classes}

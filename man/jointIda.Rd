\name{jointIda}
\alias{jointIda}
\title{Estimate Multiset of Possible Total Joint Effects}

\description{
  \code{jointIda()} estimates the multiset of possible total joint effects
  of a set of intervention variables (\code{X}) on another variable (\code{Y})
  from observational data.  This is a version of \code{\link{ida}} that
  allows multiple simultaneous interventions.
}
\usage{
jointIda(x.pos, y.pos, mcov, graphEst = NULL, all.pasets = NULL,
         technique = c("RRC", "MCD"), type = c("pdag", "cpdag", "dag"))
}

\arguments{
  \item{x.pos}{(integer vector) positions of the intervention variables
    \code{X} in the covariance matrix.}
  \item{y.pos}{(integer) position of variable \code{Y} in the covariance
    matrix. (\code{y.pos} can also be an integer vector, see Note.)}
  \item{mcov}{(estimated) covariance matrix.}
  \item{graphEst}{(graphNEL object) Estimated CPDAG or PDAG. The CPDAG is typically from \code{\link{pc}()}: If
    the result of \code{\link{pc}} is \code{pc.fit}, the estimated CPDAG
    can be obtained by \code{pc.fit@graph}. The PDAG can be obtained from CPDAG by adding bacground knowledge using \code{addBgKnowledge()}.  \code{graphEst} can only be considered if
    \code{all.pasets} is \code{NULL}.}
  \item{all.pasets}{(an optional argument and the default is
    \code{NULL}) A list where each element is a list of size
    \code{length(x.pos)}.  Each sub-list \code{all.pasets[[i]]} contains
    possible parent sets of \code{x.pos} in the same order, i.e.,
    \code{all.pasets[[i]][[j]]} is a possible parent set of
    \code{x.pos[j]}.  This option can be used if possible parent sets of
    the intervention variables are known.}
  \item{technique}{character string specifying the technique that will
    be used to estimate the total joint causal effects (given the parent
    sets), see details below.
    \describe{
      \item{\code{"RRC"}:}{Recursive regressions for causal effects.}
      \item{\code{"MCD"}:}{Modifying the Cholesky decomposition.}
    }
  }
   \item{type}{Type of graph \code{"graphEst"}; can be of type \code{"cpdag"}, \code{"pdag"} (e.g. a CPDAG with background knowledge from Meek, 1995) or \code{"dag"}.}
}

\details{
  It is assumed that we have observational data that are multivariate
  Gaussian and faithful to the true (but unknown) underlying causal DAG
  (without hidden variables).  Under these assumptions, this function
  estimates the multiset of possible total joint effects of \code{X} on
  \code{Y}. Here the total joint effect of \eqn{X = (X_1,X_2)} on
  \eqn{Y} is defined via Pearl's do-calculus as the vector
  \eqn{(E[Y|do(X_1=x_1+1,X_2=x_2)]-E[Y|do(X_1=x_1,X_2=x_2)], E[Y|do(X_1=x_1,X_2=x_2+1)]-E[Y|do(X_1=x_1,X_2=x_2)])},
  with a similar definition for more than two variables.  These values
  are equal to the partial derivatives (evaluated at \eqn{(x_1,x_2)}) of
  \eqn{E[Y|do(X=x_1',X_2=x_2')]} with respect to \eqn{x_1'} and
  \eqn{x_2'}.  Moreover, under the Gaussian assumption, these partial
  derivatives do not depend on the values at which they are evaluated.

  We estimate a \emph{multiset} of possible total joint effects instead of
  the unique total joint effect, since it is typically impossible to
  identify the latter when the true underlying causal DAG is unknown
  (even with an infinite amount of data).  
  
  Conceptually, the method
  works as follows.  First, we estimate the CPDAG or PDAG based on the data.
  The CPDAG represents the equivalence class of DAGs and can be estimated 
  from observational data with the function \code{\link{pc}} (see the help file of this function). 
  
  The PDAG contains more orientations than the CPDAG and thus, represents a 
  smaller equivalence class of DAGs, compared to the CPDAG. We can obtain a PDAG if
  we have background knowledge of, for example, certain edge orientations of
  undirected edges in the CPDAG. We obtain the PDAG by adding these
  orientations to the CPDAG using the function \code{\link{addBgKnowledge}} 
  (see the help file of this function). 
  
  Then using the CPDAG or PDAG we extract a collection of  "jointly valid" parent sets of the intervention variables from the
  estimated CPDAG. For each set of jointly valid parent sets we apply
  RRC (recursive regressions for causal effects) or MCD (modifying the
  Cholesky decomposition) to estimate the total joint effect of \code{X}
  on \code{Y} from the sample covariance matrix (see Section 3 of Nandy et. al, 2015).
}

\value{
  A matrix representing the multiset containing the estimated
  possible total joint effects of \code{X} on \code{Y}.  The number of
  rows is equal to \code{length(x.pos)}, i.e., each column represents a
  vector of possible joint causal effects.
}
\references{
  P. Nandy, M.H. Maathuis and T.S. Richardson (2017).
  Estimating the effect of joint interventions from observational data
  in sparse high-dimensional settings. In \emph{Annals of Statistics}.

  E. Perkovic, M. Kalisch and M.H. Maathuis (2017). Interpreting and using 
  CPDAGs with background knowledge. In \emph{Proceedings of UAI 2017}.
}
\author{
Preetam Nandy, Emilija Perkovic
}
\note{
  For a single variable \code{X}, \code{jointIda()} estimates the
  same quantities as \code{ida()}. 
  If \code{graphEst} is of \code{type = "cpdag"}, \code{jointIda()} obtains 
  \code{all.pasets} by using the semi-local approach described in Section 5 in Nandy et. al, (2015).
  Nandy et. al, (2015) show that \code{jointIda()} yields 
  correct multiplicities of the distinct elements of the resulting multiset (in the sense that it matches
  \code{ida()} with \code{method="global"} up to a constant factor).
  
  If \code{graphEst} is of \code{type = "pdag"},  \code{jointIda()} obtains 
  \code{all.pasets} by using the semi-local approach described in Algorithm 2, 
  Section 4.2 in Perkovic et. al (2017). For this case, \code{jointIda()} does not necessarily yield 
  the correct multiplicities of the distinct elements of the resulting multiset (it behaves similarly to
  \code{ida()} with \code{method="local"}).

  \code{jointIda()} (like \code{\link{idaFast}}) also allows direct
  computation of the total joint effect of a set of intervention
  variables \code{X} on another set of target variables \code{Y}. In
  this case, \code{y.pos} must be an integer vector containing positions
  of the target variables \code{Y} in the covariance matrix and the
  output is a list of matrices that correspond to the variables in
  \code{Y} in the same order. This method is slightly more efficient
  than looping over \code{jointIda()} with single target variables, if
  \code{all.pasets} is not specified.
}

\seealso{
  \code{\link{ida}}, the simple version;
  \code{\link{pc}} for estimating a CPDAG.
}
\examples{
## Create a weighted DAG
p <- 6
V <- as.character(1:p)
edL <- list(
  "1" = list(edges=c(3,4), weights=c(1.1,0.3)),
  "2" = list(edges=c(6),  weights=c(0.4)),
  "3" = list(edges=c(2,4,6),weights=c(0.6,0.8,0.9)),
  "4" = list(edges=c(2),weights=c(0.5)),
  "5" = list(edges=c(1,4),weights=c(0.2,0.7)),
  "6" = NULL)
myDAG <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed") ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
myPDAG <- addBgKnowledge(myCPDAG,1,3) ## true PDAG with background knowledge 1 -> 3
covTrue <- trueCov(myDAG) ## true covariance matrix

n <- 1000
## simulate Gaussian data from the true DAG
dat <- if (require("mvtnorm")) {
  set.seed(123)
  rmvnorm(n, mean=rep(0,p), sigma=covTrue)
} else readRDS(system.file(package="pcalg", "external", "N_6_1000.rds"))

## estimate CPDAG and PDAG -- see  help(pc), help(addBgKnoweldge)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p = p, alpha = 0.01, u2pd="relaxed")
pc.fit.pdag <- addBgKnowledge(pc.fit@graph,1,3)

if (require(Rgraphviz)) {
  ## plot the true and estimated graphs
  par(mfrow = c(1,3))
  plot(myDAG, main = "True DAG")
  plot(pc.fit, main = "Estimated CPDAG")
  plot(pc.fit.pdag, main = "Estimated PDAG")
}

## Suppose that we know the true CPDAG and covariance matrix
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="RRC", type = "cpdag")
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="MCD", type = "cpdag")

## Suppose that we know the true PDAG and covariance matrix
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myPDAG, technique="RRC", type = "pdag")
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myPDAG, technique="MCD", type = "pdag")

## Instead of knowing the true CPDAG or PDAG, it is enough to know only
## the jointly valid parent sets of the intervention variables
## to use RRC or MCD
## all.jointly.valid.pasets:
ajv.pasets <- list(list(5,c(3,4)),list(integer(0),c(3,4)),list(3,c(3,4)))
jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="RRC")
jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="MCD")

## From the true DAG, we can compute the true total joint effects
## using RRC or MCD
cat("Dim covTrue: ", dim(covTrue),"\n")
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myDAG, technique="RRC", type = "dag")
jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myDAG, technique="MCD", type = "dag")


## When working with data, we have to use the estimated CPDAG or PDAG
## and the sample covariance matrix
jointIda(x.pos=c(1,2), y.pos=6, cov(dat), graphEst=pc.fit@graph, technique="RRC", type = "cpdag")
jointIda(x.pos=c(1,2), y.pos=6, cov(dat), graphEst=pc.fit@graph, technique="MCD", type = "cpdag")

jointIda(x.pos=c(1,2), y.pos=6, cov(dat), graphEst=pc.fit.pdag, technique="RRC", type = "pdag")
jointIda(x.pos=c(1,2), y.pos=6, cov(dat), graphEst=pc.fit.pdag, technique="MCD", type = "pdag")
## RRC and MCD can produce different results when working with data

## jointIda also works when x.pos has length 1 and in the following example
## it gives the same result as ida() (see Note)
##
## When the CPDAG is known
jointIda(x.pos=1, y.pos=6, covTrue, graphEst=myCPDAG, technique="RRC", type = "cpdag")
ida(x.pos=1, y.pos=6, covTrue, graphEst=myCPDAG, method="global", type = "cpdag")

## When the PDAG is known
jointIda(x.pos=1, y.pos=6, covTrue, graphEst=myPDAG, technique="RRC", type = "pdag")
ida(x.pos=1, y.pos=6, covTrue, graphEst=myPDAG, method="global", type = "pdag")

## When the DAG is known
jointIda(x.pos=1, y.pos=6, covTrue, graphEst=myDAG, technique="RRC", type = "dag")
ida(x.pos=1, y.pos=6, covTrue, graphEst=myDAG, method="global")
## Note that, causalEffect(myDAG,y=6,x=1) does not give the correct value in this case,
## since this function requires that the variables are in a causal order.

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

\name{ida}
\title{Estimate Multiset of Possible Joint Total Causal Effects}
\alias{ida}
\alias{causalEffect}
\description{
  \code{ida()} estimates the multiset of possible joint total causal effects
  of variables (\code{X}) onto variables (\code{Y})
  from observational data via adjustment.
}
\usage{
ida(x.pos, y.pos, mcov, graphEst, method = c("local","optimal","global"),
    y.notparent = FALSE, verbose = FALSE, all.dags = NA, type = c("cpdag", "pdag"))
}
\arguments{
  \item{x.pos, x}{Positions of variables \code{X} in the covariance matrix.}
  \item{y.pos, y}{Positions of variables \code{Y} in the covariance matrix.}
  \item{mcov}{Covariance matrix that was used to estimate \code{graphEst}.}
  \item{graphEst}{Estimated CPDAG or PDAG. The CPDAG is typically from \code{\link{pc}()}: If
    the result of \code{\link{pc}} is \code{pc.fit}, the estimated CPDAG
    can be obtained by \code{pc.fit@graph}. A PDAG can be obtained from the CPDAG by adding background knowledge using \code{addBgKnowledge()}.}
  \item{method}{Character string specifying the method with default \code{"local"}.
    \describe{
      \item{\code{"global"}:}{The algorithm considers all undirected edges 
      in the CPDAG or PDAG, using the possible parents as adjustment sets 
      to estimate \emph{all} possible causal effects. It is hence \emph{slow}
      and can only be applied to singleton \code{X}.} 
      \item{\code{"local"}:}{The algorithm only considers edges in the
      neighborhood of \code{X} in the CPDAG or PDAG, 
      using the possible parents as adjustment sets 
      to estimate the \emph{unique} possible causal effects. 
      It is hence much \emph{faster} than \code{"global"} 
      and can only be applied to singleton \code{X}.}
      \item{\code{"optimal"}:}{The algorithm considers only those edges necessary to compute
      the possible optimal valid adjustment sets, using these as adjustment sets to estimate 
      the \emph{unique} possible causal effects. It is
      hence \emph{faster} than the \code{"global"} option but also slower than \code{"local"}. 
      It provides more efficient estimates than both alternatives but is also more sensitive to faulty graph estimates. 
      Can be applied to sets \code{X}.}
    }
    See details below.}
  \item{y.notparent}{Logical; for singleton \code{X} and \code{Y}. 
    If true, any edge between \code{X} and
    \code{Y} is assumed to be of the form \code{X->Y}. Not implemented for the \code{method="optimal"}}
  \item{verbose}{If TRUE, details on the regressions are printed.}
  \item{all.dags}{All DAGs in the equivalence class represented by the CPDAG or PDAG 
     can be precomputed by \code{\link{pdag2allDags}()} and passed via
    this argument.  In that case, \code{\link{pdag2allDags}(..)} is not called
    internally.
    This option is only relevant when using \code{method="global"}.}
  \item{type}{Type of graph \code{"graphEst"}; can be of type \code{"cpdag"} or \code{"pdag"} (e.g. a CPDAG with background knowledge from Meek, 1995)} 
}
\details{
  It is assumed that we have observational data from a multivariate Gaussian distribution
  faithful to the true (but unknown) underlying causal DAG (without hidden variables). 
  Under these assumptions, this function estimates the multiset of possible total joint effects of \code{X} on \code{Y}. 
  Here the total joint effect of \code{X} \eqn{= (X_1,X_2)} on \code{Y} is defined via Pearl's do-calculus as the vector 
  \eqn{(E[Y|do(X_1=x_1+1,X_2=x_2)]-E[Y|do(X_1=x_1,X_2=x_2)], E[Y|do(X_1=x_1,X_2=x_2+1)]-E[Y|do(X_1=x_1,X_2=x_2)])}, 
  with a similar definition for more than two variables. These values are equal to the partial derivatives 
  (evaluated at \eqn{x_1,x_2}) of \eqn{E[Y|do(X=x_1',X_2=x_2')]} with respect to \eqn{x_1}' and \eqn{x_2}'. 
  Moreover, under the Gaussian assumption, these partial derivatives do not depend on the values at which they are evaluated.


  We estimate a \emph{set} of possible joint total causal effects instead of
  the unique joint total causal effect, since it is typically impossible to
  identify the latter when the true underlying causal DAG is unknown
  (even with an infinite amount of data).  Conceptually, the method
  works as follows.  First, we estimate the equivalence class of DAGs
  that describe the conditional independence relationships in the data,
  using the function \code{\link{pc}} (see the help file of this
  function). For each DAG G in the equivalence class, we apply Pearl's
  do-calculus to estimate the total causal effect of \code{X} on
  \code{Y}.  This can be done via a simple linear regression
  adjusting for a valid adjustment set. 
  
  
  For example, if \code{X} and \code{Y} are singleton and \code{Y} 
  is not a parent of \code{X}, we can take the regression coefficient of
  \code{X} in the regression \code{lm(Y ~ X + pa(X,G))}, where
  \code{pa(X,G)} denotes the parents of \code{X} in the DAG G; if \code{Y}
  is a parent of \code{X} in G, we can set the estimated causal effect to
  zero. 

  If the equivalence class contains \code{k} DAGs, this will yield
  \code{k} estimated total causal effects.  Since we do not know which DAG
  is the true causal DAG, we do not know which estimated possible total joint causal
  effect of \code{X} on \code{Y} is the correct one.  Therefore, we return
  the entire multiset of \code{k} estimated effects (it is a multiset
  rather than a set because it can contain duplicate values).

  One can take summary measures of the multiset.  For example, the
  minimum absolute value provides a lower bound on the size of the true
  causal effect: If the minimum absolute value of all values in the
  multiset is larger than one, then we know that the size of the true
  causal effect (up to sampling error) must be larger than one.

  If \code{method="global"}, the method as described above is carried
  out, where all DAGs in the equivalene class of the estimated CPDAG or PDAG
  \code{graphEst} are computed using the function \code{\link{pdag2allDags}}.
  The parent set for each DAG is then used to estimate the corresponding possible
  total causal effect. This method is suitable for small graphs (say, up to 10 nodes) and can
  only be used for singleton \code{X}.
  
  If \code{method="local"}, we only consider all valid possible directions of undirected edges 
  that have \code{X} as an endpoint.
  
  In the case of a CPDAG, we consider all
  possible directions of undirected edges that have \code{X} as an
  endpoint, such that no new v-structure is created. 
  Maathuis, Kalisch and Buehlmann (2009) showed that there is at least one DAG in
  the equivalence class for each such local configuration. Hence, the 
  procedure is truly local in this setting.
  
  In the case of a PDAG, we need to verify for all possible directions whether they lead to an 
  amenable max. PDAG if we apply Meek's orientation rules. 
  In this setting the complexity of the \code{"local"} method is similar to the \code{"optimal"} one and it is not truly local.
  For details see Section 4.2 in Perkovic, Kalisch and Maathuis (2017).  
  
  We estimate the total causal effect of \code{X} on \code{Y} for each
  valid configuration as above, using linear regression adjusting for the correspoding possible parents.
  As we adjust for the same sets as in the \code{"global"} method, it follows that the multisets of total causal effects of
  the two methods have the same unique values. They may, however, have different multiplicities.

  Since the parents of \code{X} are usually an inefficient valid adjustment set we provide a third method, that uses 
  different adjustment sets.  

  If \code{method="optimal"}, we do not determine all DAGs in the
  equivalence class of the CPDAG or PDAG.  Instead, we only direct edges until
  obtaining an amenable PDAG, which is sufficient for computing the optimal 
  valid adjustment set. Each amenable PDAG can be obtained by 
  orienting the neighborhood of \code{X} and then applying Meek's orientation rules, similar to the \code{"local"} method 
  for PDAGs. This can be done faster than the \code{"global"} method but is slower than the \code{"local"} method, especially for    CPDAGs. For details see Witte, Henckel, Maathuis and Didelez (2019).
  
  For each amenable PDAG the corresponding optimal valid adjustment set is computed. 
  The optimal set is a valid adjustment set irrespectively of whether \code{X} is a singleton.
  Hence, as opposed to the other two, this method can be applied to sets \code{X}. Sometimes, however, 
  a joint total causal effect cannot be estimated via adjustment. In these cases we recommend use of the pcalg function \code{\link{jointIda}}.
  
  We then estimate the joint total causal effect of \code{X} on \code{Y} for each
  valid configuration with linear regression, adjusting for the possible optimal sets. 
  If the estimated graph is correct, each of these regressions is guaranteed 
  to be more efficient than the corresponding linear regression with any other valid adjustment set
  (see Henckel, Perkovic and Maathuis (2019) for more details). The estimates are, however, more sensitive to graph estimation errors than the ones obtained with the other two methods. 
  If \code{X} is a singleton, the output of this method is a multiset of the same size as the output of the \code{"local"} method. 
  
  

  For example, a CPDAG may represent eight DAGs, and the \code{"global"} method
  may produce an estimate of the multiset of possible total effects
  \{1.3, -0.5, 0.7, 1.3, 1.3, -0.5, 0.7, 0.7\}.
  The unique values in this set are -0.5, 0.7 and 1.3, and the
  multiplicities are 2, 3 and 3.  The \code{"local"} and \code{"optimal"} methods, on the other hand,
  may prodcue estimates of the set \{1.3, -0.5, -0.5, 0.7\}.  The unique values are again -0.5,
  0.7 and 1.3, but the multiplicities are now 2, 1 and 1.  The fact that
  the unique values of the multisets for all three methods
  are identical implies that summary measures of the multiset that only
  depend on the unique values (such as the minimum absolute value) can
  be estimated with all three.

}
\value{
  A list of length |\code{Y}| of matrices, each containing the possible joint total causal effect of \code{X} on one node in \code{Y}.
}
\references{
  M.H. Maathuis, M. Kalisch, P. Buehlmann (2009).
  Estimating high-dimensional intervention effects from observational data.
  \emph{Annals of Statistics} \bold{37}, 3133--3164.

  M.H. Maathuis, D. Colombo, M. Kalisch, P. Bühlmann (2010).
  Predicting causal effects in large-scale systems from observational data.
  \emph{Nature Methods} \bold{7}, 247--248.
  
  C. Meek (1995). Causal inference and causal explanation with background knowledge,
  In \emph{Proceedings of UAI 1995}, 403-410.

  Markus Kalisch, Martin Maechler, Diego Colombo, Marloes H. Maathuis,
  Peter Buehlmann (2012).
  Causal inference using graphical models with the R-package pcalg.
  \emph{Journal of Statistical Software} \bold{47}(11) 1--26,
  \url{http://www.jstatsoft.org/v47/i11/}.

  Pearl (2005). \emph{Causality. Models, reasoning and inference}.
  Cambridge University Press, New York.
  
  E. Perkovic, M. Kalisch and M.H. Maathuis (2017). Interpreting and using 
  CPDAGs with background knowledge. In \emph{Proceedings of UAI 2017.}
    
  L. Henckel, E. Perkovic and M.H. Maathuis (2019). Graphical criteria for efficient total 
  effect estimation via adjustment in causal linear models.
  \emph{Working Paper.}
    
  J. Witte, L. Henckel, M.H Maathuis and V. Didelez (2019). On efficient adjustment in causal graphs. 
  \emph{Working Paper.}
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch}), Emilija Perkovic and Leonard Henckel}
\seealso{
  \code{\link{jointIda}} for estimating the multiset of possible total
  \emph{joint} effects; \code{\link{idaFast}} for faster estimation of the multiset of
  possible total causal effects for several target variables. 

  \code{\link{pc}} for estimating a CPDAG. \code{\link{addBgKnowledge}} for obtaining a PDAG from CPDAG and background knowledge.
}
\examples{
## Simulate the true DAG
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
p <- 10
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
myPDAG <- addBgKnowledge(myCPDAG,2,3) ## true PDAG with background knowledge 2 -> 3
covTrue <- trueCov(myDAG) ## true covariance matrix

## simulate Gaussian data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)

## estimate CPDAG and PDAG -- see  help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01)
pc.fit.pdag <- addBgKnowledge(pc.fit@graph,2,3)

if (require(Rgraphviz)) {
  ## plot the true and estimated graphs
  par(mfrow = c(1,3))
  plot(myDAG, main = "True DAG")
  plot(pc.fit, main = "Estimated CPDAG")
  plot(pc.fit.pdag, main = "Max. PDAG")
}

## Supppose that we know the true CPDAG and covariance matrix
(l.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "local", type = "cpdag"))
(o.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "optimal", type = "cpdag"))
\dontrun{(g.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "global", type = "cpdag"))}
## All three methods produce the same unique values. 

## Supppose that we know the true PDAG and covariance matrix
(l.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "local", type = "pdag"))
(o.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "optimal", type = "pdag"))
\dontrun{(g.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "global", type = "pdag"))}
## All three methods produce the same unique values.

## From the true DAG, we can compute the true causal effect of 3 on 10
(ce.3.10 <- causalEffect(myDAG, 10, 3))
## Indeed, this value is contained in the values found by all methods

## When working with data we have to use the estimated CPDAG and
## the sample covariance matrix
(l.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph, method = "local", type = "cpdag"))
(o.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph, method = "optimal", type = "cpdag"))
\dontrun{(g.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph,
method = "global", type = "cpdag"))}
## The unique values of the local and the global method are still identical.
## While not identical, the values of the optimal method are very similar.
## The true causal effect is contained in all three sets, up to a small
## estimation error (0.118 vs. 0.112 with true value 0.114) 

## Similarly, when working with data and background knowledge we have to use the estimated PDAG and
## the sample covariance matrix
(l.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "local", type = "pdag"))
(o.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "optimal", type = "pdag"))
\dontrun{(g.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "global", type = "pdag"))}
## The unique values of the local and the global method are still identical.
## While not necessarily identical, the values of the optimal method will be similar.

## The true causal effect is contained in both sets, up to a small estimation error

## All three can also be applied to sets y.
(l.ida.cpdag.2 <- ida(3,c(6,10), cov(dat), pc.fit@graph, method = "local", type = "cpdag"))
(o.ida.cpdag.2 <- ida(3,c(6,10), cov(dat), pc.fit@graph, method = "optimal", type = "cpdag"))
\dontrun{(g.ida.cpdag.2 <- ida(3,c(6,10), cov(dat), pc.fit@graph,
method = "global", type = "cpdag"))}
## For the methods local and global we recommend use of idaFast in this case for better performance.

## Note that only the optimal method can be appplied to sets x.
(o.ida.cpdag.2 <- ida(c(2,3),10, cov(dat), pc.fit@graph, method = "optimal", type = "cpdag"))

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

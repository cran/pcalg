\name{ida}
\title{Estimate Multiset of Possible Total Causal Effects}
\alias{ida}
\alias{causalEffect}
\description{
  \code{ida()} estimates the multiset of possible total causal effects
  of one variable (\code{x}) onto another variable (\code{y})
  from observational data.

  \code{causalEffect(g,y,x)} computes the true causal effect
  (\eqn{\beta_{\mathrm{true}}}{beta.true}) of x on y in g.
}
\usage{
ida(x.pos, y.pos, mcov, graphEst, method = c("local","global"),
    y.notparent = FALSE, verbose = FALSE, all.dags = NA)

causalEffect(g, y, x)
}
\arguments{
  \item{x.pos, x}{(integer) position of variable \code{x} in the covariance matrix.}
  \item{y.pos, y}{(integer) position of variable \code{y} in the covariance matrix.}
  \item{mcov}{Covariance matrix that was used to estimate \code{graphEst}.}
  \item{graphEst}{Estimated CPDAG, typically from \code{\link{pc}()}: If
    the result of \code{\link{pc}} is \code{pc.fit}, the estimated CPDAG
    can be obtained by \code{pc.fit@graph}.}
  \item{method}{Character string specifying the method with default \code{"local"}.
    \describe{
      \item{\code{"global"}:}{The algorithm considers all DAGs in the
      equivalence class of the CPDAG, hence is \emph{slow}.}
      \item{\code{"local"}:}{The algorithm only considers the local
      neighborhood of \code{x} in the CPDAG, hence is \emph{faster}.}
    }
    See details below.}
  \item{y.notparent}{Logical; if true, any edge between \code{x} and
    \code{y} is assumed to be of the form \code{x->y}.}
  \item{verbose}{If TRUE, details on the regressions are printed.}
  \item{all.dags}{All DAGs in the equivalence class of the
    CPDAG can be precomputed by \code{\link{allDags}()} and passed via
    this argument.  In that case, \code{\link{allDags}(..)} is not called
    internally.
    This option is only relevant when using \code{method="global"}.}
  \item{g}{Graph in which the causal effect is sought.}
}
\details{
  It is assumed that we have observational data that are multivariate
  Gaussian and faithful to the true (but unknown) underlying causal DAG
  (without hidden variables).  Under these assumptions, this function
  estimates the multiset of possible total causal effects of \code{x} on
  \code{y}, where the total causal effect is defined via Pearl's
  do-calculus as \eqn{E(Y|do(X=z+1))-E(Y|do(X=z))} (this value does not
  depend on \eqn{z}, since Gaussianity implies that conditional
  expectations are linear).

  We estimate a \emph{set} of possible total causal effects instead of
  the unique total causal effect, since it is typically impossible to
  identify the latter when the true underlying causal DAG is unknown
  (even with an infinite amount of data).  Conceptually, the method
  works as follows.  First, we estimate the equivalence class of DAGs
  that describe the conditional independence relationships in the data,
  using the function \code{\link{pc}} (see the help file of this
  function).  For each DAG G in the equivalence class, we apply Pearl's
  do-calculus to estimate the total causal effect of \code{x} on
  \code{y}.  This can be done via a simple linear regression: if \code{y}
  is not a parent of \code{x}, we take the regression coefficient of
  \code{x} in the regression \code{lm(y ~ x + pa(x))}, where
  \code{pa(x)} denotes the parents of \code{x} in the DAG G; if \code{y}
  is a parent of \code{x} in G, we set the estimated causal effect to
  zero.

  If the equivalence class contains \code{k} DAGs, this will yield
  \code{k} estimated total causal effects.  Since we do not know which DAG
  is the true causal DAG, we do not know which estimated total causal
  effect of \code{x} on \code{y} is the correct one.  Therefore, we return
  the entire multiset of \code{k} estimated effects (it is a multiset
  rather than a set because it can contain duplicate values).

  One can take summary measures of the multiset.  For example, the
  minimum absolute value provides a lower bound on the size of the true
  causal effect: If the minimum absolute value of all values in the
  multiset is larger than one, then we know that the size of the true
  causal effect (up to sampling error) must be larger than one.

  If \code{method="global"}, the method as described above is carried
  out, where all DAGs in the equivalene class of the estimated CPDAG
  \code{graphEst} are computed using the function \code{\link{allDags}}.
  This method is suitable for small graphs (say, up to 10 nodes).

  If \code{method="local"}, we do not determine all DAGs in the
  equivalence class of the CPDAG.  Instead, we only consider the local
  neighborhood of \code{x} in the CPDAG.  In particular, we consider all
  possible directions of undirected edges that have \code{x} as
  endpoint, such that no new v-structure is created.  For each such
  configuration, we estimate the total causal effect of \code{x} on
  \code{y} as above, using linear regression.

  At first sight, it is not clear that such a local configuration
  corresponds to a DAG in the equivalence class of the CPDAG, since it
  may be impossible to direct the remaining undirected edges without
  creating a directed cycle or a v-structure.  However, Maathuis,
  Kalisch and Buehlmann (2009) showed that there is at least one DAG in
  the equivalence class for each such local configuration.
  %%
  As a result, it follows that the multisets of total causal effects of
  the "global" and the "local" method have the same unique values. They
  may, however, have different multiplicities.

  For example, a CPDAG may represent eight DAGs, and the global method
  may produce the multiset \{1.3, -0.5, 0.7, 1.3, 1.3, -0.5, 0.7, 0.7\}.
  The unique values in this set are -0.5, 0.7 and 1.3, and the
  multiplicities are 2, 3 and 3.  The local method, on the other hand,
  may yield \{1.3, -0.5, -0.5, 0.7\}.  The unique values are again -0.5,
  0.7 and 1.3, but the multiplicities are now 2, 1 and 1.  The fact that
  the unique values of the multisets of the "global" and "local" method
  are identical implies that summary measures of the multiset that only
  depend on the unique values (such as the minimum absolute value) can
  be estimate by the local method.

}
\value{
  A vector that represents the multiset containing the estimated
  possible total causal effects of \code{x} on \code{y}.
}
\references{
  M.H. Maathuis, M. Kalisch, P. Buehlmann (2009).
  Estimating high-dimensional intervention effects from observational data.
  \emph{Annals of Statistics} \bold{37}, 3133--3164.

  M.H. Maathuis, D. Colombo, M. Kalisch, P. Bühlmann (2010).
  Predicting causal effects in large-scale systems from observational data.
  \emph{Nature Methods} \bold{7}, 247--248.

  Markus Kalisch, Martin Maechler, Diego Colombo, Marloes H. Maathuis,
  Peter Buehlmann (2012).
  Causal Inference Using Graphical Models with the R Package pcalg.
  \emph{Journal of Statistical Software} \bold{47}(11) 1--26,
  \url{http://www.jstatsoft.org/v47/i11/}.

  Pearl (2005). \emph{Causality. Models, reasoning and inference}.
  Cambridge University Press, New York.
}
\author{Markus Kalisch (\email{kalisch@stat.math.ethz.ch})}
\seealso{
  \code{\link{idaFast}} for estimating the multiset of possible total
  causal effects for several target variables simultaneously.

  \code{\link{pc}} for estimating a CPDAG.
}
\examples{
## Simulate the true DAG
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

## simulate Gaussian data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)

## estimate CPDAG -- see  help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01)

if (require(Rgraphviz)) {
  ## plot the true and estimated graphs
  par(mfrow = c(1,2))
  plot(myDAG, main = "True DAG")
  plot(pc.fit, main = "Estimated CPDAG")
}

## Supppose that we know the true CPDAG and covariance matrix
(l.ida <- ida(2,5, covTrue, myCPDAG, method = "local"))
(g.ida <- ida(2,5, covTrue, myCPDAG, method = "global"))
## The global and local method produce the same unique values.
stopifnot(all.equal(sort(unique(g.ida)),
                    sort(unique(l.ida))))

## From the true DAG, we can compute the true causal effect of 2 on 5
(ce.2.5 <- causalEffect(myDAG, 5, 2))
## Indeed, this value is contained in the values found by both the
## local and global method

## When working with data, we have to use the estimated CPDAG and
## the sample covariance matrix
(l.ida <- ida(2,5, cov(dat), pc.fit@graph, method = "local"))
(g.ida <- ida(2,5, cov(dat), pc.fit@graph, method = "global"))
## The unique values of the local and the global method are still identical,
stopifnot(all.equal(sort(unique(g.ida)), sort(unique(l.ida))))
## and the true causal effect is contained in both sets, up to a small
## estimation error (0.868 vs. 0.857)
stopifnot(all.equal(ce.2.5, max(l.ida), tol = 0.02))
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

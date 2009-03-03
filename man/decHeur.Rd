\name{decHeur}
\alias{decHeur}
\title{Decision Heuristic: Should robust method for PC-algorithm be used?}
\description{
Simple Heuristic for deciding whether robust PC-algorithm should be used.
}
\usage{decHeur(dat,gam=0.05,sim.method="t",est.method="o",n.sim=100,two.sided=FALSE,verbose=FALSE)
}
\arguments{
  \item{dat}{Data matrix (cols=variables, rows=samples)}
  \item{gam}{Significance level for test}
  \item{sim.method}{Reference distribution; "n" for Normal, "t" for
    N+10\% t3}
  \item{est.method}{Estimation method of correlation matrix; "s" for
    standard, "o" for OGK using Qn (robust)}
  \item{n.sim}{Number of samples drawn from reference distribution}
  \item{two.sided}{Should a two-sided test be used?}
  \item{verbose}{Run in verbose mode}
}
\value{
  \item{tvec}{Simulated values of test statistic}
  \item{tval}{Observed value of test statistic}
  \item{outlier}{Is robust method suggested? (TRUE=Suggested)}
}
\details{
Simulation studies show that the standard PC-algorithm already is
rather insensitive to outliers, provided, they are not too severe. The
effect of very heavy outliers can be dramatically reduced by using the
robust PC-algorithm; this increases the computational burden by roughly one
order of magnitude.

We provide a simple method for deciding whether data at hand has worse
outliers than a given reference distribution. Using this, we see two
heuristics for deciding whether to use the robust version of the
PC-algorithm or not. On the one hand, one could use the normal distribution
as reference distribution and apply the robust PC-algorithm to all data
that seem to contain more outliers than an appropriate normal distribution
(Heuristic A). On the other hand, one could, inspired by the results of 
simulation studies, only want to apply the robust
method in the case where the contamination is worse than a normal
distribution with 10\% outliers from a $t_3$ distribution. Then, we would
use a normal distribution with 10\% outliers from a $t_3$ distribution as
reference distribution (Heuristic B).

In order to decide whether data has worse outliers than a given reference
distribution, we proceed as follows.  We compute a robust estimate of the
covariance matrix of the data (e.g. OGK with Qn-estimator) and simulate
(several times) data from the reference distribution with this covariance
matrix. For each dimension \eqn{i \ (1 \leq i \leq p)}, we compute the ratio of
standard deviation \eqn{\sigma_i} and a robust version of it $s_i$ (e.g.,
Qn-estimator) and compute the average over all dimensions. (Since the main
input for the PC-algorithm are correlation estimates which can be expressed
in terms of scale estimates, we base our test
statistics on scale estimates.) Thus, we obtain the distribution of this
averaged ratio \eqn{R= \frac{1}{p} \sum_{i=1}^{p} \sigma_i / s_i} under the
null hypothesis that the data can be explained by the reference
distribution with given covariance matrix. We now can test this null
hypothesis by using the ratio computed with the current data set
\eqn{r=\frac{1}{p} \sum_{i=1}^{p} \hat{\sigma_i}/\hat{s_i}} on a given
significance level.
}

\seealso{\code{\link{pcAlgo}} which can be used with standard and robust
  correlation estimate.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch})
}
\examples{
set.seed(123)
## generate a data set for the example
p <- 5
myDAG <- randomDAG(p, prob = 0.6)
n <- 1000
## data without outlier
datN <- rmvDAG(n, myDAG, errDist = "normal")
## data with severe outlier (10\% Cauchy)
datC <- rmvDAG(n, myDAG, errDist = "mix")
n.sim <- 20
gam <- 0.05
sim.method <- "t"
est.method <- "o"
decHeur(datN,gam,sim.method,est.method,n.sim=n.sim,two.sided=FALSE,verbose=TRUE)
decHeur(datC,gam,sim.method,est.method,n.sim=n.sim,two.sided=FALSE,verbose=TRUE)
}
\keyword{multivariate}
\keyword{models}
\keyword{robust}
\keyword{graphs}

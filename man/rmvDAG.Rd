\name{rmvDAG}
\alias{rmvDAG}
\title{Generate Multivariate Data according to a DAG}
\description{
  Generate multivariate data with dependency structure specified by a
  (given) DAG (\bold{D}irected \bold{A}cyclic \bold{G}raph) with nodes
  corresponding to random variables. The DAG has to be
  \bold{topologically ordered}.
}
\usage{
rmvDAG(n, dag,
             errDist = c("normal", "cauchy", "mix", "mixt3", "mixN100",
"t4"), mix = 0.1, errMat = NULL)
}
\arguments{
  \item{dag}{a graph object describing the DAG; must contain weights for
    all the edges. The nodes must be topologically sorted. (For
    topological sorting use \code{\link[RBGL]{tsort}} from the
    \pkg{RBGL} package.)}
  \item{n}{number of samples that should be drawn. (integer)}
  \item{errDist}{Specifies the distribution of each node. Currently, the
    options "normal", "t4", "cauchy", "mix", "mixt3" and "mixN100" are
    supported. The first 
    three generate standard normal-, t(df=4)- and cauchy-random
    numbers. The options containing the word "mix" create standard
    normal random variables with a mix of outliers. The outliers for the
    options "mix", "mixt3", "mixN100" are drawn from a standard cauchy,
    t(df=3) and N(0,100) distribution, respectively. The fraction of 
    cauchy points can be adjusted using the parameter "mix". (string)}
  \item{mix}{For \code{errDist="mix"} this parameter specifies the
    percentage of samples that are drawn from a standard cauchy
    distribution. (numeric)}
  \item{errMat}{numeric \eqn{n * p} matrix specifiying the error vectors
    \eqn{e_i} (see Details), instead of specifying \code{errDist} (and
    maybe \code{mix}).}
}
\value{
  A \eqn{n*p} matrix with the generated data. The \eqn{p} columns
  correspond to the nodes (i.e., random variables) and each of the
  \eqn{n} rows correspond to a sample.
}
\details{
  Each node is visited in the topological order. For each node \eqn{i} we
  generate a \eqn{p}-dimensional value \eqn{X_i} in the following way:
  Let \eqn{X_1,\ldots,X_k} denote the values of all neighbours of \eqn{i}
  with lower order.
  Let \eqn{w_1,\ldots,w_k} be the weights of the corresponding edges.
  Furthermore, generate a random vector \eqn{e_i} according to the
  specified error distribution.  Then, the value of \eqn{X_i} is
  computed as
  \deqn{X_i = w_1*X_1 + \ldots + w_k*X_k + e_i.}
  If node \eqn{i} has no neighbors with lower order, \eqn{X_i} is computed as
  \eqn{X_i = e_i}.
}
\seealso{\code{\link{randomDAG}} for generating a random DAG;
  \code{\link{skeleton}} and \code{\link{pc}} for estimating the
  skeleton and the CPDAG of a DAG that
  corresponds to the data.
}
\author{
  Markus Kalisch (\email{kalisch@stat.math.ethz.ch}) and Martin Maechler.
}
\examples{
## generate random DAG
p <- 20
rDAG <- randomDAG(p, prob = 0.2, lB=0.1, uB=1)

if (require(Rgraphviz)) {
## plot the DAG
plot(rDAG, main = "randomDAG(20, prob = 0.2, ..)")
}

## generate 1000 samples of DAG using standard normal error distribution
n <- 1000
d.normMat <- rmvDAG(n, rDAG, errDist="normal")

## generate 1000 samples of DAG using standard t(df=4) error distribution
d.t4Mat <- rmvDAG(n, rDAG, errDist="t4")

## generate 1000 samples of DAG using standard normal with a cauchy
## mixture of 30 percent
d.mixMat <- rmvDAG(n, rDAG, errDist="mix",mix=0.3)

require(MASS) ## for mvrnorm()
Sigma <- toeplitz(ARMAacf(0.2, lag.max = p - 1))
dim(Sigma)# p x p
## *Correlated* normal error matrix "e_i" (against model assumption)
eMat <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
d.CnormMat <- rmvDAG(n, rDAG, errMat = eMat)
}
\keyword{multivariate}
\keyword{datagen}

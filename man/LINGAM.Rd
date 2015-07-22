\name{LINGAM}
\title{Linear non-Gaussian Additive Models (LiNGAM)}
\alias{LINGAM}
\description{
  Fits a Linear non-Gaussian Additive Model (LiNGAM) to the data and
  returns the corresponding DAG.

  For details, see the reference below.
}
\usage{
LINGAM(X, verbose = FALSE)
}
\arguments{
  \item{X}{n x p data matrix (n: sample size, p: number of variables).}
  \item{verbose}{logical or integer indicating that increased diagnostic output is to
    be provided.}
}
%% \details{
%% }
\value{
  list with components
  \item{Adj}{a \eqn{p \times p}{p x p} 0/1 adjacency matrix \eqn{A}.
    \code{A[i,j] == 1} corresponds to a directed edge from i to j.}
  \item{B}{\eqn{p \times p}{p x p} matrix of corresponding linear
    coefficients.  Note it corresponds to the \emph{transpose} of
    \code{Adj}, i.e., \code{identical( Adj, t(B) != 0 )} is true.}
}
\references{
  S. Shimizu, P.O. Hoyer, A. Hyv\"arinen, A. Kerminen (2006)
  A Linear Non-Gaussian Acyclic Model for Causal Discovery;
  \emph{Journal of Machine Learning Research} \bold{7}, 2003--2030.
}
\author{
  Patrik Hoyer <patrik.hoyer@helsinki.fi>, Doris Entner
  <entnerd@hotmail.com>, Antti Hyttinen <antti.hyttinen@cs.helsinki.fi>
  and Jonas Peters <jonas.peters@tuebingen.mpg.de>
}
\seealso{
  \code{\link{fastICA}} from package \pkg{fastICA} is used.
}
\examples{
##################################################
## Exp 1
##################################################
set.seed(123)
n <- 500
eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
eps2 <- runif(n) - 0.5

x2 <- eps2
x1 <- 0.9*x2 + eps1

X <- cbind(x1,x2)

trueDAG <- cbind(c(0,1),c(0,0))
## x1 <- x2
## adjacency matrix:
## 0 0
## 1 0

estDAG <- LINGAM(X)

cat("true DAG:\n")
show(trueDAG)

cat("estimated DAG:\n")
show(estDAG$Adj)

##################################################
## Exp 2
##################################################
set.seed(123)
n <- 500
eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
eps2 <- runif(n) - 0.5
eps3 <- sign(rnorm(n)) * abs(rnorm(n))^(1/3)
eps4 <- rnorm(n)^2

x2 <-                eps2
x1 <-   0.9*x2     + eps1
x3 <-   0.8*x2     + eps3
x4 <- -x1  -0.9*x3 + eps4

X <- cbind(x1,x2,x3,x4)


trueDAG <- cbind(x1 = c(0,1,0,0),
                 x2 = c(0,0,0,0),
                 x3 = c(0,1,0,0),
                 x4 = c(1,0,1,0))
## x4 <- x3 <- x2 -> x1 -> x4
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0

estDAG1 <- LINGAM(X, verbose = TRUE)# details on LINGAM
estDAG2 <- LINGAM(X, verbose = 2)   # details on LINGAM and fastICA
## results are the same, of course:
stopifnot(identical(estDAG1, estDAG2))
cat("true DAG:\n")
show(trueDAG)

cat("estimated DAG:\n")
show(estDAG1$Adj)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

\name{LINGAM}
\alias{LINGAM}
\title{
Linear non-Gaussian Additive Models
}
\description{
This code fits a LiNGAM to the data and outputs the corresponding DAG.
}
\usage{
LINGAM(X, verbose = FALSE)
}
\arguments{
  \item{X}{
n x p data matrix (n: sample size, p: number of variables)
}
  \item{verbose}{
boolean indicator (TRUE/FALSE) whether output should be shown.
}
}
\details{
see the reference below
}
\value{
list of objects:
\item{Adj}{p x p adjacency matrix. The entry (i,j) being one corresponds to a directed edge from i to j.}
\item{B}{p x p matrix of corresponding linear coefficients.}
}
\references{
S. Shimizu, P.O. Hoyer, A. Hyv\"arinen, A. Kerminen: A Linear Non-Gaussian Acyclic Model for Causal Discovery, Journal of Machine Learning Research 7, pages 2003-2030, 2006.
}
\author{
Patrik Hoyer <patrik.hoyer@helsinki.fi>, Doris Entner <entnerd@hotmail.com>, Antti Hyttinen <antti.hyttinen@cs.helsinki.fi> and Jonas Peters <jonas.peters@tuebingen.mpg.de>
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

estDAG <- LINGAM(X, verbose = TRUE)

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

x2 <- eps2
x1 <- 0.9*x2 + eps1
x3 <- 0.8*x2 + eps3
x4 <- -0.9*x3 - x1 + eps4

X <- cbind(x1,x2,x3,x4)

trueDAG <- cbind(c(0,1,0,0),c(0,0,0,0),c(0,1,0,0),c(1,0,1,0))
## x4 <- x3 <- x2 -> x1 -> x4
## adjacency matrix:
## 0 0 0 1
## 1 0 1 0
## 0 0 0 1
## 0 0 0 0

estDAG <- LINGAM(X, verbose = TRUE)

cat("true DAG:\n")
show(trueDAG)

cat("estimated DAG:\n")
show(estDAG$Adj)
}
\keyword{Causality }
\keyword{Linear Non-Gaussian }
\keyword{Restricted Structural Equation Model}


library(pcalg)
## library(fastICA)
## library(clue)
## source("/u/kalischm/research/packages/LINGAM/R/lingamFuns.R")

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

stopifnot(all(as.numeric(estDAG$Adj) == trueDAG))

## cat("true DAG:\n")
## show(trueDAG)

## cat("estimated DAG:\n")
## show(estDAG$Adj)

## using pcalg
## n <- nrow(X)
## V <- as.character(1:ncol(X)) # labels aka node names
     
## ## estimate CPDAG
## pc.fit <- pc(suffStat = list(C = cor(X), n = n),
##              indepTest = gaussCItest, ## indep.test: partial correlations
##              alpha=0.01, labels = V, verbose = FALSE)
## if (require(Rgraphviz)) {
##     plot(pc.fit, main = "Estimated CPDAG")
## }

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

stopifnot(all(as.numeric(estDAG$Adj) == trueDAG))

## cat("true DAG:\n")
## show(trueDAG)

## cat("estimated DAG:\n")
## show(estDAG$Adj)

## using pcalg
## n <- nrow(X)
## V <- as.character(1:4) # labels aka node names
     
## ## estimate CPDAG
## pc.fit <- pc(suffStat = list(C = cor(X), n = n),
##              indepTest = gaussCItest, ## indep.test: partial correlations
##              alpha=0.01, labels = V, verbose = FALSE)
## if (require(Rgraphviz)) {
##     plot(pc.fit, main = "Estimated CPDAG")
## }


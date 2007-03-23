library(pcalg)
## library(pcalg,lib.loc="/u/kalisch/R/Pkgs/pcalg.Rcheck")

## This function tests randomDAG

## make results reproducable
set.seed(47)
mlB <- 0.1
muB <- 1
if(class(try(randomDAG(20,0.2, lB=muB, uB=mlB), silent=TRUE)) != "try-error")
    stop("Testing problem in test_randomDAG.R") ## should return error message
rd <- randomDAG(20,0.2,lB=mlB,uB=muB)

##- rd <- randomDAG(10,0.8)
##- plot(rd)

## The total number of edges should be Bin(p(p-1)/2;s);
## check this for some s with Chi^2-Test
set.seed(42)
p <- 30
s <- 0.2
n.rep <- 100
resO <- numeric(n.rep)
for (i in 1:n.rep)
  resO[i] <- numEdges(ugraph(randomDAG(p,s)))  ## number of edges
resE <- rbinom(n.rep, size = p*(p-1)/2, prob = s)
rng <- range(resO, resE)
mybreaks <- seq(from = rng[1], to = rng[2], length = 5+1)
histO <- tabulate(cut(resO, mybreaks, include.lowest = TRUE, labels = FALSE), 5)
histE <- tabulate(cut(resE, mybreaks, include.lowest = TRUE, labels = FALSE), 5)
(chi2res <- chisq.test(histO,histE))
stopifnot(all.equal(chi2res$p.value, 0.24143,tolerance=1e-4))


## Check whether really acyclic - will produce error on interface
##- library(RBGL)
##- n.rep <- 100
##- for (i in 1:n.rep) tsort(randomDAG(p,s))

## Check distribution of weights
library(pcalg)
p <- 10
s <- 0.2
mlB <- 0.1
muB <- 1

set.seed(42)
g <- randomDAG(p,s)
resWO <- wgtMatrix(g) ## weights
resWO <- resWO[resWO > 0]
resWE <- runif(length(resWO), min=mlB, max=muB)
(resWE)
brks <- seq(from= mlB, to=muB, length = 6+1)
(brks)
nWO <- hist(resWO, breaks=brks)$counts
nWE <- hist(resWE, breaks=brks)$counts
(chi2res <- chisq.test(nWO,nWE))

stopifnot(all.equal(chi2res$p.value, 0.55782,tolerance=1e-5))

## Test function call for several parameter settings
##- pValues <- seq(3,100,by=1)
##- sValues <- seq(0,1,by=0.1)
##- for (p in pValues) {
##-   for (s in sValues) {
##-     randomDAG(p,s)
##-   }
##- }

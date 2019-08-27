library(pcalg)

suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
nreps <- 100
res.local <- logical(nreps)
res.opt <- logical(nreps)
all.eff.true <- res.local
Rnd <- function(e) round(e, 14)## get 14 digits accuracy, as we use true (DAG, cov)
Rnd7 <- function(e) round(e, 7)## get 14 digits accuracy, as we use true (DAG, cov)
for (i in 1:nreps) {
  p <- 2 + rpois(1, lambda = 8) # ==>  p >= 2, E[p] = 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.2)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)
  
  ## x != y  in {1,2,...p} ;
  xy <- sample.int(p, 2); x <- xy[1]; y <- xy[2]
  
  ## plot(myCPDAG)
  eff.true <- Rnd(causalEffect(myDAG, y, x))
  all.eff.true[i] <- eff.true
  ## cat("x=",x," y=",y," eff=",eff.true,"\n")
  
  eff.est.local <- Rnd(ida(x,y, mcov, myCPDAG, method="local"))
  eff.est.opt <- Rnd(ida(x,y, mcov, myCPDAG, method="optimal"))
  res.local[i] <- (eff.true %in% eff.est.local)
  res.opt[i] <- (eff.true %in% eff.est.opt)
}
## cat('Time elapsed: ', (.pt <- proc.time()),"\n")

## stem(all.eff.true)
if (!all(res.local)) stop("Test ida: True effects were not recovered by local method!")
if (!all(res.opt)) stop("Test ida: True effects were not recovered by optimal method!")

## *one* test for  method="global" :
eff.g.est <- Rnd(ida(x,y, mcov, myCPDAG, method="global", verbose=TRUE))
stopifnot(eff.est.local == eff.g.est)

## cat('Time elapsed additionally: ', proc.time() - .pt,"\n")

## another special case (from Raphael Gervais)
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
amatT <- as(myDAG, "matrix") # weighted adjacency matrix of true DAG
effT <- Rnd(amatT[2,3]*amatT[3,5]) # Causal effect of 2 on 5 from true DAG weighted matrix
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix
effG <- Rnd(ida(2,5, covTrue,myCPDAG,method = "global"))

if (!(effT %in% effG)) stop("Test ida special case: True effects were not recovered!")

##################################################
## Tests for method = optimal; sets x
##################################################
set.seed(1)
V <- sample(c("C1","C2","C3","X","C4","C5","C6","C7","Y","C8","C9","C10"))
myDAG <- randomDAG(20, 0.3)
mcov <- trueCov(myDAG)
amat <- t(as(myDAG,"matrix"))
graphEst <- dag2cpdag(myDAG)
amat.cpdag <- t(as(graphEst,"matrix"))

opt <- ida(c(4,6),12,trueCov(myDAG),graphEst,method="optimal",type="cpdag") 
RRC <- jointIda(c(4,6),12,trueCov(myDAG),graphEst,technique="RRC")
stopifnot(all.equal(opt,RRC, tolerance = 0.01))

opt <- ida(c(4,6),c(10,19),trueCov(myDAG),graphEst,method="optimal",type="cpdag")
RRC <- jointIda(c(4,6),c(10,19),trueCov(myDAG),graphEst,technique="MCD")
stopifnot(all.equal(opt,RRC, tolerance = 0.01))

opt <- ida(c(5,10),c(3),trueCov(myDAG),graphEst,method="optimal",type="cpdag")
RRC <- jointIda(c(5,10),c(3),trueCov(myDAG),graphEst,technique="RRC")
stopifnot(all.equal(opt,RRC, tolerance = 0.01))

## sometimes they differ
## ida(c(5,20),c(10),trueCov(myDAG),graphEst,method="optimal",type="cpdag")
## jointIda(c(5,20),c(10),trueCov(myDAG),graphEst,technique="RRC")

##################################################
## Tests related to examples (use dontrun for slow global option there)
##################################################
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

## Supppose that we know the true CPDAG and covariance matrix
(l.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "local", type = "cpdag"))
(o.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "optimal", type = "cpdag"))
(g.ida.cpdag <- ida(3,10, covTrue, myCPDAG, method = "global", type = "cpdag"))
## All three methods produce the same unique values.
stopifnot(all.equal(sort(unique(Rnd7(g.ida.cpdag))),
                    sort(unique(Rnd7(l.ida.cpdag)))))
stopifnot(all.equal(sort(unique(Rnd7(g.ida.cpdag))),
                    sort(unique(Rnd7(as.vector(o.ida.cpdag))))))

## Supppose that we know the true PDAG and covariance matrix
(l.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "local", type = "pdag"))
(o.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "optimal", type = "pdag"))
(g.ida.pdag <- ida(3,10, covTrue, myPDAG, method = "global", type = "pdag"))
## All three methods produce the same unique values.
stopifnot(all.equal(sort(unique(Rnd7(g.ida.pdag))),
                    sort(unique(Rnd7(l.ida.pdag)))))
stopifnot(all.equal(sort(unique(Rnd7(g.ida.pdag))),
                    sort(unique(Rnd7(as.vector(o.ida.pdag))))))

## From the true DAG, we can compute the true causal effect of 3 on 10
(ce.3.10 <- causalEffect(myDAG, 10, 3))
## Indeed, this value is contained in the values found by all methods

## When working with data we have to use the estimated CPDAG and
## the sample covariance matrix
(l.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph, method = "local", type = "cpdag"))
(o.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph, method = "optimal", type = "cpdag"))
(g.ida.est.cpdag <- ida(3,10, cov(dat), pc.fit@graph, method = "global", type = "cpdag"))
## The unique values of the local and the global method are still identical.
stopifnot(all.equal(sort(unique(Rnd7(g.ida.est.cpdag))), sort(unique(Rnd7(l.ida.est.cpdag)))))
## While not identical, the values of the optimal method are very similar.
stopifnot(all.equal(sort(o.ida.est.cpdag), sort(l.ida.est.cpdag), tolerance = 0.025))
## The true causal effect is contained in all three sets, up to a small
## estimation error (0.118 vs. 0.112 with true value 0.114) 
stopifnot(all.equal(ce.3.10, min(l.ida.est.cpdag), tolerance = 0.04))
stopifnot(all.equal(ce.3.10, min(o.ida.est.cpdag), tolerance = 0.02))

## Similarly, when working with data and background knowledge we have to use the estimated PDAG and
## the sample covariance matrix
(l.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "local", type = "pdag"))
(o.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "optimal", type = "pdag"))
(g.ida.est.pdag <- ida(3,10, cov(dat), pc.fit.pdag, method = "global", type = "pdag"))
## The unique values of the local and the global method are still identical.
stopifnot(all.equal(sort(unique(Rnd7(g.ida.est.pdag))), sort(unique(Rnd7(l.ida.est.pdag)))))
## While not necessarily identical, the values of the optimal method will be similar.
stopifnot(all.equal(sort(Rnd7(o.ida.est.pdag)), sort(Rnd7(l.ida.est.pdag)), tolerance = 0.08))
## The true causal effect is contained in both sets, up to a small estimation error
stopifnot(all.equal(ce.3.10, min(l.ida.est.pdag), tolerance = 0.04))
stopifnot(all.equal(ce.3.10, min(o.ida.est.pdag), tolerance = 0.02))


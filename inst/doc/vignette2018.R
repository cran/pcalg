### R code from vignette source 'vignette2018.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
op.orig <-
options(width = 79,
        ## SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        useFancyQuotes = FALSE,
        ## prompt="R> ", continue="+  " # << JSS
           prompt="> " , continue="  "  # << nicer
        )
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")

PCversion <- packageDescription('pcalg')$Version # *not* "2.6-0" - needs updating


###################################################
### code chunk number 2: loadPackage
###################################################
library(pcalg)
stopifnot(require(Rgraphviz))# needed for all our graph plots
library(graph)


###################################################
### code chunk number 3: show-args
###################################################
showF <- function(f, width = 50) {
    ## 'width': larger than default on purpose:
    nam <- deparse(substitute(f))
    stopifnot(is.function(f))
    ll <- capture.output(str(removeSource(f), width=width))
    ll[1] <- sub("function *", nam, ll[1])
    writeLines(ll)
}


###################################################
### code chunk number 4: exIntro1
###################################################
library("pcalg")
data("gmG") ## loads data set gmG8


###################################################
### code chunk number 5: exIntroSkel
###################################################
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
varNames <- gmG8$g@nodes
skel.gmG8 <- skeleton(suffStat, indepTest = gaussCItest,
                      labels = varNames, alpha = 0.01)


###################################################
### code chunk number 6: exIntroPC
###################################################
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
varNames <- gmG8$g@nodes
pc.gmG8 <- pc(suffStat, indepTest = gaussCItest,
              labels = varNames, alpha = 0.01)


###################################################
### code chunk number 7: exIntroPlot
###################################################
stopifnot(require(Rgraphviz))# needed for all our graph plots
par(mfrow = c(1,3))
plot(gmG8$g, main = "");    box(col="gray")
plot(skel.gmG8, main = ""); box(col="gray")
plot(pc.gmG8, main = "");   box(col="gray")


###################################################
### code chunk number 8: exIntroFCI
###################################################
data("gmL")
suffStat <- list(C = cor(gmL$x), n = nrow(gmL$x))
fci.gmL <- fci(suffStat, indepTest=gaussCItest,
               alpha = 0.9999, labels = c("2","3","4","5"))


###################################################
### code chunk number 9: exIntroPlotFCI
###################################################
stopifnot(require(Rgraphviz))# needed for all our graph plots
par(mfrow = c(1,2))
plot(gmL$g) ; box(col="gray")
plot(fci.gmL); box(col="gray")


###################################################
### code chunk number 10: exIntroFCIplus
###################################################
suffStat <- list(C = cor(gmL$x), n = nrow(gmL$x))
fciPlus.gmL <- fciPlus(suffStat, indepTest=gaussCItest,
                       alpha = 0.9999, labels = c("2","3","4","5"))


###################################################
### code chunk number 11: obs-score-args (eval = FALSE)
###################################################
## score <- new("GaussL0penObsScore", data = matrix(1, 1, 1),
##     lambda = 0.5*log(nrow(data)), intercept = FALSE,
##     use.cpp = TRUE, ...)


###################################################
### code chunk number 12: ges-args
###################################################
showF(ges)


###################################################
### code chunk number 13: gesExpl-plot
###################################################
score <- new("GaussL0penObsScore", gmG8$x)
ges.fit <- ges(score)
par(mfrow=1:2)
plot(gmG8$g, main = "")          ; box(col="gray")
plot(ges.fit$essgraph, main = ""); box(col="gray")


###################################################
### code chunk number 14: int-score-args (eval = FALSE)
###################################################
## score <- new("GaussL0penIntScore", data = matrix(1, 1, 1),
##     targets = list(integer(0)),
##     target.index = rep(as.integer(1), nrow(data)),
##     lambda = 0.5*log(nrow(data)), intercept = FALSE,
##     use.cpp = TRUE, ...)


###################################################
### code chunk number 15: load-gmInt
###################################################
data(gmInt)


###################################################
### code chunk number 16: inspect-gmInt
###################################################
n.tot <- length(gmInt$target.index)
n.obs <- sum(gmInt$target.index == 1)
n3 <- sum(gmInt$target.index == 2)
n5 <- sum(gmInt$target.index == 3)


###################################################
### code chunk number 17: gies-args
###################################################
showF(gies, width = 60)


###################################################
### code chunk number 18: def-gmInt (eval = FALSE)
###################################################
## ## Used to generate the  'gmInt'  Gaussian data originally:
## set.seed(40)
## p <- 8
## n <- 5000
## gGtrue <- randomDAG(p, prob = 0.3)
## pardag <- as(gGtrue, "GaussParDAG")
## pardag$set.err.var(rep(1, p))
## targets <- list(integer(0), 3, 5)
## target.index <- c(rep(1, 0.6*n), rep(2, n/5), rep(3, n/5))
## 
## x1 <- rmvnorm.ivent(0.6*n, pardag)
## x2 <- rmvnorm.ivent(n/5, pardag, targets[[2]],
##                     matrix(rnorm(n/5, mean = 4, sd = 0.02), ncol=1))
## x3 <- rmvnorm.ivent(n/5, pardag, targets[[3]],
##                     matrix(rnorm(n/5, mean = 4, sd = 0.02), ncol=1))
## gmInt <- list(x = rbind(x1, x2, x3),
##               targets = targets,
##               target.index = target.index,
##               g = gGtrue)


###################################################
### code chunk number 19: gies-fit-plot
###################################################
score <- new("GaussL0penIntScore", gmInt$x, targets = gmInt$targets,
             target.index = gmInt$target.index)
gies.fit <- gies(score)
simy.fit <- simy(score)
par(mfrow = c(1,3))
plot(gmInt$g, main = "")           ; box(col="gray")
plot(gies.fit$essgraph, main = "") ; box(col="gray")
plot(simy.fit$essgraph, main = "") ; box(col="gray")


###################################################
### code chunk number 20: rgesExpl
###################################################
score <- new("GaussL0penObsScore", gmG8$x)
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
skel.fit <- skeleton(suffStat = suffStat, indepTest = gaussCItest,
                     alpha = 0.01, p = ncol(gmG8$x))
skel <- as(skel.fit@graph, "matrix")
ges.fit <- ges(score, fixedGaps = !skel)


###################################################
### code chunk number 21: argesExpl
###################################################
score <- new("GaussL0penObsScore", gmG8$x)
library(huge)
huge.path <- huge(gmG8$x, verbose = FALSE)
huge.fit <- huge.select(huge.path, verbose = FALSE)
cig <- as.matrix(huge.fit$refit)
ges.fit <- ges(score, fixedGaps = !cig, adaptive = "vstructures")


###################################################
### code chunk number 22: lingamExpl
###################################################
set.seed(1234)
n <- 500
## Truth: stde[1] = 0.89
eps1 <- sign(rnorm(n)) * sqrt(abs(rnorm(n)))
## Truth: stde[2] = 0.29
eps2 <- runif(n) - 0.5

## Truth: ci[2] = 3, Bpruned[1,1] = Bpruned[2,1] = 0
x2 <- 3 + eps2
## Truth: ci[1] = 7, Bpruned[1,2] = 0.9, Bpruned[2,2] = 0
x1 <- 0.9*x2 + 7 + eps1
# Truth: x1 <- x2


###################################################
### code chunk number 23: estimateLingam
###################################################
X <- cbind(x1,x2)
res <- lingam(X)
res


###################################################
### code chunk number 24: showAddBgKnowledge
###################################################
showF(addBgKnowledge)


###################################################
### code chunk number 25: explAddBgKnowledge
###################################################
amat <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3,3) ## a -- b -- c
colnames(amat) <- rownames(amat) <- letters[1:3]
## force a -> b
bg <- addBgKnowledge(gInput = amat, x = "a", y = "b")


###################################################
### code chunk number 26: explAddBgKnowledgePlot
###################################################
par(mfrow = c(1,2))
plot(as(t(amat), "graphNEL")); box(col="gray")
plot(as(t( bg ), "graphNEL")); box(col="gray")


###################################################
### code chunk number 27: ida-args
###################################################
showF(ida)


###################################################
### code chunk number 28: exIDA
###################################################
## Simulate the true DAG
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG

## simulate Gaussian data from the true DAG
n <- 10000
dat <- rmvDAG(n, myDAG)

## estimate CPDAG and PDAG -- see  help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p=p, alpha = 0.01)

## Supppose that we know the true CPDAG and covariance matrix
(l.ida.cpdag <- ida(2,5, cov(dat),
                    myCPDAG, method = "local", type = "cpdag"))
(g.ida.cpdag <- ida(2,5, cov(dat),
                    myCPDAG, method = "global", type = "cpdag"))


###################################################
### code chunk number 29: plotExIDA
###################################################
if (require(Rgraphviz)) {
  ## plot the true and estimated graphs
  par(mfrow = c(1,2))
  plot(myDAG, main = "True DAG"); box(col="gray")
  plot(pc.fit, main = "Estimated CPDAG"); box(col="gray")
}


###################################################
### code chunk number 30: jointIdaExpl1
###################################################
## Generate DAG for simulating data
p <- 6
V <- as.character(1:p)
edL <- list(
    "1" = list(edges=c(3,4), weights=c(1.1,0.3)),
    "2" = list(edges=c(6),  weights=c(0.4)),
    "3" = list(edges=c(2,4,6),weights=c(0.6,0.8,0.9)),
    "4" = list(edges=c(2),weights=c(0.5)),
    "5" = list(edges=c(1,4),weights=c(0.2,0.7)),
    "6" = NULL)
myDAG <- new("graphNEL", nodes=V, edgeL=edL,
             edgemode="directed") ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix


###################################################
### code chunk number 31: jointIdaDAG
###################################################
## Compute causal effect using true CPDAG and true cov. matrix
(resExactDAG <- jointIda(x.pos = c(1,2), y.pos = 6,
                         mcov = covTrue,
                         graphEst = myDAG,
                         technique = "RRC"))


###################################################
### code chunk number 32: jointIdaCPDAG
###################################################
(resExactCPDAG <- jointIda(x.pos = c(1,2), y.pos = 6,
                           mcov = covTrue,
                           graphEst = myCPDAG,
                           technique = "RRC"))


###################################################
### code chunk number 33: jointIdaPlot
###################################################
par(mfrow = c(1,2))
plot(myDAG)  ; box(col="gray")
plot(myCPDAG); box(col="gray")


###################################################
### code chunk number 34: backdoor-args
###################################################
showF(backdoor)


###################################################
### code chunk number 35: backdoorExCPDAG1
###################################################
p <- 6
amat <- t(matrix(c(0,0,1,1,0,1, 0,0,1,1,0,1, 0,0,0,0,1,0,
                   0,0,0,0,1,1, 0,0,0,0,0,0, 0,0,0,0,0,0), 6,6))
V <- as.character(1:6)
colnames(amat) <- rownames(amat) <- V
edL <- vector("list",length=6)
names(edL) <- V
edL[[1]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[2]] <- list(edges=c(3,4,6),weights=c(1,1,1))
edL[[3]] <- list(edges=5,weights=c(1))
edL[[4]] <- list(edges=c(5,6),weights=c(1,1))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

cov.mat <- trueCov(g)

myCPDAG <- dag2cpdag(g)
true.amat <- as(myCPDAG, "matrix")
## true.amat[true.amat != 0] <- 1


###################################################
### code chunk number 36: backdoorExpl
###################################################
par(mfrow = c(1,2))
plot(g, main = ""); box(col="gray")
plot(myCPDAG, main = ""); box(col="gray")


###################################################
### code chunk number 37: backdoorExCPDAG2
###################################################
backdoor(true.amat, 6, 3, type="cpdag")


###################################################
### code chunk number 38: gac-args
###################################################
showF(gac)


###################################################
### code chunk number 39: gacExample
###################################################
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6
## Z satisfies GAC :
gac(mFig1, x,y, z=c(2,4),    type)


###################################################
### code chunk number 40: vignette2018.Rnw:1088-1093
###################################################
mFig3a <- matrix(c(0,1,0,0, 0,0,1,1, 0,0,0,1, 0,0,1,0), 4,4)
type <- "pdag"
x <- 2; y <- 4
## Z does not satisfy GAC
str( gac(mFig3a,x,y, z=NULL, type) )## not amenable rel. to (X,Y)


###################################################
### code chunk number 41: plotGacExpl
###################################################
if (require(Rgraphviz)) {
  colnames(mFig1) <- rownames(mFig1) <- as.character(1:6)
  g1 <- as(t(mFig1), "graphNEL")
  colnames(mFig3a) <- rownames(mFig3a) <- as.character(1:4)
  g2 <- as(t(mFig3a), "graphNEL")
  ## plot the true and estimated graphs
  par(mfrow = c(1,2))
  plot(g1, main = "(a)"); box(col="gray")
  plot(g2, main = "(b)"); box(col="gray")
}


###################################################
### code chunk number 42: showAdjustment
###################################################
    showF(adjustment)


###################################################
### code chunk number 43: randDAGexplSimple
###################################################
n <- 100; d <- 3; s <- 2
myWgtFun <- function(m,mu,sd) { rnorm(m,mu,sd) }
set.seed(42)
dag1 <- randDAG(n=n, d=d, method = "er", DAG = TRUE)
dag2 <- randDAG(n=n, d=d, method = "power", DAG = TRUE)


###################################################
### code chunk number 44: evalRandDAGexpl
###################################################
w1 <- wgtMatrix(dag1)
deg1 <- apply(w1 + t(w1), 2, function(x) sum(x != 0))
max1 <- max(deg1)
mean1 <- mean(deg1)

w2 <- wgtMatrix(dag2)
deg2 <- apply(w2 + t(w2), 2, function(x) sum(x != 0))
max2 <- max(deg2)
mean2 <- mean(deg2)


###################################################
### code chunk number 45: randDAGexplExtended
###################################################
n <- 100; d <- 3; s <- 2
myWgtFun <- function(m,mu,sd) { rnorm(m,mu,sd) }
set.seed(42)
dag1 <- randDAG(n=n, d=d, method = "er", DAG = TRUE,
                weighted = TRUE, wFUN = list(myWgtFun, 0, s))
dag2 <- randDAG(n=n, d=d, method = "power", DAG = TRUE,
                weighted = TRUE, wFUN = list(myWgtFun, 0, s))


###################################################
### code chunk number 46: amatTypeCoercion
###################################################
## as(*, "amat") returns an adjacency matrix incl. its type
as(pc.gmG8, "amat")


###################################################
### code chunk number 47: amatTypeAdditional
###################################################
## INPUT: Adjacency matrix of type 'amat.cpdag'
m1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
               0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0), 6,6)
## more detailed information on the graph type needed by gac()
str( gac(m1, x=1,y=3, z=NULL, type = "dag") )


###################################################
### code chunk number 48: pcalg2dagitty
###################################################
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x) # labels aka node names

amat <- wgtMatrix(gmG8$g)
amat[amat != 0] <- 1
if(requireNamespace("dagitty")) {
    dagitty_dag1 <- pcalg2dagitty(amat,V,type="dag")
    dagitty::is.dagitty(dagitty_dag1)
}


###################################################
### code chunk number 49: Rgraphviz
###################################################
set.seed(42)
p <- 4
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)
myCPDAG <- dag2cpdag(myDAG)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
fci.fit <- fci(suffStat, indepTest=gaussCItest,
               alpha = 0.9999, p=p, doPdsep = FALSE)
if (require(Rgraphviz)) {
  par(mfrow = c(1,2))
  plot(myCPDAG); box(col="gray") ## CPDAG
  plot(fci.fit); box(col="gray") ## PAG
}


###################################################
### code chunk number 50: iplot
###################################################
data(gmG)
n <- nrow    (gmG8$ x)
V <- colnames(gmG8$ x) # labels aka node names

## estimate CPDAG
pc.fit <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = FALSE)
if (require(igraph)) {
  par(mfrow = c(1,1))
  iplotPC(pc.fit)
}


###################################################
### code chunk number 51: showAmat
###################################################
as(pc.fit,  "amat") ## Adj. matrix of type  'amat.cpdag'
as(fci.fit, "amat") ## Adj. matrix of type  'amat.pag'


###################################################
### code chunk number 52: showEdgeList
###################################################
showEdgeList(pc.fit) ## Edge list


###################################################
### code chunk number 53: simParSettings
###################################################
possible_p <- c(seq(5,10,by=1))   # paper: possible_p = c(seq(20,100,by=10))
possible_neighb_size <- c(1:3)    # paper: possible_neighb_size = c(3:10)
num_settings <-10                 # paper: num_settings = 1000
num_rep <- 2                      # paper: num_rep = 20
pb <- seq(0,1,by=0.5)             # paper: pb = seq(0,1,by=0.2)


###################################################
### code chunk number 54: simRunSim
###################################################
## helper function
revealEdge <- function(c,d,s) {
  ## cpdag, dag, selected edges to reveal
  if (!anyNA(s)) { ## something to reveal
    for (i in 1:nrow(s)) {
      c[s[i,1], s[i,2]] <- d[s[i,1], s[i,2]]
      c[s[i,2], s[i,1]] <- d[s[i,2], s[i,1]]
    }
  }
  c
}

## save results from each iteration in here:
resFin <- vector("list", num_settings)

## run simulation
for(r in 1:num_settings) {

  set.seed(r)
  ## Then we sample one setting:
  p <-  sample(possible_p,1)
  neigh <- sample(possible_neighb_size,1)
  prob <- round(neigh/(p-1),3)

  resFin[[r]] <- vector("list", num_rep)

  ## then for every setting selected we generate num_rep graphs
  for (i in 1:num_rep){

    ## get DAG
    isEmpty <-  1
    while(isEmpty){
      g <- randomDAG(p, prob)  ## true DAG
      cpdag <- dag2cpdag(g) ## true CPDAG

      ## get adjacency matrix of the CPDAG  and DAG
      cpdag.amat <- t(as(cpdag,"matrix"))
      dag.amat <- t(as(g,"matrix"))
      dag.amat[dag.amat != 0] <- 1

      ## only continue if the graph is not fully un-connected
      if (sum(dag.amat)!= 0){
        isEmpty <- 0
      }
    }

    ## choose x and y
    y <- NULL
    while(is.null(y)){
      # choose x
      x <- sample(p,1)

      ## choose y as a node connected to x but not x <- y
      skeleton <- cpdag.amat + t(cpdag.amat)
      skeleton[skeleton == 2] <- 1
      connectt <- possDe(skeleton,x, type = "pdag")
      if (length(connectt) != 1) {
        pa.x <- which(dag.amat[x,]==1 & dag.amat[,x]==0)
        ## remove x and parents of x (in the DAG) from pos.y
        pos.y <- setdiff(setdiff(connectt, pa.x), x)
        if (length(pos.y)==1){
          y <- pos.y[1]
        } else if (length(pos.y) > 0) {
          y <- sample(pos.y, 1)
        }
      }
    }

    ## calculate true effect:
    true_effect <- causalEffect(g,y,x)

    ## sample data for ida
    ## need to set nData
    nData <- 200
    dat <- rmvDAG(nData, g) ## sample data from true DAG

    ## Resulting lists, of same length as 'pb' :
    pdag.amat <-
    adjust_set <-
    result_adjust_set <-
    ida_effects <- vector("list", length(pb))

    ## for each proportion of background knowledge
    ## find a PDAG and an adjustment set relative to (x,y) in this
    ## PDAG aditionally calculate the set of possible total
    ## causal effects of x on y using ida in this PDAG
    for (j in 1:length(pb)){
      ## reveal proportion pb[j] of bg knowledge
      tmp <- ( (cpdag.amat + t(cpdag.amat))==2 ) &
        lower.tri(cpdag.amat)
      ude <- which(tmp, arr.ind = TRUE) ## undir edges
      nbg <- round(pb[j] * nrow(ude)) ## nmb of edges to reveal
      ## edges to reveal
      sele <- if (nbg>0) ude[sample(nrow(ude), nbg),,drop=FALSE] else NA
      ## reveal edges
      pdag.amat[[j]] <- revealEdge(cpdag.amat, dag.amat, sele)
      pdag.amat[[j]] <- addBgKnowledge(pdag.amat[[j]],
                                       checkInput = FALSE)

      ## find adjustment set (if it exists)
        if(requireNamespace("dagitty")) {
      adjust <- adjustment(pdag.amat[[j]],amat.type="pdag",x,y,
                           set.type="canonical")
      } else {
          adjust <- NULL
      }
      adjust_set[[j]] <- if (length(adjust) > 0) adjust$'1' else NA
      result_adjust_set[[j]] <- length(adjust) > 0
      ## ida
      ## convert to graph for ida()
      pdag.g <- as(t(pdag.amat[[j]]), "graphNEL")
      ida_effects[[j]] <- ida(x,y,cov(dat), graphEst = pdag.g,
                              method = "local", type = "pdag")

      ## for j = 1 that is when pdag.g == cpdag compare
      ## runtime of local method for CPDAGs vs. PDAGs
      if (j == 1){
        time.taken.ida <-
          system.time(ida(x,y,cov(dat), graphEst = pdag.g,
                          method = "local", type = "cpdag"))
        time.taken.bida <-
          system.time(ida(x,y,cov(dat), graphEst = pdag.g,
                          method = "local", type = "pdag"))
      }
    }

    ## save the results
    resFin[[r]][[i]] <- list(seed=r, p=p, prob=prob, neigh=neigh,
                             pb=pb, i=i, nData=nData,
                             dag.amat=dag.amat,
                             pdag.amat=pdag.amat,
                             x=x, y=y,
                             adjust_set = adjust_set,
                             result_adjust_set = result_adjust_set,
                             true_effect = true_effect,
                             ida_effects = ida_effects,
                             time.taken.ida = time.taken.ida,
                             time.taken.bida = time.taken.bida)
  }
}


###################################################
### code chunk number 55: simList2Df
###################################################
## total number of unique cpdags = num_settings*num_rep graphs
nn <- sum(sapply(resFin, length))
## make data frame with relevant summary info
nBG <- length(pb)
x <- rep(NA, nn*nBG)
df1 <- data.frame(setting=x, g=x, p=x, neigh=x, pb=x,
                  resAdj = x, idaNum = x, idaRange = x,
                  timeIda = x, timeBida = x,
                  trueEff = x)
ii <- 0
for (i1 in 1:length(resFin)) { ## settings
  nLE <- length(resFin[[i1]])
  for (i2 in 1:nLE) { ## graphs per setting
    for (i3 in 1:nBG) { ## BGK
      ii <- ii + 1
      df1[ii,"setting"] <- i1 ## List index for setting
      df1[ii,"g"] <- i2 ## List index for graph within setting
      df1[ii,"p"] <- resFin[[i1]][[i2]]$p ## Nmb nodes in graph
      ## Ave size of neighborhood
      df1[ii,"neigh"] <- resFin[[i1]][[i2]]$neigh
      ## fraction of background knowledge
      df1[ii,"pb"] <- resFin[[i1]][[i2]]$pb[i3]
      ## true if adj set exists
      df1[ii,"resAdj"] <-
        resFin[[i1]][[i2]]$result_adjust_set[[i3]]
      ## nmb unique results of ida
      df1[ii,"idaNum"] <-
        length(unique(resFin[[i1]][[i2]]$ida_effects[[i3]]))
      ## range of results of ida
      df1[ii,"idaRange"] <-
        diff(range(resFin[[i1]][[i2]]$ida_effects[[i3]]))
      ## runtime for CPDAG using option "cpdag"
      df1[ii,"timeIda"] <-
        resFin[[i1]][[i2]]$time.taken.ida[[1]]
      ## runtime for CPDAG using option "pdag"
      df1[ii,"timeBida"] <-
        resFin[[i1]][[i2]]$time.taken.bida[[1]]
      df1[ii,"trueEff"] <-
        resFin[[i1]][[i2]]$true_effect
    }
  }
}


###################################################
### code chunk number 56: simFig4Prepare
###################################################
  ## Fig 4 in paper: Fraction of identifiable effects
  ## Fraction of identifiable effects: ALL EFFECTS
  tm1 <- tapply(X=df1$resAdj, INDEX=as.factor(df1$pb), FUN = mean)
  ts1 <- tapply(X=df1$resAdj, INDEX=as.factor(df1$pb),
                FUN = function(x) sd(x)/sqrt(length(x)))
  ## Fraction of identifiable effects: add means for
  ## only NON-ZERO EFFECTS
  dfNZ <- subset(df1, subset = (trueEff!=0) )
  tm <- c(tm1,tapply(X=dfNZ$resAdj, INDEX=as.factor(dfNZ$pb),
                     FUN = mean))
  ts <- c(ts1,tapply(X=dfNZ$resAdj, INDEX=as.factor(dfNZ$pb),
                     FUN = function(x) sd(x)/sqrt(length(x))))

  dfID <- data.frame(pb = as.factor(names(tm)), fit = tm, se = ts,
                     TrueEffect =
                       as.factor(c(rep("All", length(tm)/2),
                                   rep("Non-zero", length(tm)/2))))


###################################################
### code chunk number 57: simFig4
###################################################
if(require(ggplot2)) {
  k <- ggplot(dfID, aes(pb, fit, ymin = fit-se,
                        ymax = fit+se, col = TrueEffect))
  k + geom_pointrange() +
    xlab("Proportion of background knowledge") +
    ylab("Fraction of identifiable effects via adjustment") +
    theme(legend.position = c(0.9,0.1),
          axis.text=element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text=element_text(size = 14),
          legend.title=element_text(size = 14))
}


###################################################
### code chunk number 58: simFig5prepare
###################################################
## use dfNU2: settings where effect is NOT unique given zero bg knowledge
nn <- length(pb)
idx <- rep(seq(1,nrow(df1), by = nn), each = nn) ## pb=0 rows
nmbIda <- df1$idaNum[idx]
dfNU2 <- df1[nmbIda > 1,]

bnTmp <- cut(x=dfNU2$idaNum, breaks = c(0,1,2,3,4,1e9),
             labels = c("1", "2", "3", "4", "5+"))
dfNU2$idaNumF <- factor(bnTmp, levels = levels(bnTmp)[5:1])
df3 <- dfNU2[,c("pb", "idaNumF")]
df3$idx <- 1:nrow(df3)

df3N <- aggregate(idx ~ pb + idaNumF, data = df3, FUN = length)
df3N$idaNumF <- droplevels(df3N$idaNumF)


###################################################
### code chunk number 59: simFig5
###################################################
ggplot(df3N, aes(x = pb, y=idx, fill = idaNumF)) +
  geom_bar(stat = "identity") +
  ylab("Number of simulation settings") +
  xlab("Proportion of background knowledge")+
  theme(axis.text =  element_text(size = 14),
        axis.title=  element_text(size = 14),
        legend.text= element_text(size = 14),
        legend.title=element_text(size = 14)) +
  guides(fill=guide_legend(title="#effects"))



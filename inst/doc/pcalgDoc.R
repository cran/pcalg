### R code from vignette source 'pcalgDoc.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
op.orig <-
options(SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        width = 75, digits = 5,
        ## JSS         : prompt = "R> "
        ## Good looking:
        prompt = "> ", continue = "  "
        )


###################################################
### code chunk number 2: def-gmG (eval = FALSE)
###################################################
## ## Used to generate the  'gmG'  Gaussian data originally:
## require("pcalg")
## set.seed(40)
## p <- 8
## n <- 5000
## gGtrue <- randomDAG(p, prob = 0.3)
## gmG <- list(x = rmvDAG(n, gGtrue), g = gGtrue)


###################################################
### code chunk number 3: exIntro1
###################################################
library("pcalg")
data("gmG")


###################################################
### code chunk number 4: Iplot (eval = FALSE)
###################################################
## stopifnot(require(Rgraphviz))# needed for all our graph plots
## par(mfrow = c(1,2))
## plot(gmG8$g, main = "") ; plot(pc.gmG, main = "")


###################################################
### code chunk number 5: exIntro
###################################################
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
pc.gmG <- pc(suffStat, indepTest = gaussCItest,
             p = ncol(gmG8$x), alpha = 0.01)


###################################################
### code chunk number 6: exIntroPlot
###################################################
stopifnot(require(Rgraphviz))# needed for all our graph plots
par(mfrow = c(1,2))
plot(gmG8$g, main = "") ; plot(pc.gmG, main = "")


###################################################
### code chunk number 7: exIntro2
###################################################
ida(1, 6, cov(gmG8$x), pc.gmG@graph)


###################################################
### code chunk number 8: exIntro3
###################################################
idaFast(1, c(4,5,6), cov(gmG8$x), pc.gmG@graph)


###################################################
### code chunk number 9: skeleton-args
###################################################
showF <- function(f, width = 80) {
    ## 'width': larger than default on purpose:
    nam <- deparse(substitute(f))
    stopifnot(is.function(f))
    attr(f, "source") <- NULL # if ...
    attr(f, "srcref") <- NULL
    ll <- capture.output(str(f, width=width))
    ll[1] <- sub("function *", nam, ll[1])
    writeLines(ll)
}
showF(skeleton)


###################################################
### code chunk number 10: skelExpl1Plot
###################################################
## using  data("gmG", package="pcalg")
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
skel.gmG <- skeleton(suffStat, indepTest = gaussCItest,
                   p = ncol(gmG8$x), alpha = 0.01)
par(mfrow = c(1,2))
plot(gmG8$g, main = ""); plot(skel.gmG, main = "")


###################################################
### code chunk number 11: skelExp2Plot
###################################################
data("gmD")
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
skel.gmD <- skeleton(suffStat, indepTest = disCItest,
                     p = ncol(gmD$x), alpha = 0.01)
par(mfrow= 1:2); plot(gmD$g, main = ""); plot(skel.gmD, main = "")


###################################################
### code chunk number 12: pc-args
###################################################
showF(pc)


###################################################
### code chunk number 13: pcExpl-plot
###################################################
suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
pc.fit <- pc(suffStat, indepTest=gaussCItest, p = ncol(gmG8$x), alpha = 0.01)
par(mfrow= c(1,2));  plot(gmG8$g, main = ""); plot(pc.fit, main = "")


###################################################
### code chunk number 14: obs-score-args (eval = FALSE)
###################################################
## score <- new("GaussL0penObsScore", data = matrix(1, 1, 1),
##     lambda = 0.5*log(nrow(data)), intercept = FALSE, use.cpp = TRUE, ...)


###################################################
### code chunk number 15: ges-args
###################################################
showF(ges)


###################################################
### code chunk number 16: gesExpl-plot
###################################################
score <- new("GaussL0penObsScore", gmG8$x)
ges.fit <- ges(score)
par(mfrow=1:2); plot(gmG8$g, main = ""); plot(ges.fit$essgraph, main = "")


###################################################
### code chunk number 17: fci-args
###################################################
showF(fci, width=75)


###################################################
### code chunk number 18: fciExpl-plot
###################################################
data("gmL")
suffStat1 <- list(C = cor(gmL$x), n = nrow(gmL$x))
pag.est <- fci(suffStat1, indepTest = gaussCItest,
               p = ncol(gmL$x), alpha = 0.01, labels = as.character(2:5))
par(mfrow = 1:2); plot(gmL$g, main = ""); plot(pag.est)


###################################################
### code chunk number 19: rfci-args
###################################################
showF(rfci)


###################################################
### code chunk number 20: def-rfciExpl-plot (eval = FALSE)
###################################################
## data("gmL")
## suffStat1 <- list(C = cor(gmL$x), n = nrow(gmL$x))
## pag.est <- rfci(suffStat1, indepTest = gaussCItest,
##                 p = ncol(gmL$x), alpha = 0.01, labels = as.character(2:5))


###################################################
### code chunk number 21: int-score-args (eval = FALSE)
###################################################
## score <- new("GaussL0penIntScore", data = matrix(1, 1, 1),
##     targets = list(integer(0)), target.index = rep(as.integer(1), nrow(data)),
##     lambda = 0.5*log(nrow(data)), intercept = FALSE, use.cpp = TRUE, ...)


###################################################
### code chunk number 22: gies-args
###################################################
showF(gies)


###################################################
### code chunk number 23: def-gmInt (eval = FALSE)
###################################################
## ## Used to generate the  'gmInt'  Gaussian data originally:
## set.seed(40)
## p <- 8
## n <- 5000
## gGtrue <- randomDAG(p, prob = 0.3)
## nodes(gGtrue) <- c("Author", "Bar", "Ctrl", "Goal", "V5", "V6", "V7", "V8")
## pardag <- as(gGtrue, "GaussParDAG")
## pardag$set.err.var(rep(1, p))
## targets <- list(integer(0), 3, 5)
## target.index <- c(rep(1, 0.6*n), rep(2, n/5), rep(3, n/5))
## 
## x1 <- rmvnorm.ivent(0.6*n, pardag)
## x2 <- rmvnorm.ivent(n/5, pardag, targets[[2]],
##   matrix(rnorm(n/5, mean = 4, sd = 0.02), ncol = 1))
## x3 <- rmvnorm.ivent(n/5, pardag, targets[[3]],
##   matrix(rnorm(n/5, mean = 4, sd = 0.02), ncol = 1))
## gmInt <- list(x = rbind(x1, x2, x3),
##   targets = targets,
##   target.index = target.index,
##   g = gGtrue)


###################################################
### code chunk number 24: load-gmInt
###################################################
data(gmInt)
n.tot <- length(gmInt$target.index)
n.obs <- sum(gmInt$target.index == 1)
n3 <- sum(gmInt$target.index == 2)
n5 <- sum(gmInt$target.index == 3)


###################################################
### code chunk number 25: load-gmInt2
###################################################
data(gmInt)


###################################################
### code chunk number 26: gies-fit-plot
###################################################
score <- new("GaussL0penIntScore", gmInt$x, targets = gmInt$targets,
             target.index = gmInt$target.index)
gies.fit <- gies(score)
simy.fit <- simy(score)
par(mfrow = c(1, 3)) ; plot(gmInt$g, main = "")
plot(gies.fit$essgraph, main = "")
plot(simy.fit$essgraph, main = "")


###################################################
### code chunk number 27: def-gmI (eval = FALSE)
###################################################
## set.seed(123)
## p <- 7
## n <- 10000
## myDAG <- randomDAG(p, prob = 0.2)
## datI <- rmvDAG(n, myDAG)
## gmI <- list(x = datI, g = myDAG)


###################################################
### code chunk number 28: idaExpl1
###################################################
data("gmI")
suffStat <- list(C = cor(gmI$x), n = nrow(gmI$x))
pc.gmI <- pc(suffStat, indepTest=gaussCItest,
             p = ncol(gmI$x), alpha = 0.01)


###################################################
### code chunk number 29: idaExpl2
###################################################
par(mfrow = c(1,2))
plot(gmI$g, main = "")
plot(pc.gmI, main = "")


###################################################
### code chunk number 30: idaExpl3
###################################################
am.pdag <- wgtMatrix(pc.gmI@graph)
ad <- pdag2allDags(am.pdag)$dags
gDag <- vector("list", nrow(ad))
for (i in 1:nrow(ad)) gDag[[i]] <- as(matrix(ad[i, ], 7, 7), "graphNEL")
par(mfrow = c(3,2))
for (i in 1:6) plot(gDag[[i]], main = paste("DAG",i))


###################################################
### code chunk number 31: plot-6DAGS
###################################################
sfsmisc::mult.fig(6)
for (i in 1:6) plot(gDag[[i]], main = paste("DAG",i))


###################################################
### code chunk number 32: idaExpl4
###################################################
ida(2, 5, cov(gmI$x), pc.gmI@graph, method = "global", verbose = FALSE)


###################################################
### code chunk number 33: ida-args
###################################################
showF(ida)


###################################################
### code chunk number 34: idaExpl5
###################################################
ida(2,5, cov(gmI$x), pc.gmI@graph, method = "local")


###################################################
### code chunk number 35: idaFast-args
###################################################
showF(idaFast)


###################################################
### code chunk number 36: ida-idaFast
###################################################
(eff.est1 <- ida(2,5, cov(gmI$x), pc.gmI@graph, method="local"))
(eff.est2 <- ida(2,6, cov(gmI$x), pc.gmI@graph, method="local"))
(eff.est3 <- ida(2,7, cov(gmI$x), pc.gmI@graph, method="local"))

(eff.estF <- idaFast(2, c(5,6,7), cov(gmI$x), pc.gmI@graph))


###################################################
### code chunk number 37: backdoor-args
###################################################
showF(backdoor)


###################################################
### code chunk number 38: backdoorExCPDAG1
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
### code chunk number 39: backdoorExpl
###################################################
par(mfrow = c(1,2))
plot(g, main = "")
plot(myCPDAG, main = "")


###################################################
### code chunk number 40: backdoorExCPDAG2
###################################################
backdoor(true.amat, 6, 3, type="cpdag")


###################################################
### code chunk number 41: turn-off-plus
###################################################
options(continue = " ") # MM: so we don't get the "+ " continuation lines


###################################################
### code chunk number 42: myCItest
###################################################
myCItest <- function(x,y,S, suffStat) {
  if (length(S) == 0) {
    x. <- suffStat[,x]
    y. <- suffStat[,y]
  } else {
    rxy <- resid(lm.fit(y= suffStat[,c(x,y)], x= cbind(1, suffStat[,S])))
    x. <- rxy[,1];  y. <- rxy[,2]
  }
  cor.test(x., y.)$p.value
}


###################################################
### code chunk number 43: gaussCItest-ex
###################################################
suffStat <- list(C = cor(gmG8$x), n = 5000)
pc.gmG <- pc(suffStat, indepTest=gaussCItest, p = 8, alpha = 0.01)


###################################################
### code chunk number 44: myCItest-def-plot (eval = FALSE)
###################################################
## pc.myfit <- pc(suffStat = gmG8$x, indepTest = myCItest,
##                p = 8, alpha = 0.01)
## par(mfrow = c(1,2)); plot(pc.gmG, main = ""); plot(pc.myfit, main = "")


###################################################
### code chunk number 45: myCItest-ex-plot
###################################################
pc.myfit <- pc(suffStat = gmG8$x, indepTest = myCItest,
               p = 8, alpha = 0.01)
par(mfrow = c(1,2)); plot(pc.gmG, main = ""); plot(pc.myfit, main = "")


###################################################
### code chunk number 46: time-tests (eval = FALSE)
###################################################
## system.time(for(i in 1:10)
##             pc.fit <- pc(suffStat, indepTest=gaussCItest, p = 8, alpha = 0.01))
##       ##  User      System verstrichen
##       ## 0.593       0.000       0.594
## system.time(for(i in 1:10)
##             pc.myfit <- pc(gmG8$x,  indepTest = myCItest,  p = 8, alpha = 0.01))
## ## Using  resid(lm(..)) twice:
##      ##   User      System verstrichen
##      ## 44.864       0.007      44.937
## ## Using resid(lm.fit(..)):
##      ## 10.550       0.067      10.632


###################################################
### code chunk number 47: finalizing
###################################################
options(op.orig)



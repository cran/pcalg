## Translate amat as describes in amatType to dagitty object

library(pcalg)
library(dagitty)
suppressWarnings(RNGversion("3.5.0"))
doExtras <- pcalg:::doExtras()

res <- rep(FALSE, 10)
####################
## Test DAG 1
####################
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x) # labels aka node names

amat <- wgtMatrix(gmG8$g)
amat[amat != 0] <- 1
dagitty_dag1 <- pcalg2dagitty(amat,V,type="dag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow=c(1,2))
## plot(gmG8$g, main = "True DAG")
## plot(dagitty:::graphLayout(dagitty_dag1))

res[1] <- (dagitty_dag1 == "dag {\nAuthor\nBar\nCtrl\nGoal\nV5\nV6\nV7\nV8\nAuthor -> Bar\nAuthor -> V6\nAuthor -> V8\nBar -> Ctrl\nBar -> V5\nV5 -> V6\nV5 -> V8\nV6 -> V7\n}\n")

#############
## Test DAG 2
#############
set.seed(123)
p <- 10
V <- sample(LETTERS, p)
g <- pcalg::randomDAG(p,prob=0.3, V = V)

amat <- wgtMatrix(g)
amat[amat != 0] <- 1
dagitty_dag2 <- pcalg2dagitty(amat,V,type="dag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow=c(1,2))
## plot(g, main = "True DAG")
## plot(dagitty:::graphLayout(dagitty_dag2))

res[2] <- (dagitty_dag2 == "dag {\nA\nH\nJ\nK\nQ\nT\nU\nW\nX\nZ\nA -> Q\nH -> A\nH -> K\nH -> Q\nH -> T\nH -> Z\nJ -> W\nT -> A\nT -> Q\nT -> X\nU -> Q\nU -> W\nU -> X\nW -> K\n}\n")

###############
## Test CPDAG 1
###############
data(gmG)
n <- nrow(gmG8$ x)
V <- colnames(gmG8$ x) # labels aka node names

## estimate CPDAG
pc.fit <- pc(suffStat = list(C = cor(gmG8$x), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = FALSE)
amat <- as(pc.fit, "amat")
dagitty_cpdag1 <- pcalg2dagitty(amat,V,type="cpdag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow = c(1,2))
## plot(pc.fit)
## plot(dagitty:::graphLayout(dagitty_cpdag1))

res[3] <- (dagitty_cpdag1 == "pdag {\nAuthor\nBar\nCtrl\nGoal\nV5\nV6\nV7\nV8\nAuthor -- Bar\nAuthor -> V6\nAuthor -> V8\nBar -- Ctrl\nBar -> V5\nV5 -> V6\nV5 -> V8\nV6 -> V7\n}\n")

stopifnot(all(res[1:3]))

if (doExtras) {
#############
## Test CPDAG 2
#############
set.seed(135)
p <- 10
V <- sample(LETTERS, p)
g <- dag2cpdag(pcalg::randomDAG(p,prob=0.3, V = V))

amat <- wgtMatrix(g)
amat[amat != 0] <- 1
dagitty_cpdag2 <- pcalg2dagitty(amat,V,type="cpdag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow=c(1,2))
## plot(g)
## plot(dagitty:::graphLayout(dagitty_cpdag2))

res[4] <- (dagitty_cpdag2 == "pdag {\nA\nB\nH\nI\nJ\nK\nO\nS\nV\nX\nA -- I\nA -- J\nA -- V\nA -> B\nA -> O\nH -- I\nH -> B\nI -> B\nJ -- K\nK -> B\nS -- X\nS -> O\nV -> B\n}\n")

#############
## Test MAG 1
#############
amat <- matrix(c(0,2,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
V <- LETTERS[1:4]
colnames(amat) <- rownames(amat) <- V
## plotAG(amat)
dagitty_mag1 <- pcalg2dagitty(amat,V,type="mag")
res[5] <- (dagitty_mag1 == "mag {\nA\nB\nC\nD\nA <-> B\nB -> C\nB -> D\nC -> D\n}\n")

#############
## Test MAG 2
#############
set.seed(78)
p <- 8
g <- pcalg::randomDAG(p, prob = 0.4)
## Compute the true covariance and then correlation matrix of g:
true.corr <- cov2cor(trueCov(g))

## define nodes 2 and 6 to be latent variables
L <- c(2,6)

## Find PAG
## As dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
true.pag <- dag2pag(suffStat = list(C= true.corr, n= 10^9),
                    indepTest= gaussCItest, graph=g, L=L, alpha= 0.9999)

## find a valid MAG such that no additional edges are directed into
(amat <- pag2magAM(true.pag@amat, 4)) # -> the adj.matrix of the MAG
## plotAG(amat)
V <- colnames(amat)
dagitty_mag2 <- pcalg2dagitty(amat,V,type="mag")
res[6] <- (dagitty_mag2 == "mag {\n1\n2\n3\n4\n5\n6\n1 -> 4\n1 -> 5\n1 -> 6\n2 -> 5\n3 -> 4\n3 -> 6\n4 -> 5\n4 -> 6\n5 <-> 6\n}\n")

#############
## Test PAG 1
#############
mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2,
                   0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
V <- c("V1", "V2", "X", "Y", "V4", "V3")
colnames(mFig4b) <- rownames(mFig4b) <- V
## plotAG(mFig4b)

dagitty_pag1 <- pcalg2dagitty(mFig4b,V,type="pag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow=c(1,2))
## plot(g)
## plot(dagitty:::graphLayout(dagitty_cpdag2))

res[7] <- (dagitty_pag1 == "pag {\nV1\nV2\nV3\nV4\nX\nY\nV1 @-> X\nV2 @-> X\nV3 -> Y\nV3 <-> V4\nV3 <-> X\nV4 -> Y\nX -> V4\n}\n")

#############
## Test PAG 2
#############
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- pcalg::randomDAG(p, prob = 0.4)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
fm <- fci(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p, doPdsep = FALSE)

amat <- as(fm, "amat")
V <- colnames(amat)
dagitty_pag2 <- pcalg2dagitty(amat,V,type="pag")

res[8] <- (dagitty_pag2 == "pag {\n1\n2\n3\n4\n5\n6\n7\n1 -> 7\n1 @-> 5\n1 @-> 6\n1 @-@ 3\n2 -> 7\n2 @-> 5\n2 @-> 6\n3 -> 7\n3 @-> 5\n3 @-> 6\n3 @-@ 4\n4 @-> 6\n5 @-> 7\n6 -> 7\n}\n")

#################
## Test empty DAG
#################
set.seed(123)
p <- 10
V <- sample(LETTERS, p)
g <- pcalg::randomDAG(p,prob=0, V = V)

amat <- wgtMatrix(g)
amat[amat != 0] <- 1
dagitty_dagE <- pcalg2dagitty(amat,V,type="dag")
## Use dagitty:::graphLayout instead of just graphLayout
## because Rgraphviz package that R uses has a function w the same name
## par(mfrow=c(1,2))
## plot(g, main = "True DAG")
## plot(dagitty:::graphLayout(dagitty_dagE))

res[9] <- (dagitty_dagE == "dag {\nA\nH\nJ\nK\nQ\nT\nU\nW\nX\nZ\n\n}\n")

#################
## Test empty PAG
#################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- pcalg::randomDAG(p, prob = 0)

## find skeleton and PAG using the FCI algorithm
suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
fm <- fci(suffStat, indepTest=gaussCItest,
           alpha = 0.9999, p=p, doPdsep = FALSE)

amat <- as(fm, "amat")
V <- colnames(amat)
dagitty_pagE <- pcalg2dagitty(amat,V,type="pag")

res[10] <- (dagitty_pagE == "pag {\n1\n2\n3\n4\n5\n6\n7\n\n}\n")

stopifnot(all(res))

########################################################
## Test via comparison of gac() and isAdjustmentSet() ##
########################################################
gacVSdagitty <- function(amat, x, y ,z, V, type) {
    ## x,y,z: col positions as used in GAC
    ## Result: TRUE is result is equal
    typeDG <- switch(type,
                     dag = "dag",
                     cpdag = "cpdag",
                     mag = "mag",
                     pag = "pag")
    
    dgRes <- pcalg2dagitty(amat, V, type = typeDG)
    Exp <- V[x]; Out <- V[y]; Z <- V[z]
    gacRes <- gac(amat,x,y, z, type)$gac 
    dgRes <- isAdjustmentSet(x = dgRes, Z = Z, exposure = Exp, outcome = Out)
    (gacRes == dgRes)
}

## CPDAG 1: Paper Fig 1
## mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
##                   0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
## V <- as.character(1:nrow(mFig1))
## colnames(mFig1) <- rownames(mFig1) <- V

## typeGAC <- "cpdag"
## x <- 3; y <- 6
## z <- c(2,4); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- c(4,5); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- c(4,2,1); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- c(4,5,1); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- c(4,2,5); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- c(4,2,5,1); gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- 2; gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)
## z <- NULL; gacVSdagitty(amat = mFig1, x=x, y=y, z=z, V=V, type=typeGAC)

xx <- TRUE
##################################################
## DAG / CPDAG
##################################################
## CPDAG 1: Paper Fig 1
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6

V <- as.character(1:ncol(mFig1))
rownames(mFig1) <- colnames(mFig1) <- V

xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(2,4), V=V, type)
xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(4,5), V=V, type)
xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(4,2,1), V=V, type)
xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(4,5,1), V=V, type)
xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(4,2,5), V=V, type)
xx <- xx &  gacVSdagitty(mFig1,x,y, z=c(4,2,5,1), V=V, type)
xx <- xx & gacVSdagitty(mFig1,x,y, z= 2,    V=V, type)
xx <- xx & gacVSdagitty(mFig1,x,y, z= NULL, V=V, type)

## CPDAG 2: Paper Fig 5a
mFig5a <- matrix(c(0,1,0,0,0, 1,0,1,0,0, 0,0,0,0,1, 0,0,1,0,0, 0,0,0,0,0), 5,5)
V <- as.character(1:ncol(mFig5a))
rownames(mFig5a) <- colnames(mFig5a) <- V

type <- "cpdag"
x <- c(1,5); y <- 4
xx <- xx &  gacVSdagitty(mFig5a, x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSdagitty(mFig5a, x,y, z= 2,     V=V, type)

## DAG 1 from Marloes' Talk
mMMd1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
                  0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
V <- as.character(1:ncol(mMMd1))
rownames(mMMd1) <- colnames(mMMd1) <- V

type <- "dag"
x <- 1; y <- 3
xx <- xx &  gacVSdagitty(mMMd1, x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(mMMd1, x,y, z= 2, V=V, type)
xx <- xx &  gacVSdagitty(mMMd1, x,y, z= 4, V=V, type)
xx <- xx & gacVSdagitty(mMMd1, x,y, z= 5, V=V, type)
xx <- xx & gacVSdagitty(mMMd1, x,y, z= 6, V=V, type)
xx <- xx & gacVSdagitty(mMMd1, x,y, z=c(4,5), V=V, type)

## DAG 2 from Marloes' Talk
mMMd2 <- matrix(c(0,1,0,1,0,0, 0,0,0,0,0,0, 0,1,0,0,1,0,
                  0,0,0,0,1,0, 0,0,0,0,0,1, 0,0,0,0,0,0), 6,6)
V <- as.character(1:ncol(mMMd2))
rownames(mMMd2) <- colnames(mMMd2) <- V

type <- "dag"
x <- 4; y <- 6
xx <- xx &  gacVSdagitty(mMMd2, x,y, z= 1, V=V, type)
xx <- xx &  gacVSdagitty(mMMd2, x,y, z= 3, V=V, type)
xx <- xx & gacVSdagitty(mMMd2, x,y, z= 5, V=V, type)
xx <- xx & gacVSdagitty(mMMd2, x,y, z=c(1,5), V=V, type)
xx <- xx &  gacVSdagitty(mMMd2, x,y, z=c(1,2), V=V, type)
xx <- xx &  gacVSdagitty(mMMd2, x,y, z=c(1,3), V=V, type)
xx <- xx & gacVSdagitty(mMMd2, x,y, z= 2, V=V, type)

##################################################
## PAG
##################################################
type <- "pag"
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
V <- as.character(1:ncol(mFig3a))
rownames(mFig3a) <- colnames(mFig3a) <- V
xx <- xx & gacVSdagitty(mFig3a, x=2,      y=4, z=NULL,   V=V, type)

mFig3b <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
V <- as.character(1:ncol(mFig3b))
rownames(mFig3b) <- colnames(mFig3b) <- V
xx <- xx & gacVSdagitty(mFig3b, x=2,      y=4, z=NULL,   V=V, type)

mFig3c <- matrix(c(0,3,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
V <- as.character(1:ncol(mFig3c))
rownames(mFig3c) <- colnames(mFig3c) <- V
xx <- xx &  gacVSdagitty(mFig3c, x=2,      y=4, z=NULL,   V=V, type)

mFig4a <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,3,3,2,
                   0,0,2,0,2,2, 0,0,2,1,0,2, 0,0,1,3,3,0), 6,6)
V <- as.character(1:ncol(mFig4a))
rownames(mFig4a) <- colnames(mFig4a) <- V
xx <- xx & gacVSdagitty(mFig4a, x=3,      y=4, z=NULL,   V=V, type)
xx <- xx &  gacVSdagitty(mFig4a, x=3,      y=4, z= 6,     V=V, type)
xx <- xx &  gacVSdagitty(mFig4a, x=3,      y=4, z=c(1,6), V=V, type)
xx <- xx &  gacVSdagitty(mFig4a, x=3,      y=4, z=c(2,6), V=V, type)
xx <- xx &  gacVSdagitty(mFig4a, x=3,      y=4, z=c(1,2,6), V=V, type)

mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2,
                   0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
V <- as.character(1:ncol(mFig4b))
rownames(mFig4b) <- colnames(mFig4b) <- V
xx <- xx & gacVSdagitty(mFig4b, x=3,      y=4, z=NULL,   V=V, type)
xx <- xx & gacVSdagitty(mFig4b, x=3,      y=4, z= 6,     V=V, type)
xx <- xx & gacVSdagitty(mFig4b, x=3,      y=4, z=c(5,6), V=V, type)

mFig5b <- matrix(c(0,1,0,0,0,0,0, 2,0,2,3,0,3,0, 0,1,0,0,0,0,0, 0,2,0,0,3,0,0,
                   0,0,0,2,0,2,3, 0,2,0,0,2,0,0, 0,0,0,0,2,0,0), 7,7)
V <- as.character(1:ncol(mFig5b))
rownames(mFig5b) <- colnames(mFig5b) <- V
xx <- xx & gacVSdagitty(mFig5b, x=c(2,7), y=6, z=NULL,   V=V, type)
xx <- xx &  gacVSdagitty(mFig5b, x=c(2,7), y=6, z=c(4,5), V=V, type)
xx <- xx &  gacVSdagitty(mFig5b, x=c(2,7), y=6, z=c(4,5,1), V=V, type)
xx <- xx &  gacVSdagitty(mFig5b, x=c(2,7), y=6, z=c(4,5,3), V=V, type)
xx <- xx &  gacVSdagitty(mFig5b, x=c(2,7), y=6, z=c(1,3,4,5), V=V, type)

## PAG in Marloes' talk
mMMp <- matrix(c(0,0,0,3,2,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 2,0,0,0,0,3,2,
                 3,2,2,0,0,0,3, 0,0,0,2,0,0,0, 0,0,0,2,2,0,0), 7,7)
V <- as.character(1:ncol(mMMp))
rownames(mMMp) <- colnames(mMMp) <- V

x <- c(5,6); y <- 7
xx <- xx & gacVSdagitty(mMMp, x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(mMMp, x,y, z= 1,   V=V, type)
xx <- xx & gacVSdagitty(mMMp, x,y, z= 4,   V=V, type)
xx <- xx & gacVSdagitty(mMMp, x,y, z= 2,   V=V, type)
xx <- xx & gacVSdagitty(mMMp, x,y, z= 3,   V=V, type)
xx <- xx & gacVSdagitty(mMMp, x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSdagitty(mMMp, x,y, z=c(1,4), V=V, type)
xx <- xx &  gacVSdagitty(mMMp, x,y, z=c(1,4,2), V=V, type)
xx <- xx &  gacVSdagitty(mMMp, x,y, z=c(1,4,3), V=V, type)
xx <- xx &  gacVSdagitty(mMMp, x,y, z=c(1,4,2,3), V=V, type)

##################################################
## V=V, type = "pag" -- Tests from Ema
##################################################
type <- "pag"
pag.m <- readRDS(system.file("external/gac-pags.rds", package="pcalg"))
m1 <- pag.m[["m1"]]
V <- colnames(m1)
x <- 6; y <- 9
xx <- xx & gacVSdagitty(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=1, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=3, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=4, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,7,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,5,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,5,7,8), V=V, type)

x <- c(6,8); y <- 9
xx <- xx & gacVSdagitty(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=1, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=3, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=4, V=V, type)
xx <- xx &  gacVSdagitty(m1,x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSdagitty(m1,x,y, z=c(2,3,4), V=V, type)
xx <- xx &  gacVSdagitty(m1,x,y, z=c(2,3,7), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,5), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,3,5,7), V=V, type)

x <- 3; y <- 1
xx <- xx & gacVSdagitty(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=4, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=5, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=6, V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,6), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,7,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,5,8), V=V, type)
xx <- xx & gacVSdagitty(m1,x,y, z=c(2,5,7,8), V=V, type)

m2 <- pag.m[["m2"]]
V <- colnames(m2)
x <- 3; y <-1
xx <- xx & gacVSdagitty(m2,x,y, z=NULL, V=V, type)
xx <- xx &  gacVSdagitty(m2,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=4, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=8, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=9, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=c(2,8,9), V=V, type)
xx <- xx &  gacVSdagitty(m2,x,y, z=c(2,5), V=V, type)

x <- c(3,9); y <- 1
xx <- xx & gacVSdagitty(m2,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=4, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=8, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=9, V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=c(2,8,9), V=V, type)
xx <- xx & gacVSdagitty(m2,x,y, z=c(2,5), V=V, type)

m3 <- pag.m[["m3"]]
V <- colnames(m3)
x <- 1; y <- 9
xx <- xx & gacVSdagitty(m3,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=3, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=5, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=7, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=8, V=V, type)
xx <- xx &  gacVSdagitty(m3,x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSdagitty(m3,x,y, z=c(5,7), V=V, type)

x <- 1; y <- 8
xx <- xx & gacVSdagitty(m3,x,y, z=NULL, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=2, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=3, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=5, V=V, type)
xx <- xx &  gacVSdagitty(m3,x,y, z=7, V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=9, V=V, type)
xx <- xx &  gacVSdagitty(m3,x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSdagitty(m3,x,y, z=c(5,9), V=V, type)

if (!xx) {
    stop("Problem when testing function gacVSdagitty.")
} else {
    message("OK, no issues were found.")
}
}

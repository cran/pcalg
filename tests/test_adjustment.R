library(pcalg)
(doExtras <- pcalg:::doExtras())

## Minimalistic CRAN checks

## Test 1 ############################
## Test that "no adjustment set" and "empty adjustment set" are distinguished properly
x <- 1; y <- 2
cpdag <- matrix(c(0,1,1,0),2,2) ## 1 --- 2 => no adj set
dag <- matrix(c(0,1,0,0),2,2) ## 1 --> 2 => empty adj set

adjC <- adjustment(amat = cpdag, amat.type = "cpdag", x = 1, y = 2, set.type = "canonical")
adjD <- adjustment(amat = dag, amat.type = "dag", x = 1, y = 2, set.type = "canonical")
adjP <- adjustment(amat = dag, amat.type = "pdag", x = 1, y = 2, set.type = "canonical")

stopifnot(!identical(adjC, adjD), identical(adjD, adjP) )

## Test 2 ###############################
gacVSadj <- function(amat, x, y ,z, V, type) {
  ## gac(z) is TRUE IFF z is returned by adjustment()
  ## x,y,z: col positions as used in GAC
  ## Result: TRUE is result is equal
  typeDG <- switch(type,
                   dag = "dag",
                   cpdag = "cpdag",
                   mag = "mag",
                   pag = "pag")
  gacRes <- gac(amat,x,y, z, type)$gac 
  adjRes <- adjustment(amat = amat, amat.type = typeDG, x = x, y = y, set.type = "all")
  if (gacRes) { ## z is valid adj set
    res <- any(sapply(adjRes, function(xx) setequal(z, xx)))
  } else { ## z is not valid adj set
    res <- all(!sapply(adjRes, function(xx) setequal(z, xx)))
  }
  res
}

xx <- TRUE

## CPDAG 1: Paper Fig 1
mFig1 <- matrix(c(0,1,1,0,0,0, 1,0,1,1,1,0, 0,0,0,0,0,1,
                  0,1,1,0,1,1, 0,1,0,1,0,1, 0,0,0,0,0,0), 6,6)
type <- "cpdag"
x <- 3; y <- 6

V <- as.character(1:ncol(mFig1))
rownames(mFig1) <- colnames(mFig1) <- V

xx <- xx &  gacVSadj(mFig1,x,y, z=c(2,4), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,5), V=V, type)

type <- "pag"
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
V <- as.character(1:ncol(mFig3a))
rownames(mFig3a) <- colnames(mFig3a) <- V
xx <- xx & gacVSadj(mFig3a, x=2,      y=4, z=NULL,   V=V, type)

## DAG 1 from Marloes' Talk
mMMd1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
                  0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
V <- as.character(1:ncol(mMMd1))
rownames(mMMd1) <- colnames(mMMd1) <- V

type <- "dag"
x <- 1; y <- 3
xx <- xx &  gacVSadj(mMMd1, x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(mMMd1, x,y, z= 2, V=V, type)

if (!xx) {
  stop("Problem when testing function gacVSadj.")
} else {
  message("OK, no issues were found.")
}

############################################################
## Extensive checks
############################################################
if (doExtras) {

## Test that "no adjustment set" and "empty adjustment set" are distinguished properly
x <- 1; y <- 2
cpdag <- matrix(c(0,1,1,0),2,2) ## 1 --- 2 => no adj set
dag <- matrix(c(0,1,0,0),2,2) ## 1 --> 2 => empty adj set

adjC <- adjustment(amat = cpdag, amat.type = "cpdag", x = 1, y = 2, set.type = "canonical")
adjD <- adjustment(amat = dag, amat.type = "dag", x = 1, y = 2, set.type = "canonical")
adjP <- adjustment(amat = dag, amat.type = "pdag", x = 1, y = 2, set.type = "canonical")

stopifnot(!identical(adjC, adjD), identical(adjD, adjP) )

adjCAll <- adjustment(amat = cpdag, amat.type = "cpdag", x = 1, y = 2, set.type = "all")
adjDAll <- adjustment(amat = dag, amat.type = "dag", x = 1, y = 2, set.type = "all")
adjPAll <- adjustment(amat = dag, amat.type = "pdag", x = 1, y = 2, set.type = "all")

stopifnot( !identical(adjCAll, adjDAll), identical(adjDAll, adjPAll) )

adjCMin <- adjustment(amat = cpdag, amat.type = "cpdag", x = 1, y = 2, set.type = "minimal")
adjDMin <- adjustment(amat = dag, amat.type = "dag", x = 1, y = 2, set.type = "minimal")
adjPMin <- adjustment(amat = dag, amat.type = "pdag", x = 1, y = 2, set.type = "minimal")

stopifnot( !identical(adjCMin, adjDMin), identical(adjDMin, adjPMin) )


#####################################################################################
## Test 1: Compare CPDAG and PDAG implementation and validate all sets using gac()
#####################################################################################
nreps <- 100
simRes <- data.frame(setType = rep(NA, nreps), id = rep(NA,nreps), 
                     rtCPDAG = rep(NA,nreps), rtPDAG = rep(NA, nreps),
                     nSet = rep(NA, nreps), gacCheck = rep(NA, nreps))
proc.time()
for (i in 1:nreps) {
  cat("i = ",i,"\n")
  ## generate a graph
  seed <- i
  set.seed(seed)
  p <- sample(x=5:10, size = 1)
  prob <- sample(x=3:7/10, size = 1)
  g <- pcalg:::randomDAG(p, prob)  ## true DAG
  cpdag <- dag2cpdag(g)
  cpdag.mat <- t(as(cpdag,"matrix")) ## has correct encoding

  ## define input
  amat <- cpdag.mat
  x <- sample(x = 1:p, size = 1)
  y <- sample(x = setdiff(1:p,x), size = 1)
  set.type <- sample(x = c("all", "minimal"), size = 1)
  simRes$setType[i] <- set.type
  
  ## run both methods
  simRes$rtCPDAG[i] <- system.time(res1 <- adjustment(amat = amat, amat.type = "cpdag", x = x, y = y, set.type = set.type))[3]
  simRes$rtPDAG[i] <- system.time(res2 <- adjustment(amat = amat, amat.type = "pdag", x = x, y = y, set.type = set.type))[3]
  simRes$nSet[i] <- length(res1)
  
  if (length(res1) == 0) {
    res1 <- vector("list", 0)
  }
  if (length(res2) == 0) {
    res2 <- vector("list", 0)
  }
  ## compare results
  simRes$id[i] <- identical(res1,res2)
  
  ## compare results with gac() based on "pdag"
  if (length(res2) > 0) {
    gc <- TRUE
    for (j in 1:length(res2)) {
      gc <- gc & gac(amat = amat, x = x, y = y, z = res2[[j]], type = "cpdag")$gac
    }
    simRes$gacCheck[i] <-  gc
  }
  
}
proc.time()

summary(simRes)
table(is.na(simRes$gacCheck), simRes$nSet == 0)

################################################
## Test 2: Check using predefined graphs
################################################
gacVSadj <- function(amat, x, y ,z, V, type) {
  ## gac(z) is TRUE IFF z is returned by adjustment()
  ## x,y,z: col positions as used in GAC
  ## Result: TRUE is result is equal
  typeDG <- switch(type,
                   dag = "dag",
                   cpdag = "cpdag",
                   mag = "mag",
                   pag = "pag")
  gacRes <- gac(amat,x,y, z, type)$gac 
  adjRes <- adjustment(amat = amat, amat.type = typeDG, x = x, y = y, set.type = "all")
  if (gacRes) { ## z is valid adj set
    res <- any(sapply(adjRes, function(xx) setequal(z, xx)))
  } else { ## z is not valid adj set
    res <- all(!sapply(adjRes, function(xx) setequal(z, xx)))
  }
  res
}

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

xx <- xx &  gacVSadj(mFig1,x,y, z=c(2,4), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,5), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,2,1), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,5,1), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,2,5), V=V, type)
xx <- xx &  gacVSadj(mFig1,x,y, z=c(4,2,5,1), V=V, type)
xx <- xx & gacVSadj(mFig1,x,y, z= 2,    V=V, type)
xx <- xx & gacVSadj(mFig1,x,y, z= NULL, V=V, type)

## CPDAG 2: Paper Fig 5a
mFig5a <- matrix(c(0,1,0,0,0, 1,0,1,0,0, 0,0,0,0,1, 0,0,1,0,0, 0,0,0,0,0), 5,5)
V <- as.character(1:ncol(mFig5a))
rownames(mFig5a) <- colnames(mFig5a) <- V

type <- "cpdag"
x <- c(1,5); y <- 4
xx <- xx &  gacVSadj(mFig5a, x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSadj(mFig5a, x,y, z= 2,     V=V, type)

## DAG 1 from Marloes' Talk
mMMd1 <- matrix(c(0,1,0,1,0,0, 0,0,1,0,1,0, 0,0,0,0,0,1,
                  0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6)
V <- as.character(1:ncol(mMMd1))
rownames(mMMd1) <- colnames(mMMd1) <- V

type <- "dag"
x <- 1; y <- 3
xx <- xx &  gacVSadj(mMMd1, x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(mMMd1, x,y, z= 2, V=V, type)
xx <- xx &  gacVSadj(mMMd1, x,y, z= 4, V=V, type)
xx <- xx & gacVSadj(mMMd1, x,y, z= 5, V=V, type)
xx <- xx & gacVSadj(mMMd1, x,y, z= 6, V=V, type)
xx <- xx & gacVSadj(mMMd1, x,y, z=c(4,5), V=V, type)

## DAG 2 from Marloes' Talk
mMMd2 <- matrix(c(0,1,0,1,0,0, 0,0,0,0,0,0, 0,1,0,0,1,0,
                  0,0,0,0,1,0, 0,0,0,0,0,1, 0,0,0,0,0,0), 6,6)
V <- as.character(1:ncol(mMMd2))
rownames(mMMd2) <- colnames(mMMd2) <- V

type <- "dag"
x <- 4; y <- 6
xx <- xx &  gacVSadj(mMMd2, x,y, z= 1, V=V, type)
xx <- xx &  gacVSadj(mMMd2, x,y, z= 3, V=V, type)
xx <- xx & gacVSadj(mMMd2, x,y, z= 5, V=V, type)
xx <- xx & gacVSadj(mMMd2, x,y, z=c(1,5), V=V, type)
xx <- xx &  gacVSadj(mMMd2, x,y, z=c(1,2), V=V, type)
xx <- xx &  gacVSadj(mMMd2, x,y, z=c(1,3), V=V, type)
xx <- xx & gacVSadj(mMMd2, x,y, z= 2, V=V, type)

##################################################
## PAG
##################################################
type <- "pag"
mFig3a <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,1, 0,1,1,0), 4,4)
V <- as.character(1:ncol(mFig3a))
rownames(mFig3a) <- colnames(mFig3a) <- V
xx <- xx & gacVSadj(mFig3a, x=2,      y=4, z=NULL,   V=V, type)

mFig3b <- matrix(c(0,2,0,0, 3,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
V <- as.character(1:ncol(mFig3b))
rownames(mFig3b) <- colnames(mFig3b) <- V
xx <- xx & gacVSadj(mFig3b, x=2,      y=4, z=NULL,   V=V, type)

mFig3c <- matrix(c(0,3,0,0, 2,0,3,3, 0,2,0,3, 0,2,2,0), 4,4)
V <- as.character(1:ncol(mFig3c))
rownames(mFig3c) <- colnames(mFig3c) <- V
xx <- xx &  gacVSadj(mFig3c, x=2,      y=4, z=NULL,   V=V, type)

mFig4a <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,3,3,2,
                   0,0,2,0,2,2, 0,0,2,1,0,2, 0,0,1,3,3,0), 6,6)
V <- as.character(1:ncol(mFig4a))
rownames(mFig4a) <- colnames(mFig4a) <- V
xx <- xx & gacVSadj(mFig4a, x=3,      y=4, z=NULL,   V=V, type)
xx <- xx &  gacVSadj(mFig4a, x=3,      y=4, z= 6,     V=V, type)
xx <- xx &  gacVSadj(mFig4a, x=3,      y=4, z=c(1,6), V=V, type)
xx <- xx &  gacVSadj(mFig4a, x=3,      y=4, z=c(2,6), V=V, type)
xx <- xx &  gacVSadj(mFig4a, x=3,      y=4, z=c(1,2,6), V=V, type)

mFig4b <- matrix(c(0,0,1,0,0,0, 0,0,1,0,0,0, 2,2,0,0,3,2,
                   0,0,0,0,2,2, 0,0,2,3,0,2, 0,0,2,3,2,0), 6,6)
V <- as.character(1:ncol(mFig4b))
rownames(mFig4b) <- colnames(mFig4b) <- V
xx <- xx & gacVSadj(mFig4b, x=3,      y=4, z=NULL,   V=V, type)
xx <- xx & gacVSadj(mFig4b, x=3,      y=4, z= 6,     V=V, type)
xx <- xx & gacVSadj(mFig4b, x=3,      y=4, z=c(5,6), V=V, type)

mFig5b <- matrix(c(0,1,0,0,0,0,0, 2,0,2,3,0,3,0, 0,1,0,0,0,0,0, 0,2,0,0,3,0,0,
                   0,0,0,2,0,2,3, 0,2,0,0,2,0,0, 0,0,0,0,2,0,0), 7,7)
V <- as.character(1:ncol(mFig5b))
rownames(mFig5b) <- colnames(mFig5b) <- V
xx <- xx & gacVSadj(mFig5b, x=c(2,7), y=6, z=NULL,   V=V, type)
xx <- xx &  gacVSadj(mFig5b, x=c(2,7), y=6, z=c(4,5), V=V, type)
xx <- xx &  gacVSadj(mFig5b, x=c(2,7), y=6, z=c(4,5,1), V=V, type)
xx <- xx &  gacVSadj(mFig5b, x=c(2,7), y=6, z=c(4,5,3), V=V, type)
xx <- xx &  gacVSadj(mFig5b, x=c(2,7), y=6, z=c(1,3,4,5), V=V, type)

## PAG in Marloes' talk
mMMp <- matrix(c(0,0,0,3,2,0,0, 0,0,0,0,1,0,0, 0,0,0,0,1,0,0, 2,0,0,0,0,3,2,
                 3,2,2,0,0,0,3, 0,0,0,2,0,0,0, 0,0,0,2,2,0,0), 7,7)
V <- as.character(1:ncol(mMMp))
rownames(mMMp) <- colnames(mMMp) <- V

x <- c(5,6); y <- 7
xx <- xx & gacVSadj(mMMp, x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(mMMp, x,y, z= 1,   V=V, type)
xx <- xx & gacVSadj(mMMp, x,y, z= 4,   V=V, type)
xx <- xx & gacVSadj(mMMp, x,y, z= 2,   V=V, type)
xx <- xx & gacVSadj(mMMp, x,y, z= 3,   V=V, type)
xx <- xx & gacVSadj(mMMp, x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSadj(mMMp, x,y, z=c(1,4), V=V, type)
xx <- xx &  gacVSadj(mMMp, x,y, z=c(1,4,2), V=V, type)
xx <- xx &  gacVSadj(mMMp, x,y, z=c(1,4,3), V=V, type)
xx <- xx &  gacVSadj(mMMp, x,y, z=c(1,4,2,3), V=V, type)

##################################################
## V=V, type = "pag" -- Tests from Ema
##################################################
type <- "pag"
pag.m <- readRDS(system.file(package="pcalg", "external", "gac-pags.rds"))
m1 <- pag.m[["m1"]]
V <- colnames(m1)
x <- 6; y <- 9
xx <- xx & gacVSadj(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=1, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=3, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=4, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,7,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,5,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,5,7,8), V=V, type)

x <- c(6,8); y <- 9
xx <- xx & gacVSadj(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=1, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=3, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=4, V=V, type)
xx <- xx &  gacVSadj(m1,x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSadj(m1,x,y, z=c(2,3,4), V=V, type)
xx <- xx &  gacVSadj(m1,x,y, z=c(2,3,7), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,5), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,3,5,7), V=V, type)

x <- 3; y <- 1
xx <- xx & gacVSadj(m1,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=4, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=5, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=6, V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,6), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,7,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,5,8), V=V, type)
xx <- xx & gacVSadj(m1,x,y, z=c(2,5,7,8), V=V, type)

m2 <- pag.m[["m2"]]
V <- colnames(m2)
x <- 3; y <-1
xx <- xx & gacVSadj(m2,x,y, z=NULL, V=V, type)
xx <- xx &  gacVSadj(m2,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=4, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=8, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=9, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=c(2,8,9), V=V, type)
xx <- xx &  gacVSadj(m2,x,y, z=c(2,5), V=V, type)

x <- c(3,9); y <- 1
xx <- xx & gacVSadj(m2,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=4, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=c(2,8), V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=8, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=9, V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=c(2,8,9), V=V, type)
xx <- xx & gacVSadj(m2,x,y, z=c(2,5), V=V, type)

m3 <- pag.m[["m3"]]
V <- colnames(m3)
x <- 1; y <- 9
xx <- xx & gacVSadj(m3,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=3, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=5, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=7, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=8, V=V, type)
xx <- xx &  gacVSadj(m3,x,y, z=c(2,3), V=V, type)
xx <- xx &  gacVSadj(m3,x,y, z=c(5,7), V=V, type)

x <- 1; y <- 8
xx <- xx & gacVSadj(m3,x,y, z=NULL, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=2, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=3, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=5, V=V, type)
xx <- xx &  gacVSadj(m3,x,y, z=7, V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=9, V=V, type)
xx <- xx &  gacVSadj(m3,x,y, z=c(2,3), V=V, type)
xx <- xx & gacVSadj(m3,x,y, z=c(5,9), V=V, type)

if (!xx) {
  stop("Problem when testing function gacVSadj.")
} else {
  message("OK, no issues were found.")
}

##################################################
## given same graph, type=cpdag and type=pdag
## should give same canonical set
##################################################
m <- rbind(c(0,1,0,0,0,0),
           c(1,0,1,0,0,0),
           c(0,1,0,0,0,0),
           c(0,0,0,0,0,0),
           c(0,1,1,1,0,0),
           c(1,0,1,1,1,0))
colnames(m) <- rownames(m) <- as.character(1:6)

## You can see that the current adjustment function outputs different sets
## if type = "cpdag" or type = "pdag" which shouldn't happen 
## because it is the same graph:
res1 <- adjustment(m,amat.type="cpdag",2,4,set.type="canonical")
res2 <- adjustment(m,amat.type="pdag",2,4,set.type="canonical")

if (!all.equal(res1, res2)) {
    stop("Canonical set is not the same for type=cpdag and type=pdag\n")
}

}
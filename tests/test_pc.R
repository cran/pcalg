library(pcalg)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})


##################################################
## Standard PC
##################################################
nreps <- 10

set.seed(234)
for (ii in 1:nreps) {
  p <- 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.3)
  myCPDAG <- dag2cpdag(myDAG)
  suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
  res <- pc(suffStat, indepTest=gaussCItest, p, 0.99)
  if( shd(res, myCPDAG) != 0)
    stop("Test pc wrong: CPDAG ",ii," was not found correctly")
}

showProc.time()

##################################################
## Conservative PC
##################################################
##PC algorithm sample (compared with Tetrad)
##__________________________________________________________________________________
## Example 1
p <- 10
n <- 5000
set.seed(p*37+15673)
g <- randomDAG(p,2/(p-1))

## generate n samples of DAG using standard normal error distribution
##random data

load(system.file("external", "test_conservative_pc_data1.rda", package = "pcalg"))
## set.seed(67*37)
## new.mat <- rmvnorm(n,mean=rep(0,p),sigma=trueCov(g))
## save(new.mat, file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data1.rda")
## load(file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data1.rda")

suffStat.data <- list(C=cor(new.mat1),n=n)
indepTest.data <- gaussCItest

##pcAlgo conservative sample
dag1 <- pc(suffStat.data, indepTest.data, p, alpha=0.005, u2pd="relaxed",conservative=TRUE)

##adjacency matrix
dag1.amat <- as(dag1@graph,"matrix")

##always save the transpose of the matrix to be loaded in Tetrad
##write(t(new.mat1),file="test_conservative_pc_data1.txt",ncolumns=p)

##check the output with Tetrad

amat.tetrad1 <- matrix(0,p,p)
amat.tetrad1[1,6] <- 1
amat.tetrad1[2,4] <- amat.tetrad1[2,6] <- 1
amat.tetrad1[3,4] <- 1
amat.tetrad1[4,6] <- amat.tetrad1[4,9] <- 1

correctEst1 <- all(dag1.amat == amat.tetrad1)
if (!correctEst1) stop("Test sample conservative PC wrong: example 1!")
showProc.time()


## Example 2
p <- 12
n <- 5000
set.seed(p*67+15673)
g <- randomDAG(p,2/(p-1))

## generate n samples of DAG using standard normal error distribution
##random data
load(system.file("external", "test_conservative_pc_data2.rda", package = "pcalg"))
## set.seed(67*67)
## new.mat <- rmvnorm(n,mean=rep(0,p),sigma=trueCov(g))
## save(new.mat, file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data2.rda")
## load(file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data2.rda")

suffStat.data <- list(C=cor(new.mat2),n=n)
indepTest.data <- gaussCItest

##pcAlgo conservative sample
dag2 <- pc(suffStat.data, indepTest.data, p, alpha=0.005, verbose=FALSE, u2pd="relaxed", conservative=TRUE)

##adjacency matrix
dag2.amat <- as(dag2@graph,"matrix")

##always save the transpose of the matrix to be loaded in Tetrad
##write(t(new.mat),file="test_conservative_pc_data2.txt",ncolumns=p)

##check the output with Tetrad

amat.tetrad2 <- matrix(0,p,p)
amat.tetrad2[1,2] <- amat.tetrad2[1,7] <- 1
amat.tetrad2[2,1] <- amat.tetrad2[2,11] <- 1
amat.tetrad2[3,9] <- amat.tetrad2[3,11] <- amat.tetrad2[3,12] <- 1
amat.tetrad2[4,5] <- amat.tetrad2[4,11] <- amat.tetrad2[4,12] <- 1
amat.tetrad2[5,4] <- amat.tetrad2[5,8] <- amat.tetrad2[5,9] <- 1
amat.tetrad2[7,1] <- 1
amat.tetrad2[8,5] <- amat.tetrad2[8,12] <- 1
amat.tetrad2[10,11] <- 1

correctEst2 <- all(dag2.amat == amat.tetrad2)
if (!correctEst2) stop("Test sample conservative PC wrong: example 2!")
showProc.time()


## Example 3
## new.mat <- read.table("/u/colombo/Diss/RAusw/AFCI/Conservative_algorithms/test_conservative_pc_data3.txt")
## save(new.mat, file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data3.rda")
## load(file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data3.rda")
load(system.file("external", "test_conservative_pc_data3.rda", package = "pcalg"))

##always save the transpose of the matrix to be loaded in Tetrad
##write(t(new.mat),file="test_conservative_pc_data3.txt",ncolumns=15)
suffStat.data <- list(C=cor(new.mat3),n=n)
indepTest.data <- gaussCItest

##pcAlgo conservative sample
dag3 <- pc(suffStat.data, indepTest.data, 15, alpha=0.005, verbose=FALSE, u2pd="relaxed", conservative=TRUE)

##adjacency matrix
dag3.amat <- as(dag3@graph,"matrix")

##check the output with Tetrad
amat.tetrad3 <- matrix(0,15,15)
amat.tetrad3[1,2] <- amat.tetrad3[1,3] <- amat.tetrad3[1,10] <- 1
amat.tetrad3[2,1] <- amat.tetrad3[2,5] <- amat.tetrad3[2,6] <- amat.tetrad3[2,12] <- 1
amat.tetrad3[3,1] <- amat.tetrad3[3,9] <- 1
amat.tetrad3[4,5] <- amat.tetrad3[4,7] <- 1
amat.tetrad3[5,10] <- 1
amat.tetrad3[6,2] <- amat.tetrad3[6,9] <- amat.tetrad3[6,11] <- 1
amat.tetrad3[7,4] <- amat.tetrad3[7,8] <- 1
amat.tetrad3[8,15] <- 1
amat.tetrad3[9,12] <- amat.tetrad3[9,14] <- 1
amat.tetrad3[10,1] <- amat.tetrad3[10,5] <- amat.tetrad3[10,12] <- amat.tetrad3[10,13] <- 1
amat.tetrad3[11,6] <- amat.tetrad3[11,13] <- 1
amat.tetrad3[12,2] <- amat.tetrad3[12,15] <- 1
amat.tetrad3[13,15] <- 1
amat.tetrad3[14,7] <- 1
amat.tetrad3[15,14] <- 1

correctEst3 <- all(dag3.amat == amat.tetrad3)
if (!correctEst3) stop("Test sample conservative PC wrong: example 3!")
showProc.time()


##PC algorithm population
##_____________________________________________________________________________
## Example 4
p <- 15
set.seed(15673)
g <- randomDAG(p,2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
indepTest <- gaussCItest
dag4 <- pc(suffStat, indepTest, p, alpha=0.9999, verbose=FALSE, u2pd="relaxed")
dag5 <- pc(suffStat, indepTest, p, alpha=0.9999, verbose=FALSE, u2pd="relaxed",conservative=TRUE)

##adjacency matrix
dag4.amat <- as(dag4@graph,"matrix")
dag5.amat <- as(dag5@graph,"matrix")

correctEst4 <- all(dag4.amat == dag5.amat)
if (!correctEst4) stop("Test population conservative PC wrong: example 4!")
showProc.time()


## Example 5
p <- 25
set.seed(1589873)
g <- randomDAG(p,2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
indepTest <- gaussCItest
dag6 <- pc(suffStat, indepTest, p, alpha=0.9999, verbose=FALSE, u2pd="relaxed")
dag7 <- pc(suffStat, indepTest, p, alpha=0.9999, verbose=FALSE, u2pd="relaxed",conservative=TRUE)

##adjacency matrix
dag6.amat <- as(dag6@graph,"matrix")
dag7.amat <- as(dag7@graph,"matrix")

correctEst5 <- all(dag6.amat == dag7.amat)
if (!correctEst5) stop("Test population conservative PC wrong: example 5!")
showProc.time()

## Example 6
p <- 35
set.seed(78673)
g <- randomDAG(p,2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
indepTest <- gaussCItest
dag8 <- pc(suffStat, indepTest, p, alpha=0.9999, u2pd="relaxed")
dag9 <- pc(suffStat, indepTest, p, alpha=0.9999, u2pd="relaxed", conservative=TRUE)

##adjacency matrix
dag8.amat <- as(dag8@graph,"matrix")
dag9.amat <- as(dag9@graph,"matrix")

correctEst6 <- all(dag8.amat == dag9.amat)
if (!correctEst6) stop("Test population conservative PC wrong: example 6!")

showProc.time()

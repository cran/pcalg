library(pcalg)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

########################################################
##
##       Example 1: Zhang (2008), Fig. 6, p.1882
##                  Paper with rules
##
########################################################


## create the graph g
p <- 4
amat1 <- rbind(c(0,1,0,0,1),
               c(0,0,1,0,0),
               c(0,0,0,1,0),
               c(0,0,0,0,0),
               c(0,0,0,1,0))
colnames(amat1) <- rownames(amat1) <- as.character(1:5)
L1 <- 1
V1 <- as.character(1:5)
edL1 <- vector("list",length=5)
names(edL1) <- V1
edL1[[1]] <- list(edges=c(2,4),weights=c(1,1))
edL1[[2]] <- list(edges=3,weights=c(1))
edL1[[3]] <- list(edges=5,weights=c(1))
edL1[[4]] <- list(edges=5,weights=c(1))
g1 <- new("graphNEL", nodes=V1, edgeL=edL1,edgemode="directed")

## compute the true covariance matrix of g1
cov.mat1 <- trueCov(g1)

## delete rows and columns which belong to L1
true.cov1 <- cov.mat1[-L1,-L1]

## transform it into a correlation matrix
true.corr1 <- cov2cor(true.cov1)

##PAG
suffStat1 <- list(C = true.corr1, n = 10^9)
indepTest1 <- gaussCItest

true.pag1 <- fci(suffStat1, indepTest1, p, alpha = 0.99, verbose=FALSE)

##define correct PAG
corr.pag1 <- rbind(c(0,1,1,0),
                   c(1,0,0,2),
                   c(1,0,0,2),
                   c(0,3,3,0))
correctEst1 <- all(corr.pag1 == true.pag1@amat)
if (!correctEst1) stop("Test fci wrong: example 1!")
showProc.time()

########################################################
##
##       Example 2: Zhang (2006), Fig. 5.2, p.198
##                  Dissertation
##
########################################################

## create the graph g
p <- 5
amat2 <- rbind(c(0,0,0,1,0,0),
               c(0,0,0,1,1,0),
               c(0,0,0,1,0,0),
               c(0,0,0,0,1,0),
               c(0,0,0,0,0,0),
               c(0,0,1,1,0,0))
colnames(amat2) <- rownames(amat2) <- as.character(1:6)
L2 <- 6
V2 <- as.character(1:6)
edL2 <- vector("list",length=6)
names(edL2) <- V2
edL2[[1]] <- list(edges=4,weights=1)
edL2[[2]] <- list(edges=c(4,5),weights=c(0.5,1))
edL2[[3]] <- list(edges=4,weights=1)
edL2[[4]] <- list(edges=5,weights=1)
edL2[[6]] <- list(edges=c(3,4),weights=c(1,1))
g2 <- new("graphNEL", nodes=V2, edgeL=edL2,edgemode="directed")

## compute the true covariance matrix of g
cov.mat2 <- trueCov(g2)

## delete rows and columns which belong to L
true.cov2 <- cov.mat2[-L2,-L2]

## transform it into a correlation matrix
true.corr2 <- cov2cor(true.cov2)

##PAG
suffStat2 <- list(C = true.corr2, n = 10^9)
indepTest2 <- gaussCItest

true.pag2 <- fci(suffStat2, indepTest2, p, alpha = 0.99, verbose=FALSE)

##define correct PAG
corr.pag2 <- rbind(c(0,0,0,2,0),
                   c(0,0,0,2,2),
                   c(0,0,0,2,0),
                   c(1,1,1,0,2),
                   c(0,3,0,3,0))

correctEst2 <- all(corr.pag2 == true.pag2@amat)
if (!correctEst2) stop("Test fci wrong: example 2!")
showProc.time()



########################################################
##
##             Example 3: random DAG
##
########################################################


set.seed(40)
##Random graph only R1-R10
g3 <- randomDAG(14,0.3)

##Define the latente variables
L3 <- c(8,10)

##pcAlgo.Perfect with true correlation matrix
##______________________________________________________
p <- 12
amat.g <- as(g3,"matrix")
colnames(amat.g) <- rownames(amat.g) <- graph::nodes(g3)
amat.g[amat.g!=0] <- 1

##Compute the true covariance matrix of g
cov.mat3 <- trueCov(g3)

##Delete rows and columns which belong to L
true.cov3 <- cov.mat3[-L3,-L3]

##Transform it in a correlation matrix
true.corr3 <- cov2cor(true.cov3)

##PAG
suffStat3 <- list(C = true.corr3, n = 10^9)
indepTest3 <- gaussCItest

true.pag3 <- fci(suffStat3, indepTest3, p, alpha = 0.99, verbose=FALSE)

##define correct PAG
corr.pag3 <- rbind(c(0,0,2,0,0,2,0,0,0,2,2,2),
                   c(0,0,2,0,2,0,0,0,0,2,2,2),
                   c(1,1,0,0,2,0,2,2,0,0,0,0),
                   c(0,0,0,0,0,2,0,0,2,0,2,2),
                   c(0,3,3,0,0,2,2,0,2,2,2,0),
                   c(3,0,0,1,3,0,2,0,2,0,0,0),
                   c(0,0,3,0,3,3,0,0,0,0,0,0),
                   c(0,0,3,0,0,0,0,0,0,0,0,0),
                   c(0,0,0,3,3,3,0,0,0,0,0,0),
                   c(3,3,0,0,3,0,0,0,0,0,2,0),
                   c(3,3,0,1,3,0,0,0,0,1,0,2),
                   c(3,3,0,3,0,0,0,0,0,0,3,0))

correctEst3 <- all(corr.pag3 == true.pag3@amat)
if (!correctEst3) stop("Test fci wrong: example 3!")
showProc.time()


#########################################################################################
##
##      Example 4: Spirtes 1997 p.21 DAG with latent variables and p.24 PAG
##
#########################################################################################

p <- 5
amat4 <- t(matrix(c(0,0,0,0,1,1,0, 0,0,0,1,0,0,1, 0,0,0,1,0,1,0, 0,0,0,0,1,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,1, 0,0,0,0,0,0,0),7,7))
colnames(amat4) <- rownames(amat4) <- as.character(1:7)
L4 <- c(1,2)
V4 <- as.character(1:7)
edL4 <- vector("list",length=7)
names(edL4) <- V4
edL4[[1]] <- list(edges=c(5,6),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[2]] <- list(edges=c(4,7),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[3]] <- list(edges=c(4,6),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL4[[4]] <- list(edges=5,weights=c(abs(rnorm(1))))
edL4[[6]] <- list(edges=7,weights=c(abs(rnorm(1))))
g4 <- new("graphNEL", nodes=V4, edgeL=edL4,edgemode="directed")

## compute the true covariance matrix of g1
cov.mat4 <- trueCov(g4)

## delete rows and columns which belong to L1
true.cov4 <- cov.mat4[-L4,-L4]

## transform it into a correlation matrix
true.corr4 <- cov2cor(true.cov4)

##PAG
suffStat4 <- list(C = true.corr4, n = 10^9)
indepTest4 <- gaussCItest

true.pag4 <- fci(suffStat4, indepTest4, p, alpha = 0.99, verbose=FALSE)

##define correct PAG
corr.pag4 <- rbind(c(0,2,0,2,0),
                   c(1,0,2,0,2),
                   c(0,3,0,2,0),
                   c(1,0,2,0,2),
                   c(0,2,0,3,0))

correctEst4 <- all(corr.pag4 == true.pag4@amat)
if (!correctEst4) stop("Test fci wrong: example 4!")
showProc.time()

##################################################
## Conservative FCI
##################################################
##FCI algorithm population (using biCC=TRUE to reduce running time)
##_____________________________________________________________________________
## Example 7
p <- 15
set.seed(156)
g <- randomDAG(p,2.5/(p-1))
amat <- as(g,"matrix")
amat[amat!=0] <- 1

##find root nodes
tmp <- rep(NA,p)
for (j in 1:length(tmp)) {
  tmp[j] <- length(which(amat[,j]!=0))
}
poss.rnodes <- which(tmp==0)

new.tmp <- rep(0,length(poss.rnodes))
for (k in 1:length(poss.rnodes)) {
  for (ii in 1:(p)) {
    if (amat[poss.rnodes[k],ii]!=0) {
      new.tmp[k] <- new.tmp[k]+1
    }
  }
}
rnodes <- poss.rnodes[which(new.tmp>=2)]

##define latent variables
numb <- length(rnodes)/2
set.seed(234543+numb)
L <- sample(rnodes,numb)

##Compute the true covariance matrix of g
cov.mat <- trueCov(g)

##Delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]

##Transform it in a correlation matrix
true.corr <- cov2cor(true.cov)

suffStat <- list(C=true.corr,n=10^9)
indepTest <- gaussCItest

##population version
pag1 <- fci(suffStat, indepTest, dim(true.corr)[1], alpha=0.9999, verbose=FALSE, biCC=TRUE)@amat
pag2 <- fci(suffStat, indepTest, dim(true.corr)[1], alpha=0.9999, verbose=FALSE, conservative=c(TRUE,TRUE), biCC=TRUE, cons.rules=TRUE)@amat

correctEst7 <- all(pag1 == pag2)
if (!correctEst7) stop("Test population conservative FCI wrong: example 7!")
showProc.time()

## Left out example 8 because of long runtime

## Example 9
p <- 20
set.seed(9642)
g <- randomDAG(p,2.5/(p-1))
amat <- as(g,"matrix")
amat[amat!=0] <- 1

##find root nodes
tmp <- rep(NA,p)
for (j in 1:length(tmp)) {
  tmp[j] <- length(which(amat[,j]!=0))
}
poss.rnodes <- which(tmp==0)

new.tmp <- rep(0,length(poss.rnodes))
for (k in 1:length(poss.rnodes)) {
  for (ii in 1:(p)) {
    if (amat[poss.rnodes[k],ii]!=0) {
      new.tmp[k] <- new.tmp[k]+1
    }
  }
}
rnodes <- poss.rnodes[which(new.tmp>=2)]

##define latent variables
numb <- ceiling(length(rnodes)/2)
set.seed(234543+numb)
L <- sample(rnodes,numb)

##Compute the true covariance matrix of g
cov.mat <- trueCov(g)

##Delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]

##Transform it in a correlation matrix
true.corr <- cov2cor(true.cov)

suffStat <- list(C=true.corr,n=10^9)
indepTest <- gaussCItest

##population version
pag5 <- fci(suffStat, indepTest, dim(true.corr)[1], alpha=0.9999, verbose=FALSE, biCC=TRUE)@amat
pag6 <- fci(suffStat, indepTest, dim(true.corr)[1], alpha=0.9999, verbose=FALSE, conservative=c(TRUE,TRUE), biCC=TRUE, cons.rules=TRUE)@amat

correctEst9 <- all(pag5 == pag6)
if (!correctEst9) stop("Test population conservative FCI wrong: example 9!")
showProc.time()

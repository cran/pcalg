library(pcalg)


########################################################
##
##       Example 1: Zhang (2008), Fig. 6, p.1882
##                  Paper with rules
##
########################################################


## create the graph g
amat1 <- t(matrix(c(0,1,0,0,1, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,0, 0,0,0,1,0),5,5))
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

## compute the true CPDAG
true.CPDAG1 <- pcAlgo.Perfect(true.corr1, directed=FALSE, psepset=TRUE, verbose=0)

## orient it with the AFCI algorithm
rules <- rep(TRUE,10)
true.pag1 <- udag2pag(true.CPDAG1, rules=rules, verbose=0)

##define correct PAG
corr.pag1 <- matrix(0,4,4)
corr.pag1[1,] <- c(0,1,1,0)
corr.pag1[2,] <- c(1,0,0,2)
corr.pag1[3,] <- c(1,0,0,2)
corr.pag1[4,] <- c(0,3,3,0)

correctEst1 <- FALSE
if (all(corr.pag1==true.pag1)){
  correctEst1 <- TRUE
}

if (!correctEst1) stop("Test udag2pag wrong: example 1!")


########################################################
##
##       Example 2: Zhang (2006), Fig. 5.2, p.198
##                  Dissertation
##
########################################################

## create the graph g
amat2 <- t(matrix(c(0,0,0,1,0,0, 0,0,0,1,1,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,0, 0,0,1,1,0,0),6,6))
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

## compute the true CPDAG
true.CPDAG2 <- pcAlgo.Perfect(true.corr2, directed=FALSE, psepset=TRUE, verbose=0)

## orient it with the AFCI algorithm
rules <- rep(TRUE,10)
true.pag2 <- udag2pag(true.CPDAG2, rules=rules, verbose=0)

##define correct PAG
corr.pag2 <- matrix(0,5,5)
corr.pag2[1,] <- c(0,0,0,2,0)
corr.pag2[2,] <- c(0,0,0,2,2)
corr.pag2[3,] <- c(0,0,0,2,0)
corr.pag2[4,] <- c(1,1,1,0,2)
corr.pag2[5,] <- c(0,3,0,3,0)

correctEst2 <- FALSE
if (all(corr.pag2==true.pag2)){
  correctEst2 <- TRUE
}

if (!correctEst2) stop("Test udag2pag wrong: example 2!")



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

amat.g <- as(g3,"matrix")
colnames(amat.g) <- rownames(amat.g) <- nodes(g3)
amat.g[amat.g!=0] <- 1

##Compute the true covariance matrix of g
cov.mat3 <- trueCov(g3)

##Delete rows and columns which belong to L
true.cov3 <- cov.mat3[-L3,-L3]

##Transform it in a correlation matrix
true.corr3 <- cov2cor(true.cov3)

##Compute the true CPDAG
true.CPDAG3 <- pcAlgo.Perfect(true.corr3, directed=FALSE, psepset=TRUE, verbose=0)

##Orient it with the FCI algorithm
rules <- rep(TRUE,10)
true.pag3 <- udag2pag(true.CPDAG3, rules=rules, verbose=0)
       
##define correct PAG
corr.pag3 <- matrix(0,12,12)
corr.pag3[1,] <- c(0,0,2,0,0,2,0,0,0,2,2,2)
corr.pag3[2,] <- c(0,0,2,0,2,0,0,0,0,2,2,2)
corr.pag3[3,] <- c(1,1,0,0,2,0,2,2,0,0,0,0)
corr.pag3[4,] <- c(0,0,0,0,0,2,0,0,2,0,2,2)
corr.pag3[5,] <- c(0,3,3,0,0,2,2,0,2,2,2,0)
corr.pag3[6,] <- c(3,0,0,1,3,0,2,0,2,0,0,0)
corr.pag3[7,] <- c(0,0,3,0,3,3,0,0,0,0,0,0)
corr.pag3[8,] <- c(0,0,3,0,0,0,0,0,0,0,0,0)
corr.pag3[9,] <- c(0,0,0,3,3,3,0,0,0,0,0,0)
corr.pag3[10,] <- c(3,3,0,0,3,0,0,0,0,0,2,0)
corr.pag3[11,] <- c(3,3,0,1,3,0,0,0,0,1,0,2)
corr.pag3[12,] <- c(3,3,0,3,0,0,0,0,0,0,3,0)


correctEst3 <- FALSE
if (all(corr.pag3==true.pag3)){
  correctEst3 <- TRUE
}

if (!correctEst3) stop("Test udag2pag wrong: example 3!")


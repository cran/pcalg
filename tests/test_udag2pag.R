library(pcalg)

########################################################
##
##       Example 1: Zhang (2008), Fig. 6, p.1882
##                  Paper with rules
##
########################################################

. <- 0 ## MM's trick to visualize sparse matrices on *input*

## create the graph g
amat1 <- rbind(c(.,1,.,.,1),
               c(.,.,1,.,.),
               c(.,.,.,1,.),
               c(.,.,.,.,.),
               c(.,.,.,1,.))
colnames(amat1) <- rownames(amat1) <- as.character(1:5)
L1 <- 1
V1 <- as.character(1:5)
edL1 <- list("1" = list(edges=c(2,4), weights=c(1,1)),
             "2" = list(edges= 3, weights=c(1)),
             "3" = list(edges= 5, weights=c(1)),
             "4" = list(edges= 5, weights=c(1)), "5" = NULL)
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
corr.pag1 <- rbind(c(.,1,1,.),
                   c(1,.,.,2),
                   c(1,.,.,2),
                   c(.,3,3,.))
correctEst1 <- all(corr.pag1==true.pag1)
if (!correctEst1) stop("Test udag2pag wrong: example 1!")

cat('Time elapsed: ', (.pt <- proc.time()),'\n') # "stats"

########################################################
##
##       Example 2: Zhang (2006), Fig. 5.2, p.198
##                  Dissertation
##
########################################################

## create the graph g
amat2 <- rbind(c(.,.,.,1,.,.),
               c(.,.,.,1,1,.),
               c(.,.,.,1,.,.),
               c(.,.,.,.,1,.),
               c(.,.,.,.,.,.),
               c(.,.,1,1,.,.))
colnames(amat2) <- rownames(amat2) <- as.character(1:6)
L2 <- 6
edL2 <- list(V1 = list(edges=4,weights=1),
             V2 = list(edges=c(4,5),weights=c(0.5,1)),
             V3 = list(edges=4,weights=1),
             V4 = list(edges=5,weights=1),
             V5 = NULL,
             V6 = list(edges=c(3,4),weights=c(1,1)))
g2 <- new("graphNEL", nodes = names(edL2), edgeL=edL2, edgemode="directed")

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
corr.pag2 <- rbind(c(.,.,.,2,.),
                   c(.,.,.,2,2),
                   c(.,.,.,2,.),
                   c(1,1,1,.,2),
                   c(.,3,.,3,.))
correctEst2 <- all(corr.pag2==true.pag2)
if (!correctEst2) stop("Test udag2pag wrong: example 2!")

cat('Time elapsed: ', proc.time() - .pt,'\n') # "stats"


########################################################
##
##             Example 3: random DAG
##
########################################################


set.seed(40)
##Random graph only R1-R10
g3 <- randomDAG(14,0.3)

##Define the latent variables
L3 <- c(8,10)

##pcAlgo.Perfect with true correlation matrix
##______________________________________________________

amat.g <- as(g3,"matrix")
colnames(amat.g) <- rownames(amat.g) <- graph::nodes(g3)
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
corr.pag3 <- rbind(c(.,.,2,.,.,2,.,.,.,2,2,2),
                   c(.,.,2,.,2,.,.,.,.,2,2,2),
                   c(1,1,.,.,2,.,2,2,.,.,.,.),
                   c(.,.,.,.,.,2,.,.,2,.,2,2),
                   c(.,3,3,.,.,2,2,.,2,2,2,.),
                   c(3,.,.,1,3,.,2,.,2,.,.,.),
                   c(.,.,3,.,3,3,.,.,.,.,.,.),
                   c(.,.,3,.,.,.,.,.,.,.,.,.),
                   c(.,.,.,3,3,3,.,.,.,.,.,.),
                   c(3,3,.,.,3,.,.,.,.,.,2,.),
                   c(3,3,.,1,3,.,.,.,.,1,.,2),
                   c(3,3,.,3,.,.,.,.,.,.,3,.))

correctEst3 <- all(corr.pag3==true.pag3)
if (!correctEst3) stop("Test udag2pag wrong: example 3!")

cat('Time elapsed: ', proc.time() - .pt,'\n') # "stats"

if(!interactive()) warnings()

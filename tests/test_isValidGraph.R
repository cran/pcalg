library(pcalg)
possible_p <- c(seq(10,30,by=10))
possible_neighb_size <- c(3:10)

num_setings <-10
num_rep <- 2

##amount of bg knowledge in percent
p_bg <- seq(0,1,by=0.1)

counterr <- 0

seed1 <- 1 
for (r in seed1:num_setings) {
  set.seed(r)
  ##Then we can choose one setting:
  p <-  sample(possible_p,1)
  neigh <- sample(possible_neighb_size,1)
  prob <- round(neigh/(p-1),3)

  ## cat("r=",r,", p=",p,", nb=",neigh,"\n")     
  
  ## then for every setting selected we can generate i.e. 20 graph
  for (d in 1: num_rep){
    ##i left the seed inside for now
    ## i am not sure which seed option makes it easier for us
    ## feel free to move this around
    #seed <- sample(possible_seed,1)
    #set.seed(seed)
    
    ##get graph
    g <- pcalg:::randomDAG(p, prob)  ## true DAG
    cpdag <- dag2cpdag(g) ## true CPDAG
    
    ## get adjacency matrix
    cpdag.amat <- t(as(cpdag,"matrix"))
    dag.amat <- t(as(g,"matrix"))
    dag.amat[dag.amat != 0] <- 1

    amat.cpdag <- cpdag.amat
    if (!isValidGraph(amat.cpdag,type="cpdag")){
      counterr <- counterr +1
    }
    if (!isValidGraph(dag.amat,type="dag")){
      counterr <- counterr +1
    }
  }
}
counterr

amat <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## graph::plot(as(t(amat), "graphNEL"))             
stopifnot(isValidGraph(amat = amat, type = "dag"), ## is a valid DAG
  !isValidGraph(amat = amat, type = "cpdag"), ## is a valid CPDAG 
  isValidGraph(amat = amat, type = "pdag") ) ## is a valid PDAG

## a -- b -- c
amat <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3,3)
colnames(amat) <- rownames(amat) <- letters[1:3]
## plot(as(t(amat), "graphNEL"))             
stopifnot(!isValidGraph(amat = amat, type = "dag"), ## not a valid DAG
  isValidGraph(amat = amat, type = "cpdag"), ## is a valid CPDAG
  isValidGraph(amat = amat, type = "pdag") )## is a valid PDAG

## a -- b -- c -- d -- a
amat <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), 4,4)
colnames(amat) <- rownames(amat) <- letters[1:4]
## plot(as(t(amat), "graphNEL"))             
stopifnot(!isValidGraph(amat = amat, type = "dag"),  ## not a valid DAG
  !isValidGraph(amat = amat, type = "cpdag"), ## not a valid CPDAG
  !isValidGraph(amat = amat, type = "pdag") ) ## not a valid PDAG

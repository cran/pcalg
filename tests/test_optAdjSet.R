library(pcalg)
suppressWarnings(RNGversion("3.5.0"))

set.seed(123)
myDAG <- randomDAG(20, 0.3)
mcov <- trueCov(myDAG)
amat <- t(as(myDAG,"matrix"))
amat[which(amat!=0)] <- 1
graphEst <- dag2cpdag(myDAG)
## plot(graphEst)
amat.cpdag <- t(as(graphEst,"matrix"))

stopifnot(sort(optAdjSet(amat,2,6))==c())


stopifnot(pcalg:::isAmenable(amat.cpdag,5,16,type="cpdag"))
stopifnot(sort(optAdjSet(amat.cpdag,5,16))==c(1,2,3,4,6,7,8,11,13))
stopifnot(optAdjSet(amat.cpdag,5,16)==optAdjSet(amat,5,16))

stopifnot(pcalg:::isAmenable(amat.cpdag,2,12,type="cpdag"))
stopifnot(sort(optAdjSet(amat.cpdag,2,12))==c(5,7))
stopifnot(optAdjSet(amat.cpdag,2,12)==optAdjSet(amat,2,12))

stopifnot(optAdjSet(amat.cpdag,2,c(3,7,12))==optAdjSet(amat,2,12))
          

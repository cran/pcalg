library(pcalg)

# Y-structure MAG
# Encode as adjacency matrix
p <- 10 # total number of variables
V <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10") # variable labels
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,3,0,0,0,0,0,0,0,0),
              c(3,0,3,0,0,0,0,0,0,0),
              c(0,3,0,2,0,0,0,0,0,0),
              c(0,0,3,0,2,0,0,0,0,0),
              c(0,0,0,3,0,2,0,2,2,1),
              c(0,0,0,0,3,0,2,0,0,0),
              c(0,0,0,0,0,3,0,0,0,0),
              c(0,0,0,0,2,0,0,0,0,0),
              c(0,0,0,0,1,0,0,0,0,0),
              c(0,0,0,0,1,0,0,0,0,0))
rownames(amat)<-V
colnames(amat)<-V

stopifnot(all.equal(searchAM(amat,5,type = "an"), c(3,4,5))) # ancestors of X5
stopifnot(all.equal(searchAM(amat,5,type = "de"), c(5,6,7))) # descendants of X5
stopifnot(all.equal(searchAM(amat,5,type = "ant"), c(1,2,3,4,5))) # anteriors of X5
stopifnot(all.equal(searchAM(amat,5,type = "sp"), c(8))) # spouses of X5
stopifnot(all.equal(searchAM(amat,2,type = "nb"), c(1,3))) # neighbors of X2
stopifnot(all.equal(searchAM(amat,c(4,6),type = "pa"), c(3,5))) # parents of {X4,X6}
stopifnot(all.equal(searchAM(amat,c(3,5),type = "ch"), c(4,6))) # children of {X3,X5}
stopifnot(all.equal(searchAM(amat,5,type = "pde"), c(5,6,7,9,10))) # possible descendants of X5

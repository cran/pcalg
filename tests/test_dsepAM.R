library(pcalg)

# Y-structure MAG
# Encode as adjacency matrix
p <- 4 # total number of variables
V <- c("X1","X2","X3","X4") # variable labels
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,0,2,0),
              c(0,0,2,0),
              c(3,3,0,2),
              c(0,0,3,0))
rownames(amat)<-V
colnames(amat)<-V

stopifnot(dsepAM(1,2,S=NULL,amat) == TRUE)
stopifnot(dsepAM(1,2,S=4,amat) == FALSE)
stopifnot(dsepAM(1,2,S=c(3,4),amat) == FALSE)

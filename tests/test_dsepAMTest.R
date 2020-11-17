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

suffStat<-list(g=amat,verbose=FALSE)

cat('X1 d-separated from X2? ', dsepAMTest(1,2,S=NULL,suffStat),'\n') ## d-separated
cat('X1 d-separated from X2 given X4? ', dsepAMTest(1,2,S=4,suffStat),'\n') ## not d-separated given node 3
cat('X1 d-separated from X2 given X3 and X4? ', dsepAMTest(1,2,S=c(3,4),suffStat),'\n') ## not d-separated by node 3 and 4

# Derive PAG that represents the Markov equivalence class of the MAG with the FCI algorithm
# Make use of d-separation oracle as "independence test"
indepTest <- dsepAMTest
fci.pag <- fci(suffStat,indepTest,alpha = 0.5,labels = V,verbose=FALSE)

true.pag <- rbind(c(0,0,2,0),
                  c(0,0,2,0),
                  c(1,1,0,2),
                  c(0,0,3,0))
rownames(true.pag)<-V
colnames(true.pag)<-V

cat('True MAG:\n')
print(amat)
cat('PAG output by FCI:\n')
print(fci.pag@amat)
cat('True PAG:\n')
print(true.pag)

stopifnot(identical(true.pag,fci.pag@amat))

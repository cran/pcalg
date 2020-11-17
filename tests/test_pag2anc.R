library(pcalg)

##################################################
## Mooij et al. (2020), Fig. 43(a), p. 97
##################################################

# Encode ADMG as adjacency matrix
p <- 8 # total number of variables
V <- c("Ca","Cb","Cc","X0","X1","X2","X3","X4") # 3 context variables, 5 system variables
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j
amat <- rbind(c(0,2,2,2,0,0,0,0),
              c(2,0,2,0,2,0,0,0),
              c(2,2,0,0,2,2,0,0),
              c(3,0,0,0,0,0,2,0),
              c(0,3,3,0,0,3,0,2),
              c(0,0,3,0,2,0,0,0),
              c(0,0,0,3,0,0,0,2),
              c(0,0,0,0,2,0,3,0))
rownames(amat)<-V
colnames(amat)<-V

# Make use of d-separation oracle as "independence test"
indepTest <- dsepAMTest
suffStat<-list(g=amat,verbose=FALSE)

# Derive PAG that represents the Markov equivalence class of the ADMG with the FCI algorithm
# (assuming no selection bias)
fci.pag <- fci(suffStat,indepTest,alpha = 0.5,labels = V,verbose=TRUE,selectionBias=FALSE)

# Read off causal features from the FCI PAG
cat('Identified absence (-1) and presence (+1) of ancestral causal relations from FCI PAG:\n')
fci.anc <- pag2anc(fci.pag@amat)
print(fci.anc)

true.anc <- rbind(c( 1, 0, 0, 0, 0, 0, 0, 0),
                  c( 0, 1, 0, 0, 0, 0, 0, 0),
                  c( 0, 0, 1, 0, 0, 0, 0, 0),
                  c( 0, 0, 0, 1, 0, 0, 0, 0),
                  c(-1,-1,-1,-1, 1,-1,-1,-1),
                  c( 0, 0, 0, 0, 0, 1, 0, 0),
                  c( 0, 0, 0, 0, 0, 0, 1, 0),
                  c(-1,-1,-1,-1,-1,-1,-1, 1))
rownames(true.anc)<-V
colnames(true.anc)<-V

stopifnot(identical(true.anc,fci.anc))

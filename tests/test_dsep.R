#### Compare our dsep() with  dSep() from package "ggm" :
library(pcalg)

set.seed(22)
p <- 8
nreps <- 10
ok <- rep(FALSE,nreps)
for (i in 1:nreps) {
  myDAG <- randomDAG(p, prob = 0.3)
  amat <- as(myDAG,"matrix")
  amat[amat!=0] <- 1

  x <- sample(1:p,1)
  y <- sample(setdiff(1:p,x),1)
  S <- sample(setdiff(1:p,c(x,y)),sample(1:5,1))

  dsepOld <- ggm::dSep(amat,as.character(x),as.character(y),as.character(S))
  dsepRes <- dsep	   (as.character(x),as.character(y),as.character(S),
			    myDAG)
  ok[i] <- (dsepRes == dsepOld)
}

if (!all(ok)) stop("Test dsep wrong: dsep oracle made a mistake!")

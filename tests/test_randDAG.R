library(pcalg)
## setwd("/sfs/u/kalischm/research/packages/unifDAGs/")
## source("aux_general.R")
## source("randDAG.R")

## check weights
set.seed(123)
n <- 100
g <- randDAG(n=n,d=3, wFUN=list(runif,min=0,max=1))
m <- wgtMatrix(g)
v <- as.numeric(m)
v <- v[v!=0]
ct <- cut(x=v, breaks=seq(0,1,by=0.1))
pval <- chisq.test(as.numeric(table(ct)), p = rep(0.1,10))$p.value
if (pval < 0.05) {
    stop("test_unifDAG: Edge weights don't look uniformly distributed!")
}

## check generation of negative weights (fixed Bug)
set.seed(123)
tmp1 <- randDAG(3,2,wFUN = list(runif, min = 2, max = 2))
all( unlist(tmp1@edgeData@data) == 2 ) 
set.seed(123)
tmp2 <- randDAG(3,2,wFUN = list(runif, min = -2, max = -2))
all( unlist(tmp2@edgeData@data) == -2 ) 


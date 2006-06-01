library(pcalg)

set.seed(2)

## Parameters
p <- 5
n <- 1000
alpha <- 0.05
s <- 0.6

mlb <- 0.1
mub <- 1
dag <- randomDAG(p,s,lB=mlb,uB=mub)
data <- rmvDAG(n, dag)
## TODO: improve this, once we have print(),plot() methods:
r1 <- pcAlgo(data,alpha, verbose=TRUE)

if(dev.interactive()) {
    op <- par(mfrow=c(1,2))
    plot(dag,      main = "true")
    plot(r1@graph, main = "estimate")
    par(op)
}

## check TPR, FPR and TDR for randomly generated data
dag <- randomDAG(p,s,lB=mlb,uB=mub)
data <- rmvDAG(n, dag)
r <- pcAlgo(data,alpha,verbose=TRUE)
cgr <- compareGraphs(r@graph, dag)
if (abs(round(cgr["tpr"]-0.83333333333,10))>1e-10)
    stop("Inconsistent result in Testing test_pcAlgo.R")

if(dev.interactive()) {
    c.cgr <- format(round(cgr, 3))
    op <- par(mfrow=c(1,2))
    plot(dag, main = "true")
    plot(r@graph, main = "estimated skeleton")
    mtext(paste(paste(toupper(names(c.cgr)), c.cgr, sep="= "),
		collapse=", "))
    par(op)
}



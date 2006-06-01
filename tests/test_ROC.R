library(pcalg)

### quite a few simulations -- extended by MM

set.seed(42)
p <- 10
n <- 100
s <- 0.1
trueDAG <- randomDAG(p, s, lB=0.1, uB=1)
trueG <- ugraph(trueDAG)

ncor <- length(corMethods <- c("standard", "Qn"))
## using different alpha cutoffs for "standard" and "robust"
alpha <- cbind(alphaS =c(0.01, 0.05,  0.10, 0.20,  0.35, 0.50,  0.60, 0.70),
               alphaQ =c(0.001,0.005, 0.010,0.050, 0.100,0.200, 0.300,0.500))
colnames(alpha) <- corMethods
alphaL <- nrow(alpha)

errDs <- c("normal", "cauchy")

nReps <- 5
res <- array(0, dim= c(length(errDs), alphaL, ncor, 4),
               dimnames = list(err.Dist = errDs, alphas = rep("", alphaL),
               corrMeth = corMethods,
               stat = c("m.tpr", "SE.tpr", "m.fpr", "SE.fpr")))
tprMat <- matrix(0, nReps, ncor)
fprMat <- matrix(0, nReps, ncor)
## tdrMat <- matrix(0, nReps, ncor)
rMat <- matrix(0, alphaL, 4)
mix <- 0.1
for(errDist in errDs) {
    cat("\nError Dist: ", errDist,"\n")
    ## MM {FIXME?}: why don't keep both 'trueDAG' and 'dm' fixed for all j ?
    ## --          and really make 'j in 1:alphaL' to the *inner* loop?
    for (j in 1:alphaL) {
        cat("Outer: ",j," of ",alphaL," -> ", nReps, " inner : ", sep='')
        for (i in 1:nReps) {
            cat("",i,"")
            trueDAG <- randomDAG(p, s, lB=0.1, uB=1)
            trueG <- ugraph(trueDAG)
            dm <- rmvDAG(n, trueDAG, errDist)
            for (icor in 1:ncor) {
                estG <- pcAlgo(dm, alpha = alpha[j,icor],
                               corMethod = corMethods[icor])@graph
                cmp <- compareGraphs(estG, trueG)
                tprMat[i,icor] <- cmp["tpr"]
                fprMat[i,icor] <- cmp["fpr"]
                ## tdrMat[i,icor] <- cmp["tdr"]
            }
        }
        cat("\n")
        res[errDist, j, , ] <-
            cbind(colMeans(tprMat), apply(tprMat, 2, sd) / sqrt(nReps),
                  colMeans(fprMat), apply(fprMat, 2, sd) / sqrt(nReps))
    }
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

dim(res)
aperm(res, c(2,4,1,3))

## ROC curve plot "tpr" vs "fpr"

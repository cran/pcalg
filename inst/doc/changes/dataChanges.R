## Binary
datB <- as.matrix(read.table("/sfs/u/staff/kalisch/research/packages/pcalg/inst/doc/changes/dataBin.txt"))

mat <- matrix(c(0,1,0,0,0, 0,0,0,1,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,1,0,0),5,5)
colnames(mat) <- rownames(mat) <- as.character(1:5)
gBtrue <- as(mat, "graphNEL")
plot(gBtrue)

gmB <- list(x = datB, g = gBtrue)

## save(gmB, file = "/sfs/u/staff/kalisch/research/packages/pcalg/data/gmB.rda")

##################################################

## Gaussian
data(gaussianData)
set.seed(40)
p <- 8
n <- 5000
gGtrue <- randomDAG(p, prob = 0.3) ## true DAG
datG <- rmvDAG(n, gGtrue)

gmG <- list(x = datG, g = gGtrue)
save(gmG, file = "/sfs/u/staff/kalisch/research/packages/pcalg/data/gmG.rda")

##################################################

## Discrete
data(discreteData)
gmD <- list(x = datD, g = gDtrue)
save(gmD, file = "/sfs/u/staff/kalisch/research/packages/pcalg/data/gmD.rda")

##################################################

## IDA
data(idaData)
gmI <- list(x = datI, g = myDAG)
save(gmI, file = "/sfs/u/staff/kalisch/research/packages/pcalg/data/gmI.rda")

##################################################

## Latent Variable Graph
data(latentVarGraph)
gmL <- list(x = datL, g = g1)
save(gmL, file = "/sfs/u/staff/kalisch/research/packages/pcalg/data/gmL.rda")

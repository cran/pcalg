library("pcalg", lib.loc = "/u/kalisch/research/packages/pcalg.Rcheck")

## Collider with 2 in the middle
p <- 3
amat <- matrix(c(0,0,0, 1,0,1, 0,0,0), 3,3)
colnames(amat) <- rownames(amat) <- as.character(1:3)
g <- as(amat, "graphNEL")

plot(g)
wgtMatrix(g)

## oracle gets it right
indepTest <- dsepTest
suffStat <- list(g = g, jp = johnson.all.pairs.sp(g))
alpha <- 0.01 ## value is irrelevant as dsepTest returns either 0 or 1
fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)
plot(fit)

## generated data does not correspond to the DAG because rmvDAG can only be
## used on graphs that are topologically ordered
n <- 10^6
g@edgeData@data[[1]]$weight <- 1
g@edgeData@data[[2]]$weight <- 1

dat <- rmvDAG(n, g)
cor(dat)
indepTest <- gaussCItest 
## define sufficient statistics
suffStat <- list(C = cor(dat), n = n)
## estimate skeleton
alpha <- 0.01
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)
plot(pc.fit)

## Permutation
set.seed(123)
p <- 5
n <- 10^7
g <- randomDAG(p, 0.3)

dat <- rmvDAG(n, g)

perm <- sample(p)
indepTest <- gaussCItest 
## define sufficient statistics
suffStat <- list(C = cor(dat), n = n)
## estimate skeleton
alpha <- 0.05
pc.fit <- pc(suffStat, indepTest, p, alpha, verbose = TRUE)
par(mfrow = c(1,2))
plot(dag2cpdag(g))
plot(pc.fit)


## Test file ages
library(pcalg)
suppressWarnings(RNGversion("3.5.0"))
(doExtras <- pcalg:::doExtras())

## Known example where ges and ages output a different result
bool3 <- TRUE
set.seed(77)

p <- 8
n <- 5000
## true DAG:
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
gGtrue <- randomDAG(p, prob = 0.3, V = vars)
data = rmvDAG(n, gGtrue)

## Estimate the aggregated PDAG with ages
ages.fit <- ages(data = data)


## Estimate the essential graph with ges
score <- new("GaussL0penObsScore", data)
ges.fit <- ges(score)

diff <- as(ges.fit$essgraph,"matrix") - as(ages.fit$essgraph,"matrix")

bool3 <- ( bool3 & (diff[6,2]==1) * (diff[8,2]==1) *(sum(abs( as(ges.fit$essgraph,"matrix") - as(ages.fit$essgraph,"matrix") ))==2) )
stopifnot(bool3)

if (doExtras) {
## Test 1: Need to make sure that the skeleton of ges and ages are the same

bool1 <- TRUE
for(i in 1:10){
  
  p <- 20
  n <- 5000
  ## true DAG:
  gGtrue <- randomDAG(p, prob = 0.4)
  data <- rmvDAG(n, gGtrue)
  
  
  score <- new("GaussL0penObsScore", data=data, lambda = 200)
  ges.fit <- ges(score = score, phase = c("forward", "backward"), iterate = F)
  ages.fit <- ages(data = data, lambda_min = 200)
  
  
  ges.mat <- as(ges.fit$essgraph,"matrix")
  ages.mat <- as(ages.fit$essgraph,"matrix")
  # 
  # ages1mat2 <- ages1[[3]][[length(ages1[[3]])]]
  
  
  bool1 <- ( bool1 & all((((ges.mat + t(ges.mat))!=0)*1)==(((ages.mat + t(ages.mat))!=0)*1)) )
}
stopifnot(bool1)

## Test 2: Number of CPDAGs used and number of penalty parameter has to be the same
bool2 <- TRUE
for(i in 1:10){
  
  p <- 20
  n <- 5000
  ## true DAG:
  gGtrue <- randomDAG(p, prob = 0.4)
  data = rmvDAG(n, gGtrue)
  
  agesfit <- ages(data = data, lambda_min = 100)
  
  bool2 <- ( bool2 & (length(agesfit$CPDAGsList) == length(agesfit$lambda)) )
}
stopifnot(bool2)

}

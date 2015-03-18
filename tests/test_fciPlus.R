library(pcalg)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

p <- 14 ## will get slow if p > 14; these settings take already a few minutes
prob <- 0.2

##################################################
## simulate graph; remove confounders
##################################################
error1 <- error2 <- 0
  seed <- 547 ## Seed 547: Dsep Link wird benoetigt
  set.seed(seed)
  g <- randomDAG(p, prob)  ## true DAG
  
  #get the adjacency matrix
  tmp <- as(g,"matrix")
  tmp[which(tmp != 0)]=1
  
  #first three confounders should be removed
  #series - vector of confounding nodes to be removed
  b <- tmp
  n <- length(b[1,]) ## nmb of nodes
  k <- 1
  counter <- 1
  series <- c(0,0,0)
  while(k<=n & counter<=3)
  { 
   
    if (sum(b[k,])>1) ## k has more than one child
    {
      series[counter] <- k ## mark node k; leave out later
      counter <- counter+1
    }
    k <- k+1
  }
   
  #for i = 1..p, 
  #vect[i] = new name of node i in the graph when series nodes are removed
  vect <- c(1:p)
  for(r in 1:p)
  {
    if (r>=series[1]& series[1]>0) 
      vect[r] <- vect[r]-1
    if (r>=series[2]& series[2]>0) 
      vect[r] <- vect[r]-1
    if (r>=series[3] & series[3]>0) 
      vect[r] <- vect[r]-1
  }
  vect[series] <- 0

  ##################################################
  ## prepare sufficient statistics for tests
  ##################################################
##   cat("i=",i,"\n") 
  ## define sufficient statistics (d-separation oracle)

  cov.mat1 <- trueCov(g)
  true.cov1 <- cov.mat1
  
  ## delete rows and columns belonging to confounder variables in series
  if (series[1]>0)
    true.cov1 <- cov.mat1[-(series[1]),-(series[1])]
  if (series[2]>0)
    true.cov1 <- true.cov1[-(series[2]-1),-(series[2]-1)]
  if (series[3]>0)
    true.cov1 <- true.cov1[-(series[3]-2),-(series[3]-2)]
  
## transform covariance matrix into a correlation matrix
true.corr1 <- cov2cor(true.cov1)
true.corr <- cov2cor(cov.mat1)
suffStat1 <- list(C = true.corr1, n = 10^9)
indepTest1 <- gaussCItest

## FCI+: neue Variante
fciplus.fit <- fciPlus(suffStat = suffStat1, indepTest = indepTest1, alpha = 0.99, p = p-(series[1]>0)-(series[2]>0)-(series[3]>0))
fciplus.amatNeu <- fciplus.fit@amat

## Fit FCI
fci.fit <- fci(suffStat = suffStat1, indepTest = indepTest1, alpha = 0.99,
               p = p-(series[1]>0)-(series[2]>0)-(series[3]>0),
               rules = rep(TRUE, 10), verbose = FALSE)
fcimatrix <- fci.fit@amat

## cat("Comparison of adjmat ok:", all(fciplus.amat == fcimatrix),"\n")
adjError <- as.numeric(!all(fciplus.amatNeu == fcimatrix))

if (adjError != 0) {
    stop("Test fciPlus wrong: PAG was not found correctly")
}

showProc.time()
    

library(pcalg)
suppressWarnings(RNGversion("3.5.0"))

p <- 10
ip <- seq_len(p) # 1:p = 1 2 .. p

one.equal <- function(x1, xx, tol = 1e-14) {
    stopifnot(length(x1) == 1L, length(xx) >= 1L)
    ## the one xx[] that is close to x1 :
    x. <- xx[which.min(abs(xx - x1))]
    (if(abs(x1) < 1e-100) abs(x.- x1) else abs(x./x1 - 1)) <= tol
}

for (i in 1:50) {
  cat("i=",i,"\n--==\n",sep="") # need to see where it fails (if ..)
  set.seed(i)
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.4)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)

  x  <- sample(ip, 1)
  y1 <- sample(setdiff(ip,  x),       1)
  y2 <- sample(setdiff(ip,c(x,y1)),   1)
  y3 <- sample(setdiff(ip,c(x,y1,y2)),1)
  ## plot(myCPDAG)
  eff.true1 <- causalEffect(myDAG, y1, x)
  eff.true2 <- causalEffect(myDAG, y2, x)
  eff.true3 <- causalEffect(myDAG, y3, x)
  ## cat("x=",x," y1=",y1," eff=",eff.true1,"\n")
  ## cat("x=",x," y1=",y2," eff=",eff.true2,"\n")

  ## cat("est1: "); print(eff.est1 <- ida (x, y1, mcov,myCPDAG, method="local"))
  ## cat("est2: "); print(eff.est2 <- ida (x, y2, mcov,myCPDAG, method="local"))
  ## cat("est3: "); print(eff.est3 <- ida (x, y3, mcov,myCPDAG, method="local"))
  eff.est1 <- ida (x, y1, mcov,myCPDAG, method="local")
  eff.est2 <- ida (x, y2, mcov,myCPDAG, method="local")
  eff.est3 <- ida (x, y3, mcov,myCPDAG, method="local")
  eff.estF <- idaFast(x, c(y1,y2,y3), mcov,myCPDAG)
  cat("estF:\n"); print(eff.estF)

  stopifnot(exprs = {
      one.equal(eff.true1, eff.est1)
      one.equal(eff.true2, eff.est2)
      one.equal(eff.true3, eff.est3)
      all.equal(unname(eff.estF), tol = 3e-15, # even tol=0  works on some
                rbind(eff.est1, eff.est2, eff.est3, deparse.level=0L))
  })
}

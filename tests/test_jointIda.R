library(pcalg)

## Create a weighted DAG --- Example from ../man/jointIda.Rd  (keep in sync !)
p <- 6
V <- as.character(1:p)
edL <- list(
  "1" = list(edges=c(3,4), weights=c(1.1,0.3)),
  "2" = list(edges=c(6),  weights=c(0.4)),
  "3" = list(edges=c(2,4,6),weights=c(0.6,0.8,0.9)),
  "4" = list(edges=c(2),weights=c(0.5)),
  "5" = list(edges=c(1,4),weights=c(0.2,0.7)),
  "6" = NULL)
myDAG <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed") ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

n <- 1000
## simulate Gaussian data from the true DAG
dat <- if (require("mvtnorm")) {
  set.seed(123)
  rmvnorm(n, mean=rep(0,p), sigma=covTrue)
} else readRDS(system.file(package="pcalg", "external", "N_6_1000.rds"))

## estimate CPDAG -- see  help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p = p, alpha = 0.01, u2pd="relaxed")

## Suppose that we know the true CPDAG and covariance matrix
m1 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="RRC")
m2 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="MCD") ### MCD needed a bugfix

##sorting for comparison
order.m1 <- m1[,order(m1[1,])]
order.m2 <- m2[,order(m2[1,])]
tempres1 <- all( all.equal(order.m1, order.m2) )
if (!tempres1) stop("There is a mismatch with jointIda RRC and MCD!")

## Instead of knowing the true CPDAG, it is enough to know only
## all jointly valid parent sets of the intervention variables ## to use RRC or MCD
ajv.pasets <- list(list(5, c(3,4)),
                   list(integer(0), c(3,4)),
                   list(3, c(3,4)))
m3 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="RRC")
m4 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="MCD") 

##sorting for comparison
order.m3 <- m3[,order(m3[1,])]
order.m4 <- m4[,order(m4[1,])]
tempres2 <- all( all.equal(order.m3, order.m4) )
if (!tempres2) stop("There is a mismatch with jointIda RRC and MCD!")

## From the true DAG, we can compute the true total joint effects
## using RRC or MCD
m5 <- jointIda(x.pos=c(1,2),y.pos=6,covTrue,graphEst=myDAG,technique="RRC")
m6 <- jointIda(x.pos=c(1,2),y.pos=6,covTrue,graphEst=myDAG,technique="MCD")

##sorting for comparison
order.m5 <- m5[,order(m5[1,])]
order.m6 <- m6[,order(m6[1,])]
tempres3 <- all( all.equal(order.m5, order.m6))
if (!tempres3) stop("There is a mismatch with jointIda  RRC and MCD!")

################################################################

mTrue <- rbind(c(0,0.99,0.99), c(0.4,0.4,0.4))
Rnd <- 7
res1 <- (all(round(order.m1,Rnd) == mTrue) &
         all(round(order.m2,Rnd) == mTrue) &
         all(round(order.m3,Rnd) == mTrue) &
         all(round(order.m4,Rnd) == mTrue) )

if(!res1) stop("Test in jointIda: True causal effects were not recovered!")


#############################################################################

## jointIda also works when x.pos has length 1 and in the following example
## it gives the same result as ida() (see Note)
##
## When the CPDAG is known
v1 <- sort(round(jointIda(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
                     technique="RRC"), Rnd))
v2 <- sort(round(jointIda(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
               technique="MCD"), Rnd))    ###EMA MCD problem as above
v3 <- sort(round(ida(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
          method="global"), Rnd))   ###EMA there is an issue in v3 with lm.cov ??? i get the error:

res2 <- (all(v1==v2) & all(v2==v3) & all(v3==v1))

if(!res2) stop("Test in jointIda: ida() and jointIda() don't produce the same result!")

#########################################
## test extract parent set
#########################################
# library(digest)
# 
# 
# ## compare outputs
# jidaCompare <- function(a,b) {
#   digest(a) ## has "integer (0)"
#   digest(b) ## has "named integer(0)"
#   for (i in 1:length(b)) {
#     for (j in 1:length(b[[i]])) {
#       if (length(b[[i]][[j]]) == 0) b[[i]][[j]] <- integer(0)
#     }
#   }
#   
#   for (i in 1:length(a)) {
#     for (j in 1:length(a[[i]])) {
#       if (length(a[[i]][[j]]) == 0) a[[i]][[j]] <- integer(0)
#     }
#   }
#   
#   all(unique(sort(sapply(a,digest))) == unique(sort(sapply(b,digest))))
# }
# 
# 
# ##testing
# ##b[[5]]
# ##a[[14]]
# 
# ##c <- a[[14]]
# ##c2 <- b[[5]]
# 
# ###testing
# ## possible_p <- c(seq(20,100,by=10))
# possible_p <- c(seq(20,50,by=5))
# 
# ## possible_neighb_size <- c(3:15)
# ## possible_neighb_size <- c(3:10)
# possible_neighb_size <- c(3:8)
# 
# ## or choose prob and generate neighborhood size 
# possible_prob <- seq(0.1,0.6,by=0.1)
# 
# ##for function msep
# 
# ##for joint IDA MCD:
# ##if you want to run locally set nc to 1!!
# nc <-  1
# registerDoParallel(nc) ## Request nc processors in parallel
# num_setings <-20
# num_rep <- 5
# 
# ##amount of bg knowledge in percent
# p_bg <- seq(0,1,by=0.1)
# 
# ##  start.time <- Sys.time()
# 
# seed1 <- 1 
# resFin <- foreach(r=seed1: num_setings, .packages = 'pcalg') %dopar%{
#   set.seed(r)
#   ##Then we can choose one setting:
#   
#   p <-  sample(possible_p,1)
#   neigh <- sample(possible_neighb_size,1)
#   prob <- round(neigh/(p-1),3)
#   res <- rep(FALSE, num_rep)
#   
#   cat("r=",r,", p=",p,", nb=",neigh,"\n")     
#   
#   ## then for every setting selected we can generate i.e. 20 graph
#   for (d in 1: num_rep){
#     ##i left the seed inside for now
#     ## i am not sure which seed option makes it easier for us
#     ## feel free to move this around
#     #seed <- sample(possible_seed,1)
#     #set.seed(seed)
#     
#     ##get graph
#     g <- pcalg:::randomDAG(p, prob)  ## true DAG
#     cpdag <- dag2cpdag(g) ## true CPDAG
#     
#     ## get adjacency matrix
#     cpdag.amat <- t(as(cpdag,"matrix"))
#     dag.amat <- t(as(g,"matrix"))
#     dag.amat[dag.amat != 0] <- 1
#     
#     #############################
#     ##                         ##    
#     ## NEW CODE STARTS HERE!!! ##
#     ##                         ##    
#     #############################
#     
#     amat.cpdag <- cpdag.amat
#     #plot(cpdag)
#     x <- sample(p,3)
#     
#     a <- pcalg:::extract.parent.sets(x,t(amat.cpdag)) ## extract valid joint parent sets using pcalg
#     ##      b <- extract.parent.sets(x,amat.cpdag)          ## and using my own function
#     b <- b.extract.parent.sets(x,amat.cpdag)          ## and using my own function
#     ## cat("l.a=",length(a), " - ","l.b=",length(b),"\n")
#     ## cat("jidaCompare:",jidaCompare(a,b),"\n")
#     res[d] <- jidaCompare(a,b)
#     ##cat(res[d],"\n")
#   }
#   all(res)
# }
# 
# 
# stopifnot(all(unlist(resFin)))
# 
# print(resFin)

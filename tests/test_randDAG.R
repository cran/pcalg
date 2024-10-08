## Sep 2024: Stopped using output because of volative behaviour
## of interER, power, geometric, barabasi

library(pcalg)
suppressWarnings(RNGversion("3.5.0"))
## setwd("/sfs/u/kalischm/research/packages/unifDAGs/")
## source("aux_general.R")
## source("randDAG.R")

### Check all methods: ----------------------------------------------

## MM hack: extract them from the randDAG() function definition
body. <- body(randDAG)
is.switch <- function(P) !is.symbol(P) && identical(as.symbol("switch"), P[[1]])
switchCall <- body.[vapply(body., is.switch, NA)][[1]]
stopifnot(identical(as.symbol("switch"), switchCall[[1]]))
(rDAGmeths <- names(switchCall)[-c(1:2, length(switchCall))])
rDAGall <- function(n, d, ...)
  sapply(rDAGmeths, function(meth) randDAG(n,d, method=meth, ...),
         simplify=FALSE)
set.seed(37)
rD.10.4 <- rDAGall(10, 4) # 2024-02: no "low-level warning" anymore
##     , warning = function(w) {
##         rDAG.warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
## ## with a low-level warning:
## ## IGNORE_RDIFF_BEGIN
## rDAG.warn
## ## IGNORE_RDIFF_END
## stopifnot(grepl("graph_molloy_.*Cannot shuffle graph", rDAG.warn))
##
rD.10.4 # looks ok
## Show, but ignore the package startup messages:
## IGNORE_RDIFF_BEGIN
stopifnot( require("graph") )
## IGNORE_RDIFF_END

stopifnot(vapply(rD.10.4, isDirected, NA),
          vapply(rD.10.4, inherits, NA, what="graph"))
## nice plot of all 8 :
# op <- par(mfrow=c(4,2))
# invisible(lapply(names(rD.10.4), function(nm) plot(rD.10.4[[nm]], main=nm)))
# par(op)


str(outs <- lapply(rD.10.4, leaves, "out"))
if(packageVersion("igraph") < "2.0") ## currently fails
stopifnot(identical( outs,
    list(er = "3", regular = c("1", "5", "6"), watts = c("3", "4", "6"),
         bipartite = c("1", "2", "5"), barabasi = c("4", "8"),
         geometric = c("4", "7"), power = c("4", "5", "9"),
         interEr = c("3", "7"))
))

str(ins <- lapply(rD.10.4, leaves, "in"))
if(packageVersion("igraph") < "2.0") ## currently fails
stopifnot(identical( ins,
    list(er = c("1", "4", "7"), regular = c("3", "7", "10"),
         watts = c("1", "8"), bipartite = c("4", "6"),
         barabasi = c("6", "7"), geometric = c("5", "10"),
         power = c("2", "7"), interEr = c("8", "10"))
  ))

set.seed(47)
(rD.12.2 <- rDAGall(12, 2))
if(packageVersion("igraph") < "2.0") ## currently fails
stopifnot(exprs = {
    vapply(rD.12.2, isDirected, NA)
    vapply(rD.12.2, numNodes, 1) == 12
    identical(vapply(rD.12.2, numEdges, 1),
              setNames(c(9, 12, 12, 11, 11, 11, 13, 8), rDAGmeths))
})

##---------------------------------------------------------------------------

## Use the output here
require(Matrix)
# lapply(rD.10.4, function(g) as(as(g, "Matrix"),"nMatrix"))
# lapply(rD.12.2, function(g) as(as(g, "Matrix"),"nMatrix"))

## Minimal checks on graphs generated via igraph
stopifnot( require("graph") )
set.seed(37)
dagList <- vector(mode = "list", length = 6)
dagList[[1]] <- randDAG(10, 4, "regular")
dagList[[2]] <- randDAG(10, 4, "watts")
dagList[[3]] <- randDAG(10, 4, "er")
dagList[[4]] <- randDAG(10, 4, "bipartite")
dagList[[5]] <- randDAG(10, 4, "geometric")
dagList[[6]] <- randDAG(10, 4, "interEr", par2 = 0.5)

## number of nodes
stopifnot(all.equal(
  sapply(dagList, numNodes),
  rep(10,6)
))
## number of edges
stopifnot(all.equal(
  sapply(dagList, numEdges),
  c(20,20,16,15,10,18)
))

## Use the output here -- FIXME: check more
require(Matrix)
# lapply(rD.10.4, function(g) as(as(g, "Matrix"),"nMatrix"))
# lapply(rD.12.2, function(g) as(as(g, "Matrix"),"nMatrix"))

## check weights
set.seed(123)
n <- 100
g <- randDAG(n=n,d=3, wFUN=list(runif,min=0,max=1))
m <- wgtMatrix(g)
stopifnot(sum(m != 0) == 137)
v <- as.numeric(m)
v <- v[v!=0]
## dput(as.vector(summary(v, digits=7)))
stopifnot(all.equal(as.vector(summary(v, digits=7)),
                    c(0.008103577, 0.2589966, 0.5287397,
                      0.5232445, 0.8159941, 0.9915566)))
ct <- cut(x=v, breaks=seq(0,1,by=0.1))
stopifnot(all.equal(chisq.test(as.numeric(table(ct)), p = rep(0.1,10))$p.value,
                    0.3101796548))

## check generation of negative weights (fixed Bug)
set.seed(123) ; tmp1 <- randDAG(3,2, wFUN = list(runif, min =  2, max =  2))
set.seed(123) ; tmp2 <- randDAG(3,2, wFUN = list(runif, min = -2, max = -2))
stopifnot(unlist(tmp1@edgeData@data) ==  2,
	  unlist(tmp2@edgeData@data) == -2 )

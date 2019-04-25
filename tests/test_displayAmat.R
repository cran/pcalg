library(pcalg)
suppressWarnings(RNGversion("3.5.0"))
##################################################
## pcAlgo object
##################################################
## Load predefined data
data(gmG)
n <- nrow    (gmG8$x)
V <- colnames(gmG8$x)

## define sufficient statistics
suffStat <- list(C = cor(gmG8$x), n = n)
## estimate CPDAG
skel.fit <- skeleton(suffStat, indepTest = gaussCItest,
                     alpha = 0.01, labels = V)
(amSkel <- as(skel.fit, "amat"))
str(amSkel)
stopifnot(attr(amSkel, "type") == "cpdag",
          amSkel["Author", "Bar"] == 1,
          amSkel["Bar", "Author"] == 1,
          amSkel["Ctrl","Author"] == 0)

pc.fit <- pc(suffStat, indepTest = gaussCItest,
             alpha = 0.01, labels = V)
(amPC <- as(pc.fit, "amat"))
stopifnot(attr(amPC, "type") == "cpdag",
          amPC["V5", "V8"] == 0,
          amPC["V8", "V5"] == 1,
          amPC["Goal","Author"] == 0)

##################################################
## fciAlgo object
##################################################
set.seed(42)
p <- 7
## generate and draw random DAG :
myDAG <- randomDAG(p, prob = 0.4)

## find PAG using the FCI algorithm
myC <- cov2cor(trueCov(myDAG))
suffStat <- list(C = myC, n = 10^9)
V <- LETTERS[1:p] ## labels of nodes

fmFCI <- fci(suffStat, indepTest=gaussCItest, labels = V,
             alpha = 0.9999, doPdsep = FALSE)
(amFCI <- as(fmFCI, "amat"))
stopifnot(attr(amFCI, "type") == "pag",
          amFCI["B","E"] == 2,
          amFCI["C","D"] == 1,
          amFCI["G","A"] == 3)

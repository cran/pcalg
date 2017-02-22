####' Tests the causal inference algorithms for interventional data:
####' GIES, GES, DP
####'
####' @author Alain Hauser
####' $Id: test_gies.R 409 2017-02-05 15:41:51Z alhauser $

cat("Testing the causal inference algorithms for interventional data:\n")

library(pcalg)

source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> showProc.time(), assertError(), relErrV(), ...

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed
# str(gauss.data)
p <- ncol(gauss.data)

(doExtras <- pcalg:::doExtras())
DBG <- if(doExtras) TRUE else FALSE # no debugging by default
## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps) # = default for all.equal()

## Define all test settings
settings <- expand.grid(
    fcn = c("gies", "gds"),
    cpp = c(FALSE, TRUE),
    format = c("scatter", "raw"),
    stringsAsFactors = FALSE)
nreps <- 5

for (m in seq_along(settings)) {
  cat(sprintf("Algorithm: %s, C++: %s, storage format: %s\n", 
          settings$fcn[m], settings$cpp[m], settings$format[m]))
  
  for (i in 1:nreps) {
    perm <- 1:nrow(gauss.data)
    
    ## Randomly permute data
    if (i > 1) {
      set.seed(i)
      perm <- sample(perm)
    }
    
    score <- new("GaussL0penIntScore", 
                 targets = gauss.targets, 
                 target.index = gauss.target.index[perm], 
                 data = gauss.data[perm, ],
                 format = settings$format[m],
                 use.cpp = settings$cpp[m])
    f <- get(settings$fcn[m])
    est.graph <- f(score, verbose = DBG)
    for (j in 1:p) {
      if (!isTRUE(all.equal(est.graph$essgraph$.in.edges[[j]],
                            gauss.parents[[j]], tolerance = tol)))
        stop("Parents are not estimated correctly.")
    }
      showProc.time()
  }
  cat("[Ok]\n")
}

## Test compatibility with deprecated calling conventions
cat("Compatibility with deprecated calling conventions... ")
score <- new("GaussL0penIntScore", 
    targets = gauss.targets, 
    target.index = gauss.target.index, 
    data = gauss.data)

warningIssued <- FALSE
tryCatch(est.graph <- gies(p, gauss.targets, score),
    warning = function(w) warningIssued <<- TRUE)
if (!warningIssued) {
  stop("No warning issued for old calling conventions.")
} else {
  for (j in 1:p) {
    if (!isTRUE(all.equal(est.graph$essgraph$.in.edges[[j]],
            gauss.parents[[j]], tolerance = tol)))
      stop("Parents are not estimated correctly.")
  }
}
warningIssued <- FALSE
tryCatch(est.graph <- gies(p = p, targets = gauss.targets, score = score),
    warning = function(w) warningIssued <<- TRUE)
if (!warningIssued) {
  stop("No warning issued for old calling conventions.")
}
cat("[OK]\n")


## Test stepwise execution of GIES
cat(if(doExtras)"\n\n", "GIES stepwise", if(doExtras)":\n" else ": ... ",
    if(doExtras) paste0(paste(rep("=", 14), collapse=""), "\n"),
    sep = "")
for (cpp in c(FALSE, TRUE)) {
  ## Randomly permute data
  for (i in 1:nreps) {
    perm <- 1:nrow(gauss.data)
    if (i > 1) {
      set.seed(i)
      perm <- sample(perm)
    }
    score <- new("GaussL0penIntScore",
        targets = gauss.targets,
        target.index = gauss.target.index[perm],
        data = gauss.data[perm, ],
        use.cpp = cpp)

    ## Stepwise execution
    essgraph <- new("EssGraph", nodes = as.character(1:p), score = score)
    cont <- TRUE
    while(cont) {
      cont <- FALSE
      while(essgraph$greedy.step("forward")) cont <- TRUE
      while(essgraph$greedy.step("backward")) cont <- TRUE
      while(essgraph$greedy.step("backward")) cont <- TRUE
    }
    for (i in 1:p) {
      if(doExtras) cat("  use.cpp = ", cpp,"; i = ", i, "\n", sep="")
      if (!isTRUE(all.equal(est.graph$essgraph$.in.edges[[i]],
              gauss.parents[[i]], tolerance = tol)))
        stop("Parents are not estimated correctly.")
    }
    showProc.time()
  }
}

cat(if(doExtras) "\n", "Done.\n")

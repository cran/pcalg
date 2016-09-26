#' Tests the calculation of BIC and MLE, as well as basic functions
#' the corresponding score class
#' 
#' @author Alain Hauser
#' $Id: test_bicscore.R 393 2016-08-20 09:43:47Z alhauser $

cat("Testing the calculation of the BIC score...\n")

library(pcalg)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed

## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps)

## Define all test settings
settings <- expand.grid(
  format = c("scatter", "raw"),
  cpp = c(FALSE, TRUE),
  stringsAsFactors = FALSE)
nreps <- 5

## Check data storage and calculation of scores
for (m in 1:nrow(settings)) {
  cat(sprintf("Setting: storage format = %s, C++ library = %s\n", 
          settings$format[m], settings$cpp[m]))
  
  for (i in 1:nreps) {
    perm <- 1:nrow(gauss.data)
    
    ## Randomly permute data
    if (i > 1) {
      set.seed(i)
      perm <- sample(perm)
    }    
    
    ## Create the score object with valid data
    score <- new("GaussL0penIntScore", 
        targets = gauss.targets, 
        target.index = gauss.target.index[perm], 
        data = gauss.data[perm, ],
        format = settings$format[m],
        use.cpp = settings$cpp[m],
        intercept = FALSE)
    
    # print(score$pp.dat)
  
    if (any(score$pp.dat$data.count != 1000)) {
      stop("The number of non-interventions are not calculated correctly.")
    }
  
    if (settings$format[m] == "scatter") {
      if (any(score$pp.dat$scatter.index != 1:5)) {
        stop("The indices of the scatter matrices are not calculated correctly.")
      }
    
      for (j in 1:5) {
        if (!isTRUE(all.equal(score$pp.dat$scatter[[score$pp.dat$scatter.index[j]]][1:5, 1:5], 
            gauss.scatter[[j]], 
            tolerance = tol))) {
          stop("The scatter matrices are not calculated correctly.")
        }
      }
    } # IF "scatter"
    
    for (j in 1:5) {
      if (!isTRUE(all.equal(gauss.loc.score[[j]], 
          score$local.score(j, gauss.parents[[j]]), 
          tolerance = tol))) {
        stop("The local score is not calculated correctly.")
      }
    }
    
    # print(lapply(1:5, function(i) score$local.fit(i, gauss.parents[[i]])))
    
    for (j in 1:5) {
      local.mle <- score$local.fit(j, gauss.parents[[j]])
      if (length(local.mle) != length(gauss.mle[[j]]) ||
          !isTRUE(all.equal(gauss.mle[[j]], local.mle, 
          tolerance = tol))) {
        stop("The local MLE is not calculated correctly.")
      }
    }
  }
}

## List targets in a non-unique way, check if representation is corrected for
## internal storage
temp.targets <- gauss.targets
temp.targets[[2]] <- rep(temp.targets[[2]], 4)
score <- new("GaussL0penIntScore",
    targets = temp.targets,
    target.index = gauss.target.index,
    data = gauss.data,
    format = "scatter",
    use.cpp = FALSE,
    intercept = FALSE)
stopifnot(isTRUE(all.equal(score$pp.dat$targets, gauss.targets)))
stopifnot(isTRUE(all.equal(score$pp.dat$target.index, gauss.target.index)))

## Try to create the score object with non-valid data,
## check if error is thrown
stopifnot(isTRUE(
    tryCatch(
        score <- new("GaussL0penIntScore", 
            targets = gauss.targets,
            target.index = gauss.target.index),
        error = function(e) {
          cat(paste("  Error caught:", e$message, "\n", sep = " "))
          TRUE
        }
    )))

set.seed(307)
temp.targets <- gauss.targets
temp.targets <- c(temp.targets, temp.targets[[6]])
temp.target.index <- gauss.target.index
temp.target.index[sample(which(gauss.target.index == 6), size = 20)] <- length(temp.targets)
stopifnot(isTRUE(
    tryCatch(
        score <- new("GaussL0penIntScore", 
            targets = temp.targets,
            target.index = temp.target.index,
            data = gauss.data),
        error = function(e) {
          cat(paste("  Error caught:", e$message, "\n", sep = " "))
          TRUE
        }
    )))

temp.targets <- gauss.targets
temp.targets[[2]] <- c(temp.targets[[2]], 9)
stopifnot(isTRUE(
    tryCatch(
        score <- new("GaussL0penIntScore", 
            targets = temp.targets, 
            target.index = gauss.target.index,
            data = gauss.data),
        error = function(e) {
          cat(paste("  Error caught:", e$message, "\n", sep = " "))
          TRUE
        }
    )))

temp.target.index <- gauss.target.index
temp.target.index[1] <- length(gauss.targets) + 1
stopifnot(isTRUE(
    tryCatch(score <- new("GaussL0penIntScore", 
            targets = gauss.targets, 
            target.index = temp.target.index,
            data = gauss.data),
        error = function(e) {
          cat(paste("  Error caught:", e$message, "\n", sep = " "))
          TRUE
        }
    )))

## Test calculation of BIC score for discrete data
# TODO use more sophisticated data set...
discr.data <- cbind(c(3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4),
    c(5,5,5,5,7,7,7,7,7,7,5,5,5,5,5,7,7,7,7,7),
    c(1,1,9,8,1,1,8,8,9,9,1,1,9,9,9,1,1,1,9,9))
score <- new("DiscrL0penIntScore", data = discr.data)
# score$local.score(1, integer(0))

cat("Done.\n")

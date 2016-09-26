####' Tests adaptive versions of GES (ARGES and ARGES-skeleton)
####'
####' @author Alain Hauser
####' $Id: test_arges.R 393 2016-08-20 09:43:47Z alhauser $

cat("Testing adaptive versions of GES:\n")

library(pcalg)
library(graph)

## Test with DAG of 3 vertices

# Create DAG with 3 vertices a shielded v-structure (A --> B <-- C, A --> C).
# The edge weight of A --> C should be smaller than the others.
# Only allowing edges between A and B and B and C at the beginning, 
# we can check whether ARGES also allows an edge between A and C in the
# end.
dag <- new("GaussParDAG", 
    nodes = as.character(1:3), 
    in.edges = list(integer(0), c(1, 3), 1),
    params = list(c(0.8, 0), c(0.2, 0, 0.7, 1.2), c(0.6, 0, 0.1)))
cpdag <- dag2cpdag(dag)
adjMat <- as(cpdag, "matrix")

# Simulate data
n <- 5000
set.seed(307)
X <- rmvnorm.ivent(n, dag)

# Create a score object
score <- new("GaussL0penObsScore", X)

# Estimate DAG without restriction
ges.fit <- ges(score)
stopifnot(all.equal(adjMat, as(ges.fit$essgraph, "matrix")))

# Test old calling convention of GES
warningIssued <- FALSE
tryCatch(ges.fit <- ges(3, score),
    warning = function(w) warningIssued <<- TRUE)
stopifnot(warningIssued)

# Force a gap between vertices 1 and 3
fixedGaps <- matrix(FALSE, 3, 3)
fixedGaps[1, 3] <- fixedGaps[3, 1] <- TRUE
ges.fit <- ges(score, fixedGaps = fixedGaps)
adjMat <- matrix(FALSE, 3, 3)
adjMat[1, 2] <- adjMat[3, 2] <- TRUE
stopifnot(all.equal(adjMat, as(ges.fit$essgraph, "matrix")))

# Test ARGES (adaptive = 'vstructures')
arges.fit <- ges(score, fixedGaps = fixedGaps, adaptive = "vstructures")
adjMat <- as(cpdag, "matrix")
stopifnot(all.equal(adjMat, as(arges.fit$essgraph, "matrix")))

# Checking ARGES-skeleton (adaptive = 'triples')
# Create a new DAG of the form A --> B --> C, A --> C, where the edge weight
# of A --> C is weaker than the other edge weights
dag <- new("GaussParDAG", 
    nodes = as.character(1:3), 
    in.edges = list(integer(0), 1, 1:2),
    params = list(c(0.8, 0), c(0.4, 0, 0.7), c(0.4, 0, 0.1, 0.6)))
cpdag <- dag2cpdag(dag)
adjMat <- as(cpdag, "matrix")

# Simulate data
set.seed(307)
X <- rmvnorm.ivent(n, dag)

# Make score object
score <- new("GaussL0penObsScore", X)

# Fitting with a restriction (forbid edge A -- C)
fixedGaps <- matrix(FALSE, 3, 3)
fixedGaps[1, 3] <- fixedGaps[3, 1] <- TRUE

ges.fit <- ges(score, fixedGaps = fixedGaps)
adjMat[1, 3] <- adjMat[3, 1] <- FALSE
stopifnot(all.equal(adjMat, as(ges.fit$essgraph, "matrix")))

# Test ARGES
arges.fit <- ges(score, fixedGaps = fixedGaps, adaptive = "vstructures")
stopifnot(all.equal(adjMat, as(arges.fit$essgraph, "matrix")))

# Test ARGES-skeleton: should reproduce perfect fit
arges.fit <- ges(score, fixedGaps = fixedGaps, adaptive = "triples")
adjMat <- as(cpdag, "matrix")
stopifnot(all.equal(adjMat, as(arges.fit$essgraph, "matrix")))

cat("Done.\n")

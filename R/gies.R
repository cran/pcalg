## GIES algorithm
##
## Author: Alain Hauser <alain.hauser@bfh.ch>
## $Id: gies.R 498 2019-10-20 11:19:45Z alhauser $
###############################################################################

##################################################
## Auxiliary functions for simulations
##################################################

#' Randomly generates a Gaussian causal model >>>  ../man/r.gauss.pardag.Rd
#'
#' @param p number of vertices
#' @param prob probability of inserting an edge between two given
#'                    vertices
#' @param top.sort indicates whether the produced DAG should be
#'                    topologically sorted
#' @param normalize indicates whether weights and error variances
#'                    should be normalized s.t. the diagonal of the
#'                    corresponding covariance matrix is 1. Note that
#'                    weights and error variances can then lie outside
#'                    the boundaries specified below!
#' @param lbe lower bound of edge weights. Default: 0.1
#' @param ube upper bound of edge weights. Default: 1
#' @param neg.coef indicates whether also negative edge weights should
#'                    be sampled
#' @param labels
#' @param lbv lower bound of vertex variance. Default: 0.5
#' @param ubv upper bound of vertex variance. Default: 1
#' @return  an instance of gauss.pardag
r.gauss.pardag <- function(p,
    prob,
    top.sort = FALSE,
    normalize = FALSE,
    lbe = 0.1,
    ube = 1,
    neg.coef = TRUE,
    labels = as.character(1:p),
    lbv = 0.5,
    ubv = 1)
{
  ## Error checking
  stopifnot(is.numeric(p), length(p) == 1, p >= 2,
      is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
      is.numeric(lbe), is.numeric(ube), lbe <= ube,
      is.logical(neg.coef),
      is.numeric(lbv), is.numeric(ubv), lbv <= ubv,
      is.character(labels), length(labels) == p)

  ## Create list of nodes, edges and parameters
  edL <- as.list(labels)
  names(edL) <- labels

  ## Create list of parameters; first entry: error variances
  pars <- as.list(runif(p, min = lbv, max = ubv))
  names(pars) <- labels

  ## Create topological ordering
  top.ord <- if (top.sort) 1:p else sample.int(p)

  ## Sample edges and corresponding coefficients, respecting the generated
  ## topological ordering
  for (i in 2:p) {
    ii <- top.ord[i]
    parentCount <- rbinom(1, i - 1, prob)
    edL[[ii]] <- top.ord[sample.int(i - 1, size = parentCount)]
    weights <- runif(parentCount, min = lbe, max = ube)
    if (neg.coef)
      weights <- weights * sample(c(-1, 1), parentCount, replace = TRUE)
    pars[[ii]] <- c(pars[[ii]], 0, weights)
  }
  edL[[top.ord[1]]] <- integer(0)
  pars[[top.ord[1]]] <- c(pars[[top.ord[1]]], 0)

  ## Create new instance of gauss.pardag
  result <- new("GaussParDAG", nodes = labels, in.edges = edL, params = pars)

  ## Normalize if requested
  if (normalize) {
    H <- diag(result$cov.mat())
    result$set.err.var(result$err.var() / H)
    H <- sqrt(H)
    for (i in 1:p)
      if (length(edL[[i]]) > 0)
        result$.params[[i]][-c(1, 2)] <- pars[[i]][-c(1, 2)] * H[edL[[i]]] / H[i]
  }

  ## Validate and return object
  validObject(result)
  result
}

#' Simulates independent observational or interventional data for a
#' specified interventions from a Gaussian causal model
#'
#' @param   n         number of data samples
#' @param   object    an instance of gauss.pardag
#' @param   target    intervention target
#' @param   target.value    value of intervention targets
rmvnorm.ivent <- function(n, object, target = integer(0), target.value = numeric(0))
{

  p <- object$node.count()
  ## Error checking
  stopifnot(length(target) == 0 || (1 <= min(target) && max(target) <= p))
  stopifnot((is.vector(target.value) && length(target.value) == length(target)) ||
            (is.matrix(target.value) && dim(target.value) == c(n, length(target))))

  ## Simulate error terms
  sigma <- sqrt(object$err.var())
  mu <- object$intercept()
  Y <- matrix(rnorm(n*p, mu, sigma), nrow = p, ncol = n)

  ## Insert intervention values
  Y[target, ] <- target.value

  ## Calculate matrix of structural equation system
  A <- - t(object$weight.mat(target))
  diag(A) <- 1.

  ## Solve linear structural equations
  t(solve(A, Y))
}


##################################################
## Structure learning algorithms
##################################################

##' Wrapper function for all causal inference algorithms.  It's not recommended
##' to use it directly; adapted wrapper functions for the single algorithms are
##' provided
#'
##' @param algorithm 	name of the causal inference algorithm to be used
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets 	unique list of targets. Normally determined from the scoring object
##' @param ... 		additional parameters passed to the algorithm chosen
caus.inf <- function(
    algorithm = c("GIES", "GDS", "SiMy"),
    score,
    labels = score$getNodes(),
    targets = score$getTargets(),
    ...)
{
  algorithm <- match.arg(algorithm)

  # Catching error occurring when a user called one of the causal
  # inference algorithms using the old calling conventions: try to
  # rearrange passed arguments, print a warning
  #
  # NOTE: old calling conventions were
  # (algorithm, p, targets, score) for caus.inf
  # (p, targets, score) for all functions allowing interventional data
  # (p, score) for GES
  if (is.numeric(score)) {
    # This happens when the old calling convention is used with all
    # mandatory arguments unnamed
    p <- score
    if (is.list(labels) && is(targets, "Score")) {
      score <- targets
      targets <- labels
      labels <- as.character(1:p)
      warning(paste("You are using a DEPRECATED calling convention for",
              "gies(), gds() or simy(); please refer to the documentation",
              "of these functions to adapt to the new calling conventions."))
    } else if (is(labels, "Score")) {
      score <- labels
      labels <- as.character(1:p)
      warning(paste("You are using a DEPRECATED calling convention for",
              "ges(); please refer to the documentation",
              "to adapt to the new calling convention."))
    }
  } else if (is.numeric(labels) && length(labels) == 1) {
    # This happens when the old calling convention is used with only the
    # 'score' argument named
    labels <- as.character(1:labels)
    warning(paste("You are using a DEPRECATED calling convention for",
            "gies(), ges(), gds() or simy(); please refer to the documentation",
            "of these functions to adapt to the new calling conventions."))
  }

  if (!is(score, "Score")) {
    stop("'score' must be of a class inherited from the class 'Score'.")
  }
  if (!is.character(labels)) {
    stop("'labels' must be a character vector.")
  }
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("'targets' must be a list of integer vectors.")
  }

  essgraph <- new("EssGraph", nodes = labels, targets = targets, score = score)
  if (essgraph$caus.inf(algorithm, ...)) {
    if (algorithm == "GIES") {
      ## GIES yields an essential graph; calculate a representative thereof
      list(essgraph = essgraph, repr = essgraph$repr())
    } else {
      ## GDS and SiMy yield a DAG; calculate the corresponding essential graph,
      ## although calculations may come from a model class where Markov equivalence
      ## does not hold!
      list(essgraph = dag2essgraph(essgraph$repr(), targets = targets),
           repr = essgraph$repr())
    }
  } else stop("invalid 'algorithm' or \"EssGraph\" object")
}

##' Greedy Interventional Equivalence Search - GIES --> ../man/gies.Rd
##'
##' @param score	scoring object to be used
##' @param labels	node labels
##' @param targets	unique list of targets. Normally determined from the scoring object
##' @param fixedGaps	logical matrix indicating forbidden edges
##' @param adaptive sets the behaviour for adaptiveness in the forward phase (cf. "ARGES")
##' @param phase  lists the phases that should be executed
##' @param iterate  indicates whether the phases should be iterated. iterated = FALSE
##'   means that the required phases are run just once
##' @param turning	indicates whether the turning step should be included (DEPRECATED).
##' @param maxDegree	maximum vertex degree allowed
##' @param verbose	indicates whether debug output should be printed
##' @param ...		additional parameters (currently none)
gies <- function(
    score,
    labels = score$getNodes(),
    targets = score$getTargets(),
    fixedGaps = NULL,
    adaptive = c("none", "vstructures", "triples"),
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    turning = NULL,
    maxDegree = integer(0),
    verbose = FALSE,
    ...)
{
  # Catch calling convention of previous package versions:
  # ges(p, targets, score, fixedGaps = NULL, ...)
  # If this calling convention is used, issue a warning, but adjust the
  # arguments
  if (is.numeric(score) && is.list(labels) && inherits(targets, "Score")) {
    score <- targets
    targets <- labels
    labels <- as.character(1:length(score$getNodes()))
    warning(paste("You are using a deprecated calling convention for gies()",
            "which will be disabled in future versions of the package;",
            "cf. ?gies.", sep = " "))
  }
  # If the old calling convention was used with named arguments, "p = ..."
  # would assign a numerical value to "phase" (expanding arguments...)
  if (is.numeric(phase)) {
    phase <- c("forward", "backward", "turning")
    warning(paste("You are using a deprecated calling convention for gies()",
            "which will be disabled in future versions of the package;",
            "cf. ?gies.", sep = " "))
  }

  # Issue warning if argument 'turning' was used
  if (!missing(turning)) {
    stopifnot(is.logical(turning))
    warning(paste0("The argument 'turning' is deprecated; please use 'phase'",
                   "instead (cf. ?ges)"))

    if (turning) {
      phase <- c("forward", "backward", "turning")
      iterate <- FALSE
    } else {
      phase <- c("forward", "backward")
      iterate <- FALSE
    }
  }

  # Error checks
  if (!inherits(score, "Score")) {
    stop("Argument 'score' must be an instance of a class inherited from 'Score'.")
  }
  phase <- match.arg(phase, several.ok = TRUE)
  # TODO extend...

  caus.inf(
      "GIES",
      score = score,
      labels = labels,
      targets = targets,
      fixedGaps = fixedGaps,
      adaptive = adaptive,
      phase = phase,
      iterate = iterate,
      maxDegree = maxDegree,
      verbose = verbose,
      ...)
}

##' Greedy Equivalence Search - GES --> ../man/ges.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param fixedGaps 	logical matrix indicating forbidden edges
##' @param adaptive sets the behaviour for adaptiveness in the forward phase (cf. "ARGES")
##' @param phase  lists the phases that should be executed
##' @param iterate  indicates whether the phases should be iterated. iterated = FALSE
##'   means that the required phases are run just once
##' @param turning	indicates whether the turning step should be included (DEPRECATED).
##' @param maxDegree 	maximum vertex degree allowed
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
##' @param targets 	unique list of targets. Normally determined from the scoring object
ges <- function(
    score,
    labels = score$getNodes(),
    fixedGaps = NULL,
    adaptive = c("none", "vstructures", "triples"),
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    turning = NULL,
    maxDegree = integer(0),
    verbose = FALSE,
    ...)
{
  # Catch calling convention of previous package versions:
  # ges(p, score, fixedGaps = NULL, ...)
  # If this calling convention is used, issue a warning, but adjust the
  # arguments
  if (is.numeric(score) && inherits(labels, "Score")) {
    score <- labels
    labels <- as.character(1:length(score$getNodes()))
    warning(paste("You are using a deprecated calling convention for ges()",
            "which will be disabled in future versions of the package;",
            "please refer to the help page of ges().", sep = " "))
  }
  # If the old calling convention was used with named arguments, "p = ..."
  # would assign a numerical value to "phase" (expanding arguments...)
  if (is.numeric(phase)) {
    phase <- c("forward", "backward", "turning")
    warning(paste("You are using a deprecated calling convention for ges()",
            "which will be disabled in future versions of the package;",
            "cf. ?ges.", sep = " "))
  }

  # Issue warning if argument 'turning' was used
  if (!missing(turning)) {
    stopifnot(is.logical(turning))
    warning(paste0("The argument 'turning' is deprecated; please use 'phase'",
                   "instead (cf. ?ges)"))

    if (turning) {
      phase <- c("forward", "backward", "turning")
      iterate <- FALSE
    } else {
      phase <- c("forward", "backward")
      iterate <- FALSE
    }
  }

  # Error checks
  if (!inherits(score, "Score")) {
    stop("Argument 'score' must be an instance of a class inherited from 'Score'.")
  }
  phase <- match.arg(phase, several.ok = TRUE)
  # TODO extend...

  if(min(score$pp.dat$data.count) <= score$pp.dat$vertex.count){
      warning("The data set is high-dimensional, ges might not be
able to terminate")
  }

  caus.inf(
      "GIES",
      score = score,
      labels = labels,
      targets = list(integer(0)),
      fixedGaps = fixedGaps,
      adaptive = adaptive,
      phase = phase,
      iterate = iterate,
      maxDegree = maxDegree,
      verbose = verbose,
      ...)
}

##' Greedy DAG Search - GDS : greedy search in the DAG space --> ../man/gds.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets
##' @param fixedGaps 	logical matrix indicating forbidden edges
##' @param phase  lists the phases that should be executed
##' @param iterate  indicates whether the phases should be iterated. iterated = FALSE
##'   means that the required phases are run just once
##' @param turning	indicates whether the turning step should be included (DEPRECATED).
##' @param maxDegree 	maximum vertex degree allowed
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
gds <- function(
    score,
    labels = score$getNodes(),
    targets = score$getTargets(),
    fixedGaps = NULL,
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    turning = TRUE,
    maxDegree = integer(0),
    verbose = FALSE,
    ...)
{
  # Issue warning if argument 'turning' was used
  # TODO: do not check whether 'turning' is false, but whether 'turning'
  # was provided as an argument.
  if (!turning) {
    phase <- c("forward", "backward")
    iterate <- FALSE
    warning(paste("The argument 'turning' is deprecated; please use 'phase' instead",
            "(cf. ?ges)", sep = " "))
  }

  phase <- match.arg(phase, several.ok = TRUE)

  caus.inf(
      "GDS",
      score = score,
      labels = labels,
      targets = targets,
      fixedGaps = fixedGaps,
      phase = phase,
      iterate = iterate,
      maxDegree = maxDegree,
      verbose = verbose,
      ...)
}

##' Dynamic programming approach of Silander and Myllimäki - SiMy --> ../man/simy.Rd
##'
##' @param score 	scoring object to be used
##' @param labels 	node labels
##' @param targets
##' @param verbose 	indicates whether debug output should be printed
##' @param ... 		additional parameters (currently none)
simy <- function(score, labels = score$getNodes(), targets = score$getTargets(),
                 verbose = FALSE, ...)
{
  caus.inf("SiMy", score = score, labels = labels, targets = targets, verbose = verbose, ...)
}


#' Converts a DAG to an (observational or interventional) essential graph
dag2essgraph <- function(dag, targets = list(integer(0))) {
  edgeListDAG <- inEdgeList(dag)
  edgeListEssGraph <- .Call("dagToEssentialGraph", edgeListDAG, targets)
  if (is.matrix(dag)) {
    p <- nrow(dag)
    result <- sapply(1:p, function(i) 1:p %in% edgeListEssGraph[[i]])
    rownames(result) <- rownames(dag)
    colnames(result) <- colnames(dag)
    result
  } else if (inherits(dag, "graphNEL")) {
    nodeNames <- nodes(dag)
    names(edgeListEssGraph) <- nodeNames
    result <- new("graphNEL",
        nodes = nodeNames,
        edgeL = lapply(edgeListEssGraph, function(v) nodeNames[v]),
        edgemode = "directed")
    reverseEdgeDirections(result)
  } else {
    new("EssGraph",
        nodes = dag$.nodes,
        in.edges = edgeListEssGraph,
        targets = targets)
  }
}

##################################################
## Active learning algorithms
##################################################

##' Optimal intervention targets
##'
##' @param essgraph (Observational or interventional) essential graph,
##'   represented by an EssGraph or a graphNEL object.
##' @param max.size Maximum size of intervention target; only 1 and the
##'   number of nodes of `essgraph` (the default, if not set) are supported.
##' @param use.node.names Indicates if the intervention target should be
##'   returned as a list of node names (if `TRUE`) or indices (if `FALSE`).
opt.target <- function(essgraph, max.size, use.node.names = TRUE) {
  # Test parameters.
  if (inherits(essgraph, "graphNEL")) {
    essgraph <- as(essgraph, "EssGraph")
  }
  if (!inherits(essgraph, "EssGraph")) {
    stop("`essgraph` must be an object of class EssGraph or graphNEL.")
  }
  p <- essgraph$node.count()
  if (missing(max.size)) {
    max.size <- p
  }
  if (!(max.size %in% c(1, p))) {
    stop("`max.size` must either be 1 or the number of nodes of `essgraph` (",
         p, "); actual value: ", max.size)
  }

  # Get the optimal intervention target.
  target <- essgraph$opt.target(max.size = max.size)
  if (use.node.names) {
    return(essgraph$.nodes[target])
  } else {
    return(target)
  }
}

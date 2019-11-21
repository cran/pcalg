##################################################
### Part 1 : S4 classes used by pc and r/fci
##################################################

## $Id: AllClasses.R 499 2019-11-18 20:08:36Z alhauser $

setClass("gAlgo",
         slots = c(call = "call",
                        n = "integer",
                        max.ord = "integer",
                        n.edgetests= "numeric",
                        sepset= "list",
                   pMax= "matrix"),
         contains = "VIRTUAL")


setClass("fciAlgo", contains = "gAlgo",
         slots = c(amat = "matrix", allPdsep = "list",
                   n.edgetestsPDSEP = "numeric", max.ordPDSEP = "integer"))

setClass("pcAlgo", contains = "gAlgo",
         slots = c(graph = "graph", zMin = "matrix")) ## zMin for compatibility

## Methods

##' Extract graph part of an R object:
setGeneric("getGraph", function(x) as(x,"graph"))
setMethod("getGraph", "matrix", function(x) as(x, "graphAM"))
if(FALSE) {## if we would importFrom("Matrix", ....) in NAMESPACE
    setMethod("getGraph", "sparseMatrix", function(x) as(x, "graphNEL"))
    setMethod("getGraph", "Matrix", function(x) as(x, "graphAM"))
}
setMethod("getGraph", "pcAlgo", function(x) x@graph)
setMethod("getGraph", "fciAlgo", function(x) as(x@amat, "graphAM"))

setOldClass("amat")# our Adjacency Matrics -- are S3 classes (but want some S4 methods)

##' as(*, "matrix") methods --- give the adjacency matrices   with a  "type"  attribute
##' as(*, "amat")   methods --- adjacency matrix class "amat" with a  "type"  attribute
setAs("pcAlgo", "matrix",
      function(from) structure(wgtMatrix(from@graph), type = "amat.cpdag"))
setAs("pcAlgo", "amat",
      function(from) structure(wgtMatrix(from@graph), class = "amat", type = "cpdag"))

setAs("fciAlgo", "matrix",
      function(from) structure(from@amat, type = "amat.pag"))
setAs("fciAlgo", "amat",
      function(from) structure(from@amat, class = "amat", type = "pag"))


##' auxiliary, hidden
show.pc.amat <- function(amat, zero.print, ...) {
    cat("\nAdjacency Matrix G:\n")
    print.table(amat, zero.print=zero.print, ...)
}

##' auxiliary, hidden
show.fci.amat <- function(amat, zero.print, ...) {
    cat("\nAdjacency Matrix G:",
        "G[i,j] = 1/2/3 if edge mark of edge i-j at j is circle/head/tail.",
        "", sep="\n")
    print.table(amat, zero.print=zero.print, ...)
}

print.amat <- function(x, zero.print = ".", ...) {
    stopifnot(is.character(typ <- attr(x, "type")), length(typ) == 1,
	      is.matrix(x), (d <- dim(x))[1] == d[2])
    cat(sprintf("Adjacency Matrix 'amat' (%d x %d) of type %s:\n",
		d[1], d[2], sQuote(typ)))
    ## TODO: if dimension is too large, e.g. use  Matrix::printSpMatrix2(x)
    print.table(x, zero.print=zero.print, ...)
    invisible(x)
}

setMethod("summary", "pcAlgo",
          function(object, amat = TRUE, zero.print = ".", ...) {
 	    cat("Object of class 'pcAlgo', from Call:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
 		"\n\nNmb. edgetests during skeleton estimation:\n", sep = "")
            cat("===========================================\n")
            cat("Max. order of algorithm: ", object@max.ord,
                "\nNumber of edgetests from m = 0 up to m =", object@max.ord,
                ": ", object@n.edgetests)
            g <- object@graph
            nbrs <- vapply(g@edgeL, function(x) length(x$edges), 1L)
            cat("\n\nGraphical properties of skeleton:\n")
            cat("=================================\n")
            cat("Max. number of neighbours: ", max(nbrs),
                "at node(s)", which(nbrs==max(nbrs)),
                "\nAvg. number of neighbours: ", mean(nbrs),"\n")
            if(amat)
                show.pc.amat(as(g, "matrix"), zero.print=zero.print)
          })

setMethod("summary", "fciAlgo",
	  function(object, amat = TRUE, zero.print = ".", ...) {
	    cat("Object of class 'fciAlgo', from Call:\n",
		paste(deparse(object@call), sep = "\n", collapse = "\n"), sep="")
            ## NB: fciPlus() result has *none* of these {apart from adj.mat}:
	    if(length(o.max <- object@max.ord)) {
		cat("\n\nNmb. edgetests during skeleton estimation:\n",
		    "===========================================\n", sep="")
		cat("Max. order of algorithm: ", o.max,
		    "\nNumber of edgetests from m = 0 up to m =", o.max,
		    ": ", object@n.edgetests)
	    }
	    if(length(o.maxP <- object@max.ordPDSEP)) {
		cat("\n\nAdd. nmb. edgetests when using PDSEP:\n=====================================")
		cat("\nMax. order of algorithm: ", o.maxP,
		    "\nNumber of edgetests from m = 0 up to m =", o.maxP,
		    ": ", object@n.edgetestsPDSEP)
	    }
	    if(length(object@sepset)) {
		myLength <- function(x) if(is.null(x)) NA_integer_ else length(x)
		cat("\n\nSize distribution of SEPSET:")
		myTab <- table(sapply(object@sepset,
				      function(x) vapply(x, myLength, 1L)),
			       useNA = "always")
		print(myTab)
	    }
	    if(length(object@allPdsep)) {
		cat("\nSize distribution of PDSEP:")
		print(table(vapply(object@allPdsep, length, 1L)))
	    }
	    ##
	    if(amat)
		show.fci.amat(object@amat, zero.print=zero.print)
	  })


print.pcAlgo <- function(x, amat = FALSE, zero.print = ".", ...) {
    cat("Object of class 'pcAlgo', from Call:\n",
        paste(deparse(x@call), sep = "\n", collapse = "\n"),
        "\n", sep="")
    A <- as(x@graph, "matrix")
    if(amat)
        show.pc.amat(A, zero.print=zero.print, ...)
    amat2 <- A + 2*t(A)
    ude <- sum(amat2 == 3)/2
    de <- sum(amat2 == 1)
    cat("Number of undirected edges: ", ude, "\n")
    cat("Number of directed edges:   ", de, "\n")
    cat("Total number of edges:      ", de + ude, "\n")
    invisible(x)
}
setMethod("show", "pcAlgo", function(object) print.pcAlgo(object))


print.fciAlgo <- function(x, amat = FALSE, zero.print = ".", ...) {
    cat("Object of class 'fciAlgo', from Call:\n",
        paste(deparse(x@call), sep = "\n", collapse = "\n"),
        "\n", sep="")
    if(amat)
	show.fci.amat(x@amat, zero.print=zero.print, ...)
    invisible(x)
}
setMethod("show", "fciAlgo", function(object) print.fciAlgo(object))

## -> ../man/pcAlgo-class.Rd
setMethod("plot", signature(x = "pcAlgo"),
          function(x, y, main = NULL, zvalue.lwd = FALSE,
                   lwd.max = 7, labels = NULL, ...)
	{
          check.Rgraphviz()

          if(is.null(main))
              main <- deparse(x@call)
          attrs <- nodeAttrs <- list()
          p <- numNodes(G <- x@graph)
          if (!is.null(labels)) {
              attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
              names(labels) <- nodes(G)
              nodeAttrs$label <- labels
          }

          if (zvalue.lwd && numEdges(G) != 0) {
	      lwd.mat <-
		  if(is.matrix(Z <- x@zMin) && all(dim(Z) == p)) Z
		  else ## from newer pc(): 'zMin' is deprecated there, but pMax corresponds:
		      qnorm(x@pMax/2, lower.tail=FALSE)
	      lwd.mat <- lwd.max * lwd.mat/max(lwd.mat)
	      z <- Rgraphviz::agopen(G, name = "lwdGraph",
				     nodeAttrs = nodeAttrs, attrs = attrs)
	      for (i in seq_along(z@AgEdge)) {
		  z@AgEdge[[i]]@lwd <- lwd.mat[as.integer(z@AgEdge[[i]]@head),
					       as.integer(z@AgEdge[[i]]@tail)]
	      }
              Rgraphviz::plot(z, main = main, ...)
          } else {
              Rgraphviz::plot(G, nodeAttrs = nodeAttrs, main = main,
                              attrs = attrs, ...)
          }
      })

setMethod("plot", signature(x = "fciAlgo"),
          function(x, y, main = NULL, ...)
      {
          check.Rgraphviz()

          if(is.null(main))
	      main <- deparse(x@call)
	  else ## see also below
	      warning("main title cannot *not* be set yet [Rgraphviz::plot() deficiency]")
          amat <- x@amat
          g <- as(amat,"graphNEL")
          nn <- nodes(g)
          p <- numNodes(g)
          ## n.edges <- numEdges(g) -- is too large:
          ## rather count edges such that  "<-->" counts as 1 :
          n.edges <- numGedges(amat)
          ahs <- ats <- rep("none", n.edges)
          nms <- character(n.edges)
          cmat <- array(c("0" = "none",   "1" = "odot",
                          "2" = "normal", "3" = "none")[as.character(amat)],
                        dim = dim(amat), dimnames = dimnames(amat))
          iE <- 0L
          for (i in seq_len(p-1)) {
              x <- nn[i]
              for (j in (i+1):p) {
                  y <- nn[j]
                  if (amat[x,y] != 0) {
                      iE <- iE + 1L
                      ahs[[iE]] <- cmat[x,y]
                      ats[[iE]] <- cmat[y,x]
                      nms[[iE]] <- paste0(x,"~",y)
                  }
              }
          }
          names(ahs) <- names(ats) <- nms
	  edgeRenderInfo(g) <- list(arrowhead = ahs, arrowtail = ats)
          ## XXX Sep/Oct 2010  --- still current -- FIXME ??
          ## XXX undid change by MM, since edge marks didn't work anymore
          ## XXX "known bug in Rgraphviz, but not something they may fix soon"
	  ## Rgraphviz::plot(g, main = main, ...)
          Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
      })

#######################################################
### Part 2 : Reference classes and Methods used by GIES
#######################################################

#' Auxiliary function bringing targets in a standard format.
#'
#' At the same time, the function checks if the targets are valid; if not,
#' it throws an exception.
#'
#' @param 	p				number of vertices
#' @param 	targets			list of (unique) targets
#' @param 	target.index	vector of target indices, or NULL
#' @return  depends on arguments:
#'   if target.index == NULL: list of sorted targets
#'   if target.index != NULL: list with two entries, "targets" and "target.index"
.tidyTargets <- function(p, targets, target.index = NULL) {
  stopifnot((p <- as.integer(p)) > 0)

  # Check and convert targets
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("Argument 'targets' must be a list of integer vectors.")
  }
  rawTargets <- lapply(targets, function(v) unique(sort(as.integer(v))))
  targets <- unique(rawTargets)
  if (length(targets) < length(rawTargets)) {
    stop("List of targets must be unique.")
  }
  allTargets <- unlist(targets)
  if (length(allTargets) > 0) {
    if (any(is.na(allTargets))) {
      stop("Argument 'targets' must not contain NAs.")
    }
    min.max <- range(allTargets)
    if (min.max[1] <= 0 || min.max[2] > p) {
      stop("Targets are out of range.")
    }
  }

  # Check validity of target index, if provided
  if (!is.null(target.index)) {
    if (!is.numeric(target.index)) {
      stop("Argument 'target.index' must be an integer vector.")
    }
    target.index <- as.integer(target.index)
    min.max <- range(target.index)
    if (min.max[1] <= 0 || min.max[2] > length(targets)) {
      stop("Target index is out of range.")
    }
    # target.index <- match(rawTargets, targets)[target.index]
  }

  # Return value
  if (is.null(target.index)) {
    targets
  } else {
    list(targets = targets, target.index = target.index)
  }
}

#' Create a list of targets and a vector of target indices out of a
#' matrix indicating interventions
#'
#' @param 	A		a n x p boolean matrix; A[i, j] is TRUE iff vertex j is intervened
#' 							in data point i
#' @return 	list with two entries, "targets" and "target.index".
#' 					targets is a list of unique intervention targets
#' 					target.index is a vector of size n; the intervention target of data point
#' 					i is given by targets[[target.index[i]]].
mat2targets <- function(A)
{
  stopifnot(is.matrix(A) && is.logical(A) && all(dim(A) > 0))

  targets.raw <- as.list(apply(A, 1, which))
  targets <- unique(targets.raw)
  list(targets = targets, target.index = match(targets.raw, targets))
}

#' Create a boolean "intervention matrix" out of a list of targets
#' and a vector of target indices.  Can be seen as the "inverse function"
#' of "mat2targets"
#'
#' @param 	p				number of vertices
#' @param 	targets			list of (unique) targets
#' @param 	target.index	vector of target indices
targets2mat <- function(p, targets, target.index)
{
  ## Check validity of targets :  targetList <-
  .tidyTargets(p, targets, target.index)

  res <- matrix(FALSE, nrow = length(target.index), ncol = p)
  for (i in seq_along(target.index))
    res[i, targets[[target.index[i]]]] <- TRUE
  res
}

#' Auxiliary function reading an edge list (as used in the constructors
#' of DAGs) out of an adjacency matrix or a graphNEL object
#' @param from adjacency matrix, graphNEL object, or object inherited
#'  from ParDAG
#' @return list of in-edges; length of list = number of vertices,
#' entries for i-th vertex = indices sources of in-edges
inEdgeList <- function(from)
{
  if (is.matrix(from)) {
    p <- nrow(from)
    stopifnot(p == ncol(from))
    lapply(1:p, function(i) which(from[, i] != 0))
  } else if(inherits(from, "graphNEL")) {
    nodeNames <- graph::nodes(from)
    edgeList <- lapply(graph::inEdges(from), function(v) match(v, nodeNames))
    names(edgeList) <- NULL
    edgeList
  } else if (length(grep(".*ParDAG", class(from)) == 1)) {
    from$.in.edges
  }else {
    stop(sprintf("Input of class '%s' is not supported.", class(from)))
  }
}

# TODO: for all reference classes, make sure the constructor also works
# without arguments; or find another way to make the $copy() method work...
# (The default implementation of the $copy() method calls the constructor
# without arguments)

##' Virtual base class for all parametric causal models.
##' The meaning of the "params" depends on the model used.
setRefClass("ParDAG",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        .params = "list"),

    validity = function(object) {
      if (anyDuplicated(object$.nodes))
        return("The node names must be unique")
      if (!all(sapply(object$.in.edges, is.numeric)))
        return("The vectors in 'in.edges' must contain numbers.")

      ## edgeRange <- range(unlist(object$.in.edges))
      ## if (object$edge.count() > 0 &&
      ##     (edgeRange[1] < 1 || edgeRange[2] > object$node.count()))
      ##   return("Invalid range of edge sources.")

      if (object$edge.count() > 0) {
          edgeRange <- range(unlist(object$.in.edges))
          if (edgeRange[1] < 1 || edgeRange[2] > object$node.count())
              return("Invalid range of edge sources.")
      }
      ## else return
      TRUE
    },

    methods = list(
        #' Constructor
        initialize = function(
            nodes,
            in.edges = replicate(length(nodes), integer(0), simplify = FALSE),
            params = list()) {
          .nodes <<- nodes

          .in.edges <<- lapply(in.edges, as.integer)
          names(.in.edges) <<- NULL
          for (i in seq_along(nodes)) {
            names(.in.edges[[i]]) <<- NULL
          }

          .params <<- params
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },

        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(sapply(.in.edges, length))
        },

        #' Yields the variable types.
        #' Function must be overridden in inherited classes!
        #'
        #' @param vertex vector of indices of the vertices for which the
        #' variable types should be reported. If vertex == NULL, the types of
        #' all variables are returned.
        var.type = function(vertex = NULL) {
          stop("var.type() is not implemented in this class.")
        },

        #' Yields the levels of the factor variables.
        #' Function must be overridden in inherited classes!
        #'
        #' @param vertex vector of indices of the vertices for which the
        #' variable types should be reported. If vertex == NULL, the types of
        #' all variables are returned.
        levels = function(vertex = NULL) {
          stop("var.type() is not implemented in this class.")
        },

        #' Simulates (draws a sample of) interventional (or observational) data
        simulate = function(n, target = integer(0), int.level = numeric(0)) {
          stop("simulate() is not implemented in this class.")
        },

        #' Fits parameters by MLE using a scoring object
        fit = function(score) {
          .params <<- score$global.fit(.self)
        }
        ),

    contains = "VIRTUAL")

#' Coercion to a graphNEL instance
setAs("ParDAG", "graphNEL",
      function(from) {
      edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
      names(edgeList) <- from$.nodes
      result <- new("graphNEL",
          nodes = from$.nodes,
          edgeL = edgeList,
          edgemode = "directed")
      return(reverseEdgeDirections(result))
    })

#' Coercion to a (logical) matrix
setAs("ParDAG", "matrix",
      function(from) {
          i.p <- seq_len(p <- from$node.count())
          in.edge <- from$.in.edges
	  vapply(i.p, function(i) i.p %in% in.edge[[i]], logical(p))
      })

#' Plot method (needs Rgraphviz to work!!)
setMethod("plot", "ParDAG",
    function(x, y, ...) {
      if (!validObject(x))
        stop("The parametric DAG model to be plotted is not valid")

      if (missing(y))
        y <- "dot"
      invisible(plot(as(x, "graphNEL"), y, ...))
    })


#' Virtual base class for all scoring classes
setRefClass("Score",
    contains = "VIRTUAL",

    fields = list(
        .nodes = "character",
        decomp = "logical",
        c.fcn = "character",
        pp.dat = "list",
        .pardag.class = "character"),

    validity = function(object) {
      ## Check if targets are valid (i.e., unique)
      targets.tmp <- object$pp.dat$targets
      for (i in seq(along = targets.tmp)) {
        targets.tmp[[i]] <- sort(unique(targets.tmp[[i]]))
        if (length(targets.tmp[[i]]) != length(object$pp.dat$targets[[i]]))
          return("Target variables must not be listed multiple times.")
      }
      if (length(unique(targets.tmp)) != length(targets.tmp)) {
        return("Targets must not be listed multiple times.")
      }

      ## Check node names
      if (anyDuplicated(object$.nodes)) {
        return("The node names must be unique")
      }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(
            targets = list(integer(0)),
            nodes = character(0),
            ...) {
          .nodes <<- nodes
          pp.dat$targets <<- .tidyTargets(length(nodes), targets)
        },

        #' Yields a vector of node names
        getNodes = function() {
          .nodes
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },

        #' Checks whether a vertex is valid
        #' @param vertex vector of vertex indices
        validate.vertex = function(vertex) {
          if (length(vertex) > 0) {
            stopifnot(all(is.whole(vertex)))
            min.max <- range(vertex)
            stopifnot(1 <= min.max[1] && min.max[2] <= node.count())
          }
        },

        #' Checks whether a vector is a valid list of parents
        validate.parents = function(parents) {
          validate.vertex(parents)
          stopifnot(anyDuplicated(parents) == 0L)
        },

        #' Creates an instance of the corresponding ParDAG class
        create.dag = function() {
          new(.pardag.class, nodes = .nodes)
        },

        #' Getter and setter function for the targets
        getTargets = function() {
          pp.dat$targets
        },

        setTargets = function(targets) {
          pp.dat$targets <<- lapply(targets, sort)
        },

        #' Creates a list of options for the C++ functions for the internal
        #' calculation of scores and MLEs
        c.fcn.options = function(DEBUG.LEVEL = 0) {
          list(DEBUG.LEVEL = DEBUG.LEVEL)
        },

        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents, ...) {
          stop("local.score is not implemented in this class.")
        },

        #' Calculates the global score of a DAG which is only specified
        #' by its list of in-edges
        global.score.int = function(edges, ...) {
          if (c.fcn == "none") {
            ## Calculate score in R
            sum(sapply(1:pp.dat$vertex.count,
                    function(i) local.score(i, edges[[i]], ...)))
          } else {
            ## Calculate score with the C++ library
            .Call("globalScore", c.fcn, pp.dat, edges, c.fcn.options(...), PACKAGE = "pcalg")
          }
        },

        #' Calculates the global score of a DAG
        global.score = function(dag, ...) {
          global.score.int(dag$.in.edges, ...)
        },

        #' Calculates a local model fit for a vertex and its parents
        local.fit = function(vertex, parents, ...) {
          if (!decomp) {
            stop("local.fit can only be calculated for decomposable scores.")
          } else {
            stop("local.fit is not implemented in this class.")
          }
        },

        #' Calculates a global model fit
        global.fit = function(dag, ...) {
          if (c.fcn == "none") {
            ## Calculate score in R
            if (decomp) {
              in.edge <- dag$.in.edges
              lapply(1:pp.dat$vertex.count,
                  function(i) local.fit(i, in.edge[[i]], ...))
            } else {
              stop("global.fit is not implemented in this class.")
            }
          } else {
            ## Calculate score with the C++ library
            .Call("globalMLE", c.fcn, pp.dat, dag$.in.edges, c.fcn.options(...),
                PACKAGE = "pcalg")
          }
        }
    )
)

setRefClass("DataScore",
    contains = "Score",

    validity = function(object) {
      ## Check whether data is available from all intervention targets
      if (!isTRUE(all.equal(sort(unique(object$pp.dat$target.index)),
                            seq_along(object$pp.dat$targets)))) {
        return("Data from all intervention targets must be available")
      }

      ## Check if dimensions of target.index and data conincide
      if (length(object$pp.dat$target.index) != nrow(object$pp.dat$data))
        return("Length of target index vector does not coincide with sample size.")

      return(TRUE)
    },

    methods = list(
        #' Constructor
        #'
        #' @param 	data 			data set, jointly interventional and observational.
        #' 							Can either be a matrix or a data frame (this might
        #' 							be different for inherited classes!)
        #' @param	targets 		unique list of targets represented in the data
        #' @param	target.index	index vector for targets of data rows
        #' @param	nodes			node labels
        #' Note: all arguments must have a default value for inheritance,
        #' see ?setRefClass; apart from that, the default values are meaningless
        initialize = function(data = matrix(1, 1, 1),
            targets = list(integer(0)),
            target.index = rep(as.integer(1), nrow(data)),
            nodes = colnames(data),
            ...) {
          ## Node names (stored in constructor of "Score"):
          ## if data has no column names, correct them
          if (is.null(nodes)) {
            nodes <- as.character(1:ncol(data))
          }
          targetList <- .tidyTargets(ncol(data), targets, target.index)
          callSuper(targets = targetList$targets, nodes, ...)

          ## Order by ascending target indices (necessary for certain scoring objects)
          if (is.unsorted(targetList$target.index)) {
            perm <- order(targetList$target.index)
          } else {
            perm <- seq_along(targetList$target.index)
          }

          ## Store pre-processed data
          # pp.dat$targets <<- lapply(targets, sort)
          pp.dat$target.index <<- targetList$target.index[perm]
          pp.dat$data <<- data[perm, ]
          pp.dat$vertex.count <<- ncol(data)


          ## Store list of index vectors of "non-interventions": for each vertex k,
          ## store the indices of the data points for which k has NOT been intervened
          A <- !targets2mat(pp.dat$vertex.count, pp.dat$targets, pp.dat$target.index)
          pp.dat$non.int <<- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
          # apply() cannot be used since we need a list, not a matrix.
          pp.dat$data.count <<- as.integer(colSums(A))
          pp.dat$total.data.count <<- as.integer(nrow(data))

          ## Declare scores as not decomposable "by default"
          decomp <<- FALSE

          ## No C++ scoring object by default
          c.fcn <<- "none"

          ## R function objects
          pp.dat$local.score <<- function(vertex, parents) local.score(vertex, parents)
          pp.dat$global.score <<- function(edges) global.score(vertex, parents)
          pp.dat$local.fit <<- function(vertex, parents) local.fit(vertex, parents)
          pp.dat$global.fit <<- function(edges) global.fit(vertex, parents)
        }
    )
)

#' l0-penalized log-likelihood for Gaussian models, with freely
#' choosable penalty lambda.
#' Special case: BIC where \lambda = 1/2 \log n (default value for lambda)
setRefClass("GaussL0penIntScore",
    contains = "DataScore",

    fields = list(
        .format = "character"),

    validity = function(object) {
      if (!is.null(object$pp.dat$scatter)) {
        ## Data storage with precalculated scatter matrices
        if (!isTRUE(all.equal(unique(object$pp.dat$scatter.index),
                              seq_along(object$pp.dat$scatter)))) {
          return("The index list of distinct scatter matrices has an invalid range.")
        }
        p <- ncol(object$pp.dat$data)
        if (any(sapply(object$pp.dat$scatter,
                       function(mat) !isTRUE(all.equal(dim(mat), c(p + 1, p + 1)))))) {
          return("The scatter matrices have invalid dimensions.")
        }
      }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(data = matrix(1, 1, 1),
            targets = list(integer(0)),
            target.index = rep(as.integer(1), nrow(data)),
            nodes = colnames(data),
            lambda = 0.5*log(nrow(data)),
            intercept = TRUE,
            format = c("raw", "scatter"),
            use.cpp = TRUE,
            ...) {
          ## Store supplied data in sorted form. Make sure data is a matrix for
          ## linear-Gaussian data
          if (!is.matrix(data)) {
            data <- as.matrix(data)
          }
          callSuper(data = data, targets = targets, target.index = target.index, nodes = nodes, ...)

          ## Number of variables
          p <- ncol(data)

          ## l0-penalty is decomposable
          decomp <<- TRUE

          ## Underlying causal model class: Gaussian
          .pardag.class <<- "GaussParDAG"

          ## Store different settings
          pp.dat$lambda <<- lambda
          pp.dat$intercept <<- intercept

          ## Store data format. Currently supporting scatter matrices
          ## and raw data only (recommended for high-dimensional data)
          .format <<- match.arg(format, several.ok = TRUE)

          ## If format not specified by user, choose it based on dimensions
          ## TODO: check if this choice is reasonable...
          if (length(.format) > 1) {
            .format <<- if(p >= nrow(data) && length(pp.dat$targets) > 1) "raw" else "scatter"
          }

          ## Use C++ functions if requested
          if (use.cpp) { ## now .format is of length one: if() is __much__ faster than ifelse()
            c.fcn <<- if(.format == "scatter") "gauss.l0pen.scatter" else "gauss.l0pen.raw"
          }

          ## Preprocess data if storage format is "scatter"; for "raw" format,
          ## everything is already available in pp.dat
          if (.format == "scatter") {
            ## Add column of ones to data matrix to calculate scatter matrices;
            ## this allows the computation of an intercept if requested
            data <- cbind(pp.dat$data, 1)# take matrix that is already pre-processed,
            # having reordered rows!

            ## Create scatter matrices for different targets
            ti.lb <- c(sapply(seq_along(pp.dat$targets), function(i) match(i, pp.dat$target.index)),
                length(pp.dat$target.index) + 1)
            scatter.mat <- lapply(seq_along(pp.dat$targets),
                function(i) crossprod(data[ti.lb[i]:(ti.lb[i + 1] - 1), , drop = FALSE]))

            ## Find all interventions in which the different variables
            ## are _not_ intervened
            non.ivent <- matrix(FALSE, ncol = p, nrow = length(pp.dat$targets))
            pp.dat$scatter.index <<- integer(p)
            max.si <- 0
            for (i in 1:p) {
              ## Generate indices of (distinct) scatter matrices
              non.ivent[ , i] <- sapply(seq_along(pp.dat$targets),
                                        function(j) i %nin% pp.dat$targets[[j]])
              pp.dat$scatter.index[i] <<- max.si + 1
              j <- 1
              while (j < i) {
                if (all(non.ivent[, i] == non.ivent[, j])) {
                  pp.dat$scatter.index[i] <<- pp.dat$scatter.index[j]
                  j <- i
                }
                j <- j + 1
              }
              if (pp.dat$scatter.index[i] == max.si + 1)
                max.si <- max.si + 1
            }

            ## Calculate the distinct scatter matrices for the
            ## "non-interventions"
            pp.dat$scatter <<- lapply(1:max.si,
               function(i) Reduce("+", scatter.mat[non.ivent[, match(i, pp.dat$scatter.index)]]))
          } # IF "scatter"
        },

        #' Calculates the local score of a vertex and its parents
        local.score = function(vertex, parents, ...) {
          ## Check validity of arguments
          validate.vertex(vertex)
          validate.parents(parents)

          if (c.fcn == "none") {
            ## Calculate score in R
            if (.format == "raw") {
              ## calculate score from raw data matrix
              ## Response vector for linear regression
              Y <- pp.dat$data[pp.dat$non.int[[vertex]], vertex]
              sigma2 <- sum(Y^2)

              if (length(parents) + pp.dat$intercept != 0) {
                ## Get data matrix on which linear regression is based
                Z <- pp.dat$data[pp.dat$non.int[[vertex]], parents, drop = FALSE]
                if (pp.dat$intercept)
                  Z <- cbind(1, Z)

                ## Calculate the scaled error covariance using QR decomposition
                Q <- qr.Q(qr(Z))
                sigma2 <- sigma2 - sum((Y %*% Q)^2)
              }
            }
            else if (.format == "scatter") {
              ## Calculate the score based on pre-calculated scatter matrices
              ## If an intercept is allowed, add a fake parent node
              parents <- sort(parents)
              if (pp.dat$intercept)
                parents <- c(pp.dat$vertex.count + 1, parents)

              pd.scMat <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]]
              sigma2 <- pd.scMat[vertex, vertex]
              if (length(parents) != 0) {
                b <- pd.scMat[vertex, parents]
                sigma2 <- sigma2 - as.numeric(b %*% solve(pd.scMat[parents, parents], b))
              }
            }

            ## Return local score
            return(-0.5*pp.dat$data.count[vertex]*(1 + log(sigma2/pp.dat$data.count[vertex])) - pp.dat$lambda*(1 + length(parents)))
          } else {
            ## Calculate score with the C++ library
            return(.Call("localScore", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
          } # IF c.fcn
        },

        #' Calculates the local MLE for a vertex and its parents
        #'
        #' @param 	vertex		vertex whose parameters shall be fitted
        #' @param 	parents		parents of the vertex
        #' @param 	...				ignored; for compatibility with the base class
        local.fit = function(vertex, parents, ...) {
          ## Check validity of arguments
          validate.vertex(vertex)
          validate.parents(parents)

          if (c.fcn == "none") {
            ## Calculate score in R
            if (.format == "raw") {
              ## Calculate MLE from raw data matrix
              ## Response vector for linear regression
              Y <- pp.dat$data[pp.dat$non.int[[vertex]], vertex]
              beta <- numeric(0)
              sigma2 <- sum(Y^2)

              ## Calculate regression coefficients
              if (length(parents) + pp.dat$intercept != 0) {
                ## Get data matrix on which linear regression is based
                Z <- pp.dat$data[pp.dat$non.int[[vertex]], parents, drop = FALSE]
                if (pp.dat$intercept)
                  Z <- cbind(1, Z)

                ## Calculate regression coefficients
                qrZ <- qr(Z)
                beta <- solve(qrZ, Y)

                ## Calculate the scaled error covariance using QR decomposition
                sigma2 <- sigma2 - sum((Y %*% qr.Q(qrZ))^2)
              }
            } else if (.format == "scatter") {
              ## Calculate MLE based on pre-calculated scatter matrices
              ## If an intercept is allowed, add a fake parent node
              parents <- sort(parents)
              if (pp.dat$intercept)
                parents <- c(pp.dat$vertex.count + 1, parents)

              pd.scMat <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]]
              sigma2 <- pd.scMat[vertex, vertex]
              if (length(parents) != 0) {
                beta <- solve(pd.scMat[parents, parents],
                    pd.scMat[vertex, parents])
                sigma2 <- sigma2 - pd.scMat[vertex, parents] %*% beta
              }
              else
                beta <- numeric(0)
            } # IF .format

            if (pp.dat$intercept) {
              return(c(sigma2/pp.dat$data.count[vertex], beta))
            } else {
              return(c(sigma2/pp.dat$data.count[vertex], 0, beta))
            }
          } else {
            ## Calculate score with the C++ library
            return(.Call("localMLE", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
          } # IF c.fcn
        }
        )
    )

##' Observational score as special case
setRefClass("GaussL0penObsScore", contains = "GaussL0penIntScore",

    methods = list(
        #' Constructor
        initialize = function(data = matrix(1, 1, 1),
            nodes = colnames(data),
            lambda = 0.5*log(nrow(data)),
            intercept = TRUE,
            format = c("raw", "scatter"),
            use.cpp = TRUE,
            ...) {
          callSuper(data = data,
              targets = list(integer(0)),
              target.index = rep(as.integer(1), nrow(data)),
              nodes = nodes,
              lambda = lambda,
              intercept = intercept,
              format = format,
              use.cpp = use.cpp,
              ...)
          }
        )
    )

#' Interventional essential graph
setRefClass("EssGraph",
    fields = list(
        .nodes = "vector",
        .in.edges = "list",
        .targets = "list",
        .score = "Score"
    ),

    validity = function(object) {
      ## Check in-edges
      if (!all(sapply(object$.in.edges, is.numeric))) {
        return("The vectors in 'in.edges' must contain numbers.")
      }
      if (!all(unique(unlist(object$.in.edges)) %in% 1:object$node.count())) {
        return(sprintf("Invalid edge source(s): edge sources must be in the range 1:%d.",
          object$node.count()))
      }

      ## Check targets
      if (anyDuplicated(object$.targets)) {
        return("Targets are not unique.")
      }
      if (!all(unique(unlist(object$.targets)) %in% 1:object$node.count())) {
        return(sprintf("Invalid target(s): targets must be in the range 1:%d.",
          object$node.count()))
      }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(nodes,
            in.edges = replicate(length(nodes), integer(0), simplify = FALSE),
            targets = list(integer(0)),
            score = NULL) {
          ## Store nodes names
          if (missing(nodes)) {
            stop("Argument 'nodes' must be specified.")
          }
          .nodes <<- as.character(nodes)

          ## Store in-edges
          stopifnot(is.list(in.edges) && length(in.edges) == length(nodes))
          # More error checking is done in validity check
          .in.edges <<- in.edges
          names(.in.edges) <<- NULL

          ## Store targets
          setTargets(targets)

          ## Store score
          setScore(score)
        },

        #' Yields the number of nodes
        node.count = function() {
          length(.nodes)
        },

        #' Yields the total number of edges in the graph
        edge.count = function() {
          sum(vapply(.in.edges, length, 1L))
        },

        #' Getter and setter functions for score object
        getScore = function() {
          .score
        },

        setScore = function(score) {
          if (!is.null(score)) {
            .score <<- score
          }
        },

        #' Getter and setter functions for targets list
        getTargets = function() {
          .targets
        },

        setTargets = function(targets) {
          .targets <<- lapply(targets, sort)
        },

        #' Creates a list of options for the C++ function "causalInference";
        #' internal function
        causal.inf.options = function(
            caching = TRUE,
            phase = c("forward", "backward", "turning"),
            iterate = length(phase) > 1,
            maxDegree = integer(0),
            maxSteps = 0,
            childrenOnly = integer(0),
            fixedGaps = NULL,
            adaptive = c("none", "vstructures", "triples"),
            verbose = 0,
            p = 0) {
          # Check for deprecated calling convention and issue a warning
          if (p > 0) {
            warning(paste("Argument 'p' is deprecated in calls of ges() or gies",
                    "and will be disabled in future package versions;",
                    "please refer to the corresponding help page.", sep = " "))
          }

          # Error checks for supplied arguments
          # TODO extend!
          if (is.logical(adaptive)) {
            adaptive <- ifelse(adaptive, "vstructures", "none")
            warning(paste("The parameter 'adaptive' should not be provided as logical anymore;",
                    "cf. ?ges or gies", sep = " "))
          }
          phase <- match.arg(phase, several.ok = TRUE)
          stopifnot(is.logical(iterate))
          adaptive <- match.arg(adaptive)
          if (is.null(fixedGaps)) {
            adaptive <- "none"
          }
          list(caching = caching,
              phase = phase,
              iterate = iterate,
              maxDegree = maxDegree,
              maxSteps = maxSteps,
              childrenOnly = childrenOnly,
              fixedGaps = fixedGaps,
              adaptive = adaptive,
              DEBUG.LEVEL = as.integer(verbose))
        },

        #' Performs one greedy step
        greedy.step = function(direction = c("forward", "backward", "turning"), verbose = FALSE, ...) {
          stopifnot(!is.null(score <- getScore()))

          ## Cast direction
          direction <- match.arg(direction)
          alg.name <- switch(direction,
              forward = "GIES-F",
              backward = "GIES-B",
              turning = "GIES-T")

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              alg.name,
              score$c.fcn,
              causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, ...),
              PACKAGE = "pcalg")
          if (identical(new.graph, "interrupt"))
            return(FALSE)

          if (new.graph$steps > 0) {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
          }

          return(new.graph$steps == 1)
        },

        greedy.search = function(direction = c("forward", "backward", "turning")) {
          stopifnot(!is.null(score <- getScore()))

          ## Cast direction
          direction <- match.arg(direction)
          alg.name <- switch(direction,
              forward = "GIES-F",
              backward = "GIES-B",
              turning = "GIES-T")

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              alg.name,
              score$c.fcn,
              causal.inf.options(caching = FALSE),
              PACKAGE = "pcalg")
          if (identical(new.graph, "interrupt"))
            return(FALSE)

          if (new.graph$steps > 0) {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
          }

          return(new.graph$steps)
        },

        #' Performs a causal inference from an arbitrary start DAG
        #' with a specified algorithm
        caus.inf = function(algorithm = c("GIES", "GIES-F", "GIES-B", "GIES-T", "GIES-STEP",
                "GDS", "SiMy"), ...) {
          stopifnot(!is.null(score <- getScore()))
          algorithm <- match.arg(algorithm)

          new.graph <- .Call("causalInference",
              .in.edges,
              score$pp.dat,
              algorithm,
              score$c.fcn,
              causal.inf.options(...),
              PACKAGE = "pcalg")

          if (identical(new.graph, "interrupt"))
            return(FALSE)
          else {
            .in.edges <<- new.graph$in.edges
            names(.in.edges) <<- .nodes
            return(TRUE)
          }
        },

        #' Yields a representative (estimating parameters via MLE)
        repr = function() {
          stopifnot(!is.null(score <- getScore()))

          result <- score$create.dag()
          result$.in.edges <- .Call("representative", .in.edges, PACKAGE = "pcalg")
          result$.params <- score$global.fit(result)

          return(result)
        },

        #' Calculates an optimal intervention target
        #'
        #' @param   max.size    maximum target size; allowed values: 1, p (= # nodes)
        opt.target = function(max.size) {
          .Call("optimalTarget", .in.edges, max.size, PACKAGE = "pcalg")
        }
        ))

##' Coercion to a graphNEL instance
.ess2graph <- function(from) {
    edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
    names(edgeList) <- from$.nodes
    reverseEdgeDirections(new("graphNEL",
                              nodes = from$.nodes,
                              edgeL = edgeList,
                              edgemode = "directed"))
}

setAs("EssGraph", "graphNEL", .ess2graph)
setAs("EssGraph", "graph", .ess2graph)

## NOTE: Coercion to SparseMatrix is more efficient via
## ----  via "graphNEL" for larger p :
##' Coercion to a (logical) matrix
setAs("EssGraph", "matrix",
      function(from) {
      ip <- seq_len(p <- from$node.count())
      vapply(ip, function(i) ip %in% from$.in.edges[[i]],
             logical(p))
    })

##' Coercion from a graph object to EssGraph (assuming an observational
##' essential graph).
setAs("graph", "EssGraph",
      function(from) {
        node.names <- nodes(from)
        new("EssGraph",
            nodes = node.names,
            in.edges = lapply(graph::inEdges(from),
                              function(v) match(v, node.names)))
      })

#' Plot method (needs Rgraphviz to work!!)
## TODO maybe adapt method to make sure that undirected edges are not plotted as
## bidirected
setMethod("plot", "EssGraph",
    function(x, y, ...) {
      if (!validObject(x))
        stop("Invalid parametric DAG model (\"EssGraph\")")
      if (missing(y))
        y <- "dot"
      invisible(plot(.ess2graph(x), y, ...))
    })

########################################################################

#' Gaussian causal model
setRefClass("GaussParDAG", contains = "ParDAG",

    validity = function(object) {
      if (!isTRUE(all.equal(sapply(object$.params, length),
                            sapply(object$.in.edges, length) + 2,
                            check.names = FALSE))) {
        return("The number of parameters does not match the number of in-edges.")
      }

      return(TRUE)
    },

    methods = list(
        #' Constructor
        initialize = function(
            nodes,
            in.edges = replicate(length(nodes), integer(0), simplify = FALSE),
            params = lapply(in.edges, function(l) numeric(length(l) + 2))) {
          callSuper(nodes, in.edges, params)
        },

        #' Yields the variable types.
        #'
        #' @param vertex vector of indices of the vertices for which the
        #' variable types should be reported. If vertex == NULL, the types of
        #' all variables are returned.
        var.type = function(vertex = NULL) {
          rep("numeric", if(is.null(vertex)) node.count else length(vertex))
        },

        #' Yields the levels of the factor variables. Always NULL in a Gaussian
        #' model.
        #'
        #' @param vertex vector of indices of the vertices for which the
        #' variable types should be reported. If vertex == NULL, the types of
        #' all variables are returned.
        levels = function(vertex = NULL) {
          vector("list", if(is.null(vertex)) node.count() else length(vertex))
        },

        #' Yields the intercept
        intercept = function() {
          sapply(.params, function(par.vec) par.vec[2])
        },

        #' Sets the intercept
        set.intercept = function(value) {
          for (i in 1L:node.count())
            .params[[i]][2] <<- value[i]
        },

        #' Yields the error variances
        err.var = function() {
          sapply(.params, function(par.vec) par.vec[1])
        },

        #' Sets the error variances
        set.err.var = function(value) {
          for (i in 1:node.count())
            .params[[i]][1] <<- value[i]
        },

        #' Yields the weight matrix w.r.t. an intervention target
        #'
        #' TODO add a method for sparse matrices...
        weight.mat = function(target = integer(0)) {
          ## Fill in weights
          p <- node.count()
          target <- as.integer(sort(target))
          result <- matrix(0, p, p)
          for (i in 1:p)
            if (as.integer(i) %nin% target)
              result[.in.edges[[i]], i] <- .params[[i]][-(1:2)]

          ## Set row and column names
          rownames(result) <- .nodes
          colnames(result) <- .nodes

          return(result)
        },

        #' Yields an observational or interventional covariance matrix
        #'
        #' @param   target    intervention target
        #' @param   ivent.var variances of the intervention variables
        #' @return  (observational or interventional) covariance matrix
        cov.mat = function(target = integer(0), ivent.var = numeric(0)) {
          A <- -weight.mat()
          A[, target] <- 0
          diag(A) <- 1
          A <- solve(A)

          all.var <- err.var()
          all.var[target] <- ivent.var

          return(t(A) %*% diag(all.var) %*% A)
        },

        #' Simulates (draws a sample of) interventional (or observational)
        #' data
        #'
        #' @param   n
        #' @param   target
        #' @param   int.level   intervention level: values of the intervened
        #'                      variables. Either a vector of the same length
        #'                      as "target", or a matrix with dimensions
        #'                      n x length(target)
        #' @return  a vector with the simulated values if n = 1, or a matrix
        #'          with rows corresponding to different samples if n > 1
        simulate = function(n, target = integer(0), int.level = numeric(0)) {
          ## Error terms, intercepts, and intervention levels
          if (n == 1) {
            Y <- rnorm(node.count(), mean = intercept(), sd = sqrt(err.var()))
          } else {
            Y <- matrix(rnorm(n*node.count(), mean = intercept(), sd = sqrt(err.var())), ncol = n)
          }
          if (length(target) > 0) {
            if (length(int.level) %nin% c(length(target), n*length(target))) {
              stop("int.level must either be a vector of the same length as target, or a matrix of dimension n x length(target)")
            }
            if (is.matrix(int.level)) {
              int.level <- t(int.level)
            }
            if (n == 1) {
              Y[target] <- int.level
            } else {
              Y[target, ] <- int.level
            }
          }

          ## Modified weight matrix (w.r.t. intervention target)
          D <- - t(weight.mat(target))
          diag(D) <- 1.

          ## Calculate results: simulation samples
          result <- solve(D, Y)
          if (n == 1) {
            result
          } else {
            t(result)
          }
        }
        )
    )

#' Coercion from a weight matrix
setAs("matrix", "GaussParDAG",
    def = function(from) {
      p <- nrow(from)
      stopifnot(p == ncol(from))

      if (!isAcyclic(from))
        stop("Input matrix does not correspond to an acyclic DAG.")

      nodes <- rownames(from)
      if (any(duplicated(nodes))) {
        warning("Row names are not unique; will reset node names.")
        nodes <- as.character(1:p)
      }
      if (is.null(nodes)) {
        nodes <- as.character(1:p)
      }

      edgeList <- inEdgeList(from)
      new("GaussParDAG",
          nodes = nodes,
          in.edges = edgeList,
          param = lapply(1:p, function(i) c(0, 0, from[edgeList[[i]], i])))
    })

#' Coercion from a "graphNEL" object
setAs("graphNEL", "GaussParDAG",
    def = function(from) {
      ## Perform coercion via weight matrix
      A <- as(from, "matrix")
      as(A, "GaussParDAG")
    })

#' Predict interventional or observational data points.  Intervention values
#' must be provided, the value of all non-intervened variables is calculated
#'
#' @param   object    an instance of GaussParDAG
#' @param   newdata   list with two entries:
#'                    target:     list of intervention targets (or single
#'                                intervention target)
#'                    int.level:  list of intervention levels (or single
#'                                vector of intervention levels)
#' @return  a matrix with rows containing the predicted values, or a vector,
#'          if a single prediction is requested
setMethod("predict", "GaussParDAG",
    function(object, newdata) {
      ## Check validity of parameters
      if (!validObject(object))
        stop("The parametric DAG model to be plotted is not valid")
      if (!is.list(newdata$target)) {
        if (is.list(newdata$int.level))
          stop("The two entries of newdata must both be vectors or both be lists.")
        newdata$target <- list(newdata$target)
        newdata$int.level <- list(newdata$int.level)
      }
      stopifnot(is.list(newdata$target),
          is.list(newdata$int.level),
          length(newdata$target) == length(newdata$int.level),
          all(sapply(newdata$target, length) == sapply(newdata$int.level, length)))

      if (length(newdata$target > 1))
        fit <- matrix(0, nrow = length(newdata$target), ncol = object$node.count())
      for (i in seq_along(newdata$target)) {
        ## Calculate predition for i-th target
        y <- object$intercept()
        y[newdata$target[[i]]] <- newdata$int.level[[i]]
        D <- -object$weight.mat(newdata$target[[i]])
        diag(D) <- 1.
        if (length(newdata$target > 1))
          fit[i, ] <- solve(D, y)
        else
          fit <- solve(D, y)
      }

      fit
    })


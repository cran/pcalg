fci <- function(suffStat, indepTest, alpha, labels, p,
                skel.method = c("stable", "original", "stable.fast"),
                type = c("normal", "anytime", "adaptive"),
                fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                maj.rule = FALSE, numCores = 1, selectionBias = TRUE,
                jci = c("0","1","12","123"), contextVars = NULL, 
                verbose = FALSE)
{
  ## Purpose: Perform FCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximum size of conditioning set
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - rules: array of length 10 wich contains TRUE or FALSE corresponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - doPdsep: compute possible dsep
  ## - biCC: TRUE or FALSE variable containing if biconnected components are
  ##         used to compute pdsep
  ## - conservative: TRUE or FALSE defining if
  ##          the v-structures after the pdsep
  ##          have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes
  ## - type: it specifies the version of the FCI that has to be used.
  ##         Per default it is normal, the normal FCI algorithm. It can also be
  ##         anytime for the Anytime FCI and in this cas m.max must be specified;
  ##         or it can be adaptive for Adaptive Anytime FCI and in this case
  ##         m.max must not be specified.
  ## - numCores: handed to skeleton(), used for parallelization
  ## - selectionBias: allow for selection bias (default: TRUE)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "12"  JCI assumptions 1 and 2 (i.e., no system variable causes any context variable,
  ##           and no system variable is confounded with any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## - contextVars: subset of variable indices that will be treated as context variables
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009; update: Diego Colombo, 2012; Martin Maechler, 2013; Joris Mooij, 2020

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")

  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")

  ## Check that jci background knowledge is valid
  jci <- match.arg(jci)
  ## Check whether contextVars is valid
  if( !is.null(contextVars) && length(contextVars) > 0 ) {
    if( !is.numeric(contextVars) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( !all(sapply(contextVars, function(i) i == as.integer(i))) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( min(contextVars) < 1 || max(contextVars) > p )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
  }
  ## Set fixed edges from JCI assumption 3 if asked for
  if( jci == "123" && length(contextVars) > 0 ) {
    if( any(is.null(fixedEdges)) )
      fixedEdges<-matrix(FALSE,p,p)
    fixedEdges[contextVars,contextVars]<-TRUE
  }

  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")

  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                     NAdelete=NAdelete, m.max=m.max, numCores=numCores, verbose=verbose)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL

  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
#    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
#                            alpha = alpha, version.unf = c(1,1),
#                            maj.rule = FALSE, verbose = verbose)
    pc.ci <- list(unfTripl = c(), sk = skel)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest = indepTest, p = p,
                      sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                      m.max = if (type == "adaptive") max.ordSKEL else m.max,
                      pdsep.max = pdsep.max, NAdelete = NAdelete,
                      unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                      biCC = biCC, fixedEdges = fixedEdges, verbose = verbose)

    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose = verbose, version.unf = c(1, 1),
                            maj.rule = maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose = verbose, version.unf = c(2, 1),
                                maj.rule = maj.rule)
      tripleList <- nopdsep$unfTripl
      ## update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if( !selectionBias )
    rules[5:7] <- FALSE
  else {
    if( jci != "0" )
      stop( 'The current JCI implementation does not support selection bias (use selectionBias=FALSE instead).' )
  }
  if (verbose)
    cat("\nDirect edges:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
                  jci = jci, contextVars = contextVars, verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} ## {fci}

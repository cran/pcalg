## only called in fci() [by default:  doPdsep=TRUE]
pdsep <- function (skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
                   pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
                   biCC = FALSE, fixedEdges = NULL, verbose = FALSE) ## FIXME: verbose : 2 --> qreach(verbose)
{
  ## Purpose: Compute Possible-D-SEP for each node, perform the condittional
  ##          independent tests and adapt graph accordingly
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - skel: Graph object returned by function skeleton
  ## - suffStat, indepTest: infofor the independence tests
  ## - p: number of nodes in the graph
  ## - sepset: Sepset that was used for finding the skeleton
  ## - alpha: niveau for the tests
  ## - pMax: Maximal p-values during estimation of skeleton
  ## - m.max: maximal size of the conditioning sets
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - unfVect: vector containing the unfaithful triples, used for the
  ##   conservative orientation of the v-structures
  ## - biCC: if the biconnected components have to be used
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G: Updated boolean adjacency matrix
  ## - sepset: Updated sepsets
  ## - pMax: Updated pMax
  ## - allPdsep: Possible d-sep for each node [list]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  9 Dec 2009
  ## Modification: Diego Colombo; Martin Maechler; Joris Mooij

  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  ord <- 0L
  allPdsep.tmp <- vector("list", p)
  if(biCC)
    conn.comp <- lapply(biConnComp(skel), as.numeric)

  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")

  if (any(G)) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    if (verbose) cat("\nCompute collider:\n")
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      for (z in allZ) {
        if (amat[x, z] == 0 &&
            !(y %in% sepset[[x]][[z]] ||
              y %in% sepset[[z]][[x]])) {

          if (length(unfVect) == 0) { ## normal version -------------------
            amat[x, y] <- amat[z, y] <- 2
            if (verbose) cat("\n",x,"*->", y, "<-*", z, "\n")
          }
          else { ## conservative version : check if x-y-z is faithful
            if (!any(unfVect == triple2numb(p,x,y,z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p,z,y,x), na.rm = TRUE)) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose)
                cat("\n",x,"*->", y, "<-*", z, "\n")
            }
          }
        }
      } ## for( z )
    } ## for( i  )
    allPdsep <- lapply(1:p, qreach, amat = amat)# verbose = (verbose >= 2)
    allPdsep.tmp <- vector("list", p)
    for(x in seq_len(p)) {
      if(verbose) cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]], "\n")
      if (any(an0 <- amat[x, ] != 0)) {
        tf1 <- setdiff(allPdsep[[x]], x)
        adj.x <- which(an0)
        for (y in adj.x)
          if( !fixedEdges[x,y] ) {
            if(verbose) cat(sprintf("\ny = %3d\n.........\n", y))
            tf <- setdiff(tf1, y)
            diff.set <- setdiff(tf, adj.x)
            ## bi-connected components
            if (biCC) {
              for(cci in conn.comp) {
                if (x %in% cci && y %in% cci)
                  break ## found it
              }
              bi.conn.comp <- setdiff(cci, c(x,y))
              tf <- intersect(tf, bi.conn.comp)
              if (verbose) {
                cat("There is an edge between",x,"and",y,"\n")
                cat("Possible D-Sep of", x,
                    "intersected with the biconnected component of",x,"and",y,
                    "is:", tf, "\n")
              }
            } ## if(biCC)
            allPdsep.tmp[[x]] <- c(tf,y) ## you must add y to the set
            ## for the large scale simulations, we need to stop the algorithm if
            ## it takes to much time, i.e. sepset>25
            if (length(tf) > pdsep.max) {
              if(verbose)
                cat("Size of Possible-D-SEP bigger than",pdsep.max,
                    ". Break the search for the edge between", x,"and",y,"\n")
            } else if (length(diff.set) > 0) {
              done <- FALSE
              ord <- 0L
              while (!done && ord < min(length(tf), m.max)) {
                ord <- ord + 1L
                if(verbose) cat("ord = ", ord, "\n")
                if (ord == 1) {
                  for (S in diff.set) {
                    pval <- indepTest(x, y, S, suffStat)
                    n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                    if (is.na(pval))
                      pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                    if (pval > pMax[x, y])
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                      done <- TRUE
                      if (verbose)
                        cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                      break
                    }
                  }
                }
                else { ## ord > 1
                  tmp.combn <- combn(tf, ord) ## has  choose( |tf|, ord ) columns
                  if (ord <= length(adj.x)) {
                    for (k in seq_len(ncol(tmp.combn))) {
                      S <- tmp.combn[, k]
                      if (!all(S %in% adj.x)) {
                        n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                        pval <- indepTest(x, y, S, suffStat)
                        if (is.na(pval))
                          pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                        if(pMax[x, y] < pval)
                          pMax[x, y] <- pval
                        if (pval >= alpha) {
                          amat[x, y] <- amat[y, x] <- 0
                          sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                          done <- TRUE
                          if (verbose)
                            cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                          break
                        }
                      }
                    } ## for(k ..)
                  }
                  else { ## ord > |adj.x| :
                    ## check all combinations; no combination has been tested before
                    for (k in seq_len(ncol(tmp.combn))) {
                      S <- tmp.combn[, k]
                      n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                      pval <- indepTest(x, y, S, suffStat)
                      if (is.na(pval))
                        pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                      if(pMax[x, y] < pval)
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                        done <- TRUE
                        if (verbose)
                          cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                        break
                      }
                    } ## for(k ..)
                  } ## else: { ord > |adj.x| }
                } ## else

              } ## while(!done ..)
            }

          } ## for(y ..)

      } ## if(any( . ))

    } ## for(x ..)
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE

  } ## if(any(G))

  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp,
       max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
} ## {pdsep}

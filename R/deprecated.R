pcAlgo <- function(dm = NA, C = NA, n = NA, alpha, corMethod = "standard",
                   verbose = FALSE, directed = FALSE,
                   G = NULL, datatype = 'continuous', NAdelete = TRUE,
                   m.max = Inf, u2pd = "rand", psepset = FALSE) {
  ## !!! DEPRECATED !!!
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdagu
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo

  .Deprecated(msg = "pcAlgo() is deprecated and only kept for backward compatibility.
 Please use skeleton, pc, or fci instead\n")
  cl <- match.call()

  if (any(is.na(dm))) {
    stopifnot(all(!is.na(C)),!is.na(n), (p <- ncol(C)) > 0)
  } else {
    n <- nrow(dm)
    p <- ncol(dm)
  }
  n <- as.integer(n)

  if (is.null(G)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else if (!(identical(dim(G),c(p,p))))
      stop("Dimensions of the dataset and G do not agree.")

  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  zMin <- matrix(Inf, p,p)
  n.edgetests <- numeric(1)# final length = max { ord}
  done <- FALSE
  ord <- 0

  if (datatype == 'continuous') {
    diag(zMin) <- 0
    if (any(is.na(C))) C <- mcor(dm, method = corMethod)
    cutoff <- qnorm(1 - alpha/2)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord+1] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[,1]), ]
      remEdges <- nrow(ind)
      if(verbose)
        cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
      for (i in 1:remEdges) {
        if(verbose && i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n")
        x <- ind[i,1]
        y <- ind[i,2]
        if (G[y,x]) {
          nbrsBool <- G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE
            S <- seq(length = ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord+1] <- n.edgetests[ord+1]+1
              z <- zStat(x,y, nbrs[S], C,n)
              if (verbose) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("z:",z,"\n"))
              if(abs(z) < zMin[x,y]) zMin[x,y] <- abs(z)
              if (abs(z) <= cutoff) {
                G[x,y] <- G[y,x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if(nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            }
          }
        } ## end if(!done)

      } ## end for(i ..)
      ord <- ord+1
      ##    n.edgetests[ord] <- remEdges
    } ## while

    for (i in 1:(p-1)) {
      for (j in 2:p) {
        zMin[i,j] <- zMin[j,i] <- min(zMin[i,j],zMin[j,i])
      }
    }
  }
  else {
    ##
    ##
    ## DISCRETE DATA ######################################################
    ##
    if (datatype == 'discrete') {
      dm.df <- as.data.frame(dm)
      while (!done && any(G) && ord <= m.max) {
        n.edgetests[ord+1] <- 0
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ## For comparison with C++ sort according to first row
        ind <- ind[order(ind[,1]), ]
        remEdges <- nrow(ind)
        if(verbose)
          cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
        for (i in 1:remEdges) {
          if(verbose) { if(i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n") }
          x <- ind[i,1]
          y <- ind[i,2]
          if (G[y,x]) {
            nbrsBool <- G[,x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord) done <- FALSE
              S <- seq(length = ord)
              repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
                n.edgetests[ord+1] <- n.edgetests[ord+1]+1
                prob <- ci.test(x,y, nbrs[S], dm.df)
                if (verbose) cat("x=",x," y=",y," S=",nbrs[S],":",prob,"\n")
                if (is.na(prob)) prob <- if(NAdelete) 1 else 0
                if(prob >= alpha) { # independent
                  G[x,y] <- G[y,x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                } else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if(nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          } ## end if(!done)

        } ## end for(i ..)
        ord <- ord+1
        ##    n.edgetests[ord] <- remEdges
      } ## while
    } else
      stop("Datatype must be 'continuous' or 'discrete'.")
  }

  if (psepset) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,] == 1),x) ## x-y-z

      for (z in allZ) {
        if (amat[x,z] == 0 &&
            !((y %in% sepset[[x]][[z]]) ||
                (y %in% sepset[[z]][[x]]))) {
          if (verbose >= 2) {
            cat("\n",x,"*->",y,"<-*",z,"\n")
            cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
          }

          ## x o-> y <-o z
          amat[x,y] <- amat[z,y] <- 2

        } ## for
      } ## if
    } ## for

    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,] != 0)) {
        tf1 <- setdiff(reach(x,-1,-1,amat), x)
        for (y in seq_p[amat[x,] != 0]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf) > 0) {
            az <- abs(zStat(x,y,tf,C,n))
            if (az < zMin[x,y]) zMin[x,y] <- az
            if (az <= cutoff) {
              ## delete x-y
              amat[x, y] <- amat[y, x] <- 0
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
            }
            if (verbose >= 2)
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - |z| = ",az,"\n")
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
  } ## end if(psepset)

  if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

  ## transform matrix to graph object (if not deprecated anyway: FIX to use correct node names!)
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = as.character(seq_p))
  } else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    as(G,"graphNEL")
  }

  res <- new("pcAlgo", graph = Gobject,
             call = cl, n = n, max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = zMin)
  if (directed)
    switch (u2pd,
            "rand"    = udag2pdag       (res),
            "retry"   = udag2pdagSpecial(res)$pcObj,
            "relaxed" = udag2pdagRelaxed(res))
  else
    res
} ## {pcAlgo} __ deprecated __

## DEPRECATED! -- use  ida() --
beta.special <- function(dat = NA, x.pos, y.pos, verbose = 0, a = 0.01,
                         myDAG = NA, myplot = FALSE, perfect = FALSE,
                         method = "local", collTest = TRUE, pcObj = NA, all.dags = NA, u2pd = "rand")
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj and all DAGs
  ## can be precomputed
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dat: data
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - verbose: 0=no comments, 1=progress in BB, 2=detail on estimates
  ## - a: significance level of tests for finding CPDAG
  ## - myDAG: needed if bootstrp==FALSE
  ## - myplot: plot estimated graph
  ## - perfect: True cor matrix is calculated from myDAG
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - collTest: True - Exclude orientations of undirected edges that
  ##   introduce a new collider
  ## - pcObj: Fit of PC Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  ## Covariance matrix: Perfect case / standard case
  if (perfect) {
    if(!is(myDAG, "graphNEL")) stop("For perfect-option the true DAG is needed!")
    mcov <- trueCov(myDAG)
    mcor <- cov2cor(mcov)
  } else {
    mcov <- cov(dat)
  }

  ## estimate skeleton and CPDAG of given data
  res <-
    if (is(pcObj, "pcAlgo"))
      pcObj
    else if(perfect)
      pcAlgo.Perfect(mcor, corMethod = "standard",directed = TRUE,u2pd = u2pd)
    else
      pcAlgo(dat, alpha = a, corMethod = "standard",directed = TRUE,u2pd = u2pd)

  ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
  amat <- ad.res <- wgtMatrix(res@graph)
  amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1

  if (method == "local") {
##############################
    ## local method
    ## Main Input: mcov
##############################
    ## find unique parents of x
    wgt.est <- ad.res
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x.pos,] != 0)
    if (y.pos %in% pa1) {
      ## x is parent of y -> zero effect
      beta.hat <- 0
    } else { ## y.pos not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x.pos,] != 0)
      if (verbose == 2) {
        cat("\n\nx=",x.pos,"y=",y.pos,"\n")
        cat("pa1=",pa1,"\n")
        cat("pa2=",pa2,"\n")
      }

      ## estimate beta
      if (length(pa2) == 0) {
        beta.hat <- lm.cov(mcov, y.pos, c(x.pos,pa1))
        if (verbose == 2)
          cat("Fit - y:",y.pos, "x:",c(x.pos,pa1), "|b.hat=", beta.hat)
      } else {
        beta.hat <- NA
        ii <- 1
        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        ## check for new collider
        if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
          beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
          if (verbose == 2)
            cat("\ny:",y.pos,"x:",c(x.pos,pa1),"|b.hat=", beta.hat[ii])
        }## else {
        ##   cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
        ## }
        ## exactly one member of pa2
        for (i2 in seq_along(pa2)) {
          ## check for new collider
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
            ii <-  ii+1
            if (y.pos %in% pa2.t) {
              ## cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
              beta.hat[ii] <- 0
            } else {
              beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
            }
            if (verbose == 2) { cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa2[i2]),
                  "|b.hat=",beta.hat[ii])
}
          } else {
            ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
          }
        }
        ## higher order subsets
        if (length(pa2) > 1) {
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))) {
              pa2.t <- pa.tmp[,j]
              pa2.f <- setdiff(pa2,pa2.t)
              ## teste auf neuen collider
              if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
                ii <- ii+1
                if (y.pos %in% pa2.t) {
                  cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                  beta.hat[ii] <- 0
                } else {
                  beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2.t))
                }
                if (verbose == 2) { cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa2.t),
                      "|b.hat=",beta.hat[ii])
}
              } else {
                ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
              }
            }
          }
        }
      } ## if pa2
    } ## if y in pa1
  } else {
##############################
    ## global method
    ## Main Input: mcov
##############################
    p <- numNodes(res@graph)
    am.pdag <- ad.res
    am.pdag[am.pdag != 0] <- 1
    ## find all DAGs if not provided externally
    if (is.na(all.dags)) {
      ## allDags(am.pdag,am.pdag,NULL)
      ad <- pdag2allDags(am.pdag)$dags
    } else {
      ad <- all.dags
    }
    n.dags <- nrow(ad)
    beta.hat <- rep.int(NA,n.dags)
    if (n.dags > 0) {
      if (myplot) {
        ## x11()
        par(mfrow = c(ceiling(sqrt(n.dags)), round(sqrt(n.dags)) ))
      }
      for (i in 1:n.dags) {
        ## compute effect for every DAG
        gDag <- as(matrix(ad[i,],p,p),"graphNEL")
        if (myplot) Rgraphviz::plot(gDag)
        ## path from y to x
        rev.pth <- RBGL::sp.between(gDag,as.character(y.pos),
                                    as.character(x.pos))[[1]]$path
        if (length(rev.pth) > 1) {
          ## if reverse path exists, beta=0
          beta.hat[i] <- 0
        } else {
          ## path from x to y
          pth <- RBGL::sp.between(gDag,as.character(x.pos),
                                  as.character(y.pos))[[1]]$path
          if (length(pth) < 2) {
            ## sic! There is NO path from x to y
            beta.hat[i] <- 0
          } else {
            ## There is a path from x to y
            wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
            pa1 <- which(wgt.unique[x.pos,] != 0)
            if (y.pos %in% pa1) {
              cat("Y in Parents: ",y.pos," in ",pa1,"\n")
              beta.hat[i] <- 0
            } else {
              beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
            }
            if (verbose == 2)
              cat("Fit - y:",y.pos,"x:",c(x.pos,pa1), "|b.hat=",beta.hat,"\n")
          } ## if length(pth)
        } ## if rev.pth
      } ## for n.dags
    } ## if n.dags
  } ## if method
  beta.hat
} ## {beta.special}

## DEPRECATED! -- use  ida() / idafast() --
beta.special.pcObj <- function(x.pos,y.pos,pcObj,mcov = NA,amat = NA,amatSkel = NA,
                               t.amat = NA)
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj has to be
  ## precomputed. This method is intended to be a fast version of
  ##
  ## beta.special(dat=NA,x.pos,y.pos,verbose=0,a=NA,myDAG=NA,myplot=FALSE,
  ## perfect=FALSE,method="local",collTest=TRUE,pcObj=pcObj,all.dags=NA,u2pd="relaxed")
  ##
  ## Thus, this is a faster version for the local method given a
  ## precomputed PC-Algo Object (relaxed udag2pdag, so CPDAG might not
  ## be a real CPDAG; this does not matter, since we try not to extend).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - pcObj: Fit of pc Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - mcov: covariance matrix of pcObj fit
  ## - amat,amatSkel,g2,t.amat are variants of the adjacency matrix that
  ##   are used internally but can be precomputed; the relevant code
  ##   is commented out
  ## ----------------------------------------------------------------------
  ## Value: List with two elements
  ## - beta.res: beta.causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  if (is.na(amat) | is.na(amatSkel) | is.na(t.amat)) {
    ## Code for computing precomputable variables
    ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
    amat <- wgtMatrix(pcObj@graph)
    amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
    t.amat <- t(amat)
    amatSkel <- amat + t.amat
    amatSkel[amatSkel != 0] <- 1
  }

  ## find unique parents of x
  tmp <- amat-t.amat
  tmp[which(tmp < 0)] <- 0
  wgt.unique <- tmp
  pa1 <- which(wgt.unique[x.pos,] != 0)
  if (y.pos %in% pa1) {
    cat("Y in Parents: ",y.pos," in ",pa1,"\n")
    beta.hat <- 0
  } else { ## y.pos not in pa1
    ## find ambiguous parents of x
    wgt.ambig <- amat-wgt.unique
    pa2 <- which(wgt.ambig[x.pos,] != 0)
    pa2 <- setdiff(pa2,y.pos)
    ## estimate beta
    if (length(pa2) == 0) {
      beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
    } else {
      beta.hat <- NA
      ii <- 1
      ## no member of pa2
      ## check for new collider
      pa2.f <- pa2
      pa2.t <- NA
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
      }
      ## exactly one member of pa2
      for (i2 in seq_along(pa2)) {
        ## check for new collider
        pa2.f <- pa2[-i2]
        pa2.t <- pa2[i2]
        if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
          ii <-  ii+1
          beta.hat[ii] <-
            if (y.pos %in% pa2.t) {
              cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
              0
            } else lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
        }
      }
      ## higher order subsets
      if (length(pa2) > 1) {
        for (i in 2:length(pa2)) {
          pa.tmp <- combn(pa2, i, simplify = TRUE)
          for (j in seq_len(ncol(pa.tmp))) {
            ## teste auf neuen collider
            pa2.t <- pa.tmp[,j]
            pa2.f <- setdiff(pa2,pa2.t)
            if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
              ii <- ii+1
              beta.hat[ii] <-
                if (y.pos %in% pa2.t) {
                  cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                  0
                } else lm.cov(mcov,y.pos,c(x.pos,pa1,pa2.t))
            }
          }
        }
      } ## if pa2
    } ## length(pa2)
  } ## y.pos %in% pa2
  beta.hat
} ## {beta.special.pcObj}

allDags <- function(gm,a,tmp, verbose = FALSE)
{
  ## Purpose: Find all DAGs for a given PDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix of initial PDAG; only 0-1 entries
  ##   i -> j iff gm(j,i)=1
  ## - a: copy of gm
  ## - tmp: "current set of DAGs", initially NULL
  ## ----------------------------------------------------------------------
  ## Value:
  ## - one 0/1 adj.matrix per row
  ## Reversion to graph: as(matrix(res[i,],p,p),"graphNEL")
  ## Reversion to wgtMatrix (i->j iff a[j,i]=1): t(matrix(res[i,],p,p))
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Apr 2008, 14:08
  .Deprecated(msg = "allDags() is deprecated and only kept for backward compatibility. Please use pdag2allDags() instead\n")
  if (sum(a) == 0) {
    if (verbose) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("#################### \n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink(a)
    if (verbose) {
      cat("Main Call: ################## \n")
        print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    for(x in sinks) {
      if (verbose) cat("Try removing", x," in a.\n")
      gm2 <- gm
      a2 <- a
      if (adj.check(a,x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(rownames(a)[inc.to.x])
          real.x <- as.numeric(rownames(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x,-x]
        if (verbose) {
          cat("Removed sink",as.numeric(rownames(a)[x]),
              "in g (", x,"in a).\n")
          cat("New graphs: \n")
          print(gm2)
          print(a)
        }
        tmp <- allDags(gm2, a2, tmp, verbose)
        ##     ------- *recursively*
      }
    }
  }
  tmp
}

## ___DEPRECATED__  rather plot(<fciAlgo>) --> setMethod("plot", "fciAlgo") in ./AllClasses.R
plotAG <- function(amat)
{
  ## Purpose: Plot ancestral graph
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix
  ##   amat[i,j]=3 & amat[j,i]=1 iff i 1-3 j
  ##   "0": no edge; "1": circle; "2": arrow; "3": tail
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 16 Feb 2009, 18:01
  check.Rgraphviz()

  g <- as(amat,"graphNEL")
  nn <- nodes(g)
  p <- numNodes(g)
  n.edges <- numEdges(g)
  ah.list <- at.list <- vector("list", n.edges)
  l.names <- character(n.edges)
  amat[amat == 1] <- "odot"
  amat[amat == 2] <- "normal"
  amat[amat == 3] <- "none"
  iE <- 0
  for (i in 1:(p-1)) {
    x <- nn[i]
    for (j in (i+1):p) {
      y <- nn[j]
      if (amat[x,y] != 0) {
        iE <- iE + 1
        ah.list[[iE]] <- amat[x,y]
        at.list[[iE]] <- amat[y,x]
        l.names[[iE]] <- paste0(x,"~",y)
      }
    }
  }
  names(ah.list) <- names(at.list) <- l.names

  edgeRenderInfo(g) <- list(arrowhead = ah.list, arrowtail = at.list)
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}

showAmat <- function(object) {
    .Deprecated(msg = "showAmat() is deprecated and only kept for backward compatibility.
 Please use as(*, \"amat\") instead\n")
  g <- getGraph(object)
  cat("\nAdjacency Matrix G:",
      "G[i,j] = 1/2 if edge mark of edge i-j at j is head/tail.",
      "", sep = "\n")
  wm <- wgtMatrix(g)
  mTmp <- t(wm - 2*t(wm))
  mTmp[ mTmp < 0 ] <- 2
  mTmp
}

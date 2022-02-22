udag2pag <- function(pag, sepset, rules = rep(TRUE,10), unfVect = NULL, jci = c("0","1","12","123"), contextVars = NULL, verbose = FALSE, orientCollider = TRUE)
{
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - pag: adjacency matrix of size pxp
  ## - sepset: list of all separation sets
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - unfVect: Vector with ambiguous triples (coded as number using triple2numb)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "12"  JCI assumptions 1 and 2 (i.e., no system variable causes any context variable,
  ##           and no system variable is confounded with any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## - contextVars: subset of variable indices that will be treated as context variables
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 6 Mar 2009; cleanup: Martin Maechler, 2010
  ## updates: Diego Colombo, 2012; Joris Mooij (2020)

  ## Notation:
  ## ----------------------------------------------------------------------
  ## 0: no edge
  ## 1: -o
  ## 2: -> (arrowhead)
  ## 3: - (tail)
  ## a=alpha
  ## b=beta
  ## c=gamma
  ## d=theta

  stopifnot(is.logical(rules), length(rules) == 10)
  if(!is.numeric(pag))  storage.mode(pag) <- "numeric"

  jci <- match.arg(jci)
  if( jci != "0" && any(rules[5:7]) )
    warning('FCI orientation rules 5-7 should be disabled if JCI functionality is used.')

  if (any(pag != 0)) {
    p <- as.numeric(dim(pag)[1])

    ## JCI: special treatment for context variables
    if (!is.null(contextVars) && (jci == "123" || jci == "12" || jci == "1")) {
      for (x in seq_len(p)) {
        for (y in seq_len(p)) {
          if (jci == "123" && x != y && x %in% contextVars && y %in% contextVars ) {
            if (verbose)
              cat("\nJCI Assumption 3:",
                  "\nOrient:", x, "*->", y, " because they are contextVars\n")
            pag[x,y]<-2
          } 
          if (pag[x,y] && (x %in% contextVars) && !(y %in% contextVars)) {
            if( jci == "123" || jci == "12" ) {
              if (verbose)
                cat("\nJCI Assumptions 1,2:",
                    "\nOrient:", x, "-->", y, " because ", x, " is a context variable, whereas ", y, " is not.\n")
              pag[y,x]<-3
              pag[x,y]<-2
            } else if( jci == "1" ) {
              if (verbose)
                cat("\nJCI Assumption 1:",
                    "\nOrient:", x, "*->", y, " because ", x, " is a context variable, whereas ", y, " is not.\n")
              pag[x,y]<-2
            }
          }
        }
      }
    }

    ## orient collider
    if (orientCollider) {
      ind <- which(pag == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        allZ <- setdiff(which(pag[y, ] != 0), x)
        if( !(jci == "123" && y %in% contextVars) ) { ## JCI: special treatment for context variables
          for (z in allZ) {
            if (pag[x, z] == 0 && !((y %in% sepset[[x]][[z]]) ||
                     (y %in% sepset[[z]][[x]]))) {
              if (length(unfVect) == 0) {
                if (verbose) {
                  cat("\n", x, "*->", y, "<-*", z, "\n")
                  cat("Sxz=", sepset[[z]][[x]], "and",
                      "Szx=", sepset[[x]][[z]], "\n")
                }
                pag[x, y] <- pag[z, y] <- 2
              }
              else {
                if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                    !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
                  if (verbose) {
                    cat("\n", x, "*->", y, "<-*", z, "\n")
                    cat("Sxz=", sepset[[z]][[x]], "and",
                        "Szx=", sepset[[x]][[z]], "\n")
                  }
                  pag[x, y] <- pag[z, y] <- 2
                }
              }
            }
          }
        }
      }
    } ## end: Orient collider

    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R1--R4
      old_pag1 <- pag
      ##-- R1 ------------------------------------------------------------------
      if (rules[1]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2] # pag[a,b] == 2, pag[b,a] != 0
          indC <- which((pag[b, ] != 0 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a) # pac[b,C] != 0, pag[C,b] ==1, pag[a,C] == 0, pag[C,a] == 0
          if (length(indC) > 0) {
            if (length(unfVect) == 0) {
              if( any(pag[b, indC] == 3) && verbose )
                cat('Contradiction in Rule 1!\n')
              pag[b, indC] <- 2
              pag[indC, b] <- 3
              if (verbose)
                cat("\nRule 1",
                    "\nOrient:", a, "*->", b, "o-*", indC,
                    "as:", b, "->", indC, "\n")
            }
            else {
              for (c in indC) {
                if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                    !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                  if( pag[b, c] == 3 && verbose )
                    cat('Contradiction in Rule 1\n')
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                  if (verbose)
                    cat("\nRule 1",
                        "\nConservatively orient:", a, "*->", b, "o-*",
                        c, "as:", b, "->", c, "\n")
                }
              } ## for( c )
            }
          }
        } ## for( i )
      }
      ##-- R2 ------------------------------------------------------------------
      if (rules[2]) {
        ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2] # pag[a,c] == 1, pag[c,a] != 0
          indB <- which((pag[a, ] == 2 & pag[, a] == 3 & pag[c, ] != 0 & pag[, c] == 2) | (pag[a, ] == 2 & pag[, a] != 0 & pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) > 0) {
            pag[a, c] <- 2
            if (verbose) {
              cat("\nRule 2","\n")
              cat("Orient:", a, "->", indB, "*->", c, "or", a, "*->", indB, "->", c, "with", a, "*-o", c, "as:", a, "*->", c, "\n")
            }
          }
        }
      }
      ##-- R3 ------------------------------------------------------------------
      if (rules[3]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          d <- ind[i, 2] # pag[b,d] != 0, pag[d,b] == 1
          indAC <- which((pag[b, ] != 0 & pag[, b] == 2) & (pag[, d] == 1 & pag[d, ] != 0))
          if (length(indAC) >= 2) {
            if (length(unfVect) == 0) {
              counter <- 0
              while ((counter < (length(indAC) - 1)) && (pag[d, b] != 2)) {
                counter <- counter + 1
                ii <- counter
                while (ii < length(indAC) && pag[d, b] != 2) {
                  ii <- ii + 1
                  if (pag[indAC[counter], indAC[ii]] == 0 && pag[indAC[ii], indAC[counter]] == 0) {
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Orient:", d, "*->", b, "\n")
                    }
                    pag[d, b] <- 2
                  }
                }
              }
            }
            else {
              comb.indAC <- combn(indAC, 2)
              for (j in seq_len(dim(comb.indAC)[2])) {
                a <- comb.indAC[1, j]
                c <- comb.indAC[2, j]
                if (pag[a, c] == 0 && pag[c, a] == 0 && c != a) {
                  if (!any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                    pag[d, b] <- 2
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Conservatively orient:", d, "*->", b, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      ##-- R4 ------------------------------------------------------------------
      if (rules[4]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2] # pag[b,c] != 0, pag[c,b] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all a s.t. a -> c and a <-* b
          indA <- which((pag[b, ] == 2 & pag[, b] != 0) &
                        (pag[c, ] == 3 & pag[, c] == 2))
          # pag[b,a] == 2, pag[a,b] != 0
          # pag[c,a] == 3, pag[a,c] == 2
          ## chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while (length(indA) > 0 && pag[c,b] == 1) {
            a <- indA[1]
            indA <- indA[-1]
            ## path is the initial triangle
            ## abc <- c(a, b, c)
            ## Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
### MM: FIXME?? Isn't  Done  set to TRUE in *any* case inside the following
### while(.), the very first time already ??????????
            while (!Done && pag[a,b] != 0 && pag[a,c] != 0 && pag[b,c] != 0) {
              ## find a minimal discriminating path for a,b,c
              md.path <- minDiscrPath(pag, a,b,c, verbose = verbose)
              ## if the path doesn't exists, we are done with this triangle
              if ((N.md <- length(md.path)) == 1) {
                Done <- TRUE
              }
              else {
                ## a path exists
                ## if b is in sepset
                if ((b %in% sepset[[md.path[1]]][[md.path[N.md]]]) ||
                    (b %in% sepset[[md.path[N.md]]][[md.path[1]]])) {
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is in Sepset of",
                        c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                }
                else {
                  ## if b is not in sepset
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is not in Sepset of",
                        c, "and", md.path[1], ". Orient:", a, "<->", b, "<->",
                        c, "\n")
                  pag[b,c] <- pag[c,b] <- 2
                  if( pag[a,b] == 3 ) { # contradiction with earlier orientation!
                    if( verbose )
                      cat('\nContradiction in Rule 4b!\n')
                    if( jci == "0" ) {
                      pag[a,b] <- 2 # backwards compatibility
                    }
                  } else { # no contradiction
                    pag[a,b] <- 2
                  }
                }
                Done <- TRUE
              }
            }
          }
        }
      }
    } ## R1-R4
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R5-R7
      old_pag1 <- pag
      ##-- R5 ------------------------------------------------------------------
      if (rules[5]) {
        ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE) ## a o-o b
        while (length(ind) > 0) {
          a <- ind[1, 1]
          b <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all c s.t. a o-o c and c is not connected to b
          indC <- which((pag[a, ] == 1 & pag[, a] == 1) & (pag[b, ] == 0 & pag[, b] == 0))
          ## delete b since it is surely in indC
          indC <- setdiff(indC, b)
          ## find all d s.t. b o-o d and d is not connected to a
          indD <- which((pag[b, ] == 1 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          ## delete a since it is surely in indD
          indD <- setdiff(indD, a)
          if (length(indC) > 0 && length(indD) > 0) {
            counterC <- 0
            while ((counterC < length(indC)) && pag[a, b] == 1) {
              counterC <- counterC + 1
              c <- indC[counterC]
              counterD <- 0
              while ((counterD < length(indD)) && pag[a, b] == 1) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if (pag[c, d] == 1 && pag[d, c] == 1) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[a, b] <- pag[b, a] <- 3
                    pag[a, c] <- pag[c, a] <- 3
                    pag[c, d] <- pag[d, c] <- 3
                    pag[d, b] <- pag[b, d] <- 3
                    if (verbose)
                      cat("\nRule 5",
                          "\nThere exists an uncovered circle path between", a, "and", b,
                          ". Orient:", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                  }
                  else { ## conservative: check that every triple on the circle is faithful
                    path2check <- c(a,c,d,b)
                    if (faith.check(path2check, unfVect, p)) {
                      pag[a, b] <- pag[b, a] <- 3
                      pag[a, c] <- pag[c, a] <- 3
                      pag[c, d] <- pag[d, c] <- 3
                      pag[d, b] <- pag[b, d] <- 3
                      if (verbose)
                        cat("\nRule 5",
                            "\nThere exists a faithful uncovered circle path between",
                            a, "and", b, ". Conservatively orient:",
                            a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                    }
                  }
                }
                ## search with a breitensuche a minimal uncovered circle path
                else {
                  ## Find a minimal uncovered circle path for these a,b,c, and d.
                  ## This path has already been checked to be uncovered and
                  ## to be faithful for the conservative case
                  ucp <- minUncovCircPath(p, pag = pag, path = c(a,c,d,b),
                  unfVect = unfVect, verbose = verbose)
                  ## there is a path ---> orient
                  if (length(ucp) > 1) {
                    ## orient every edge on the path as --
                    n <- length(ucp)
                    pag[ucp[1], ucp[n]] <- pag[ucp[n], ucp[1]] <- 3 ## a--b
                    for (j in seq_len((length(ucp)-1))) ## each edge on the path --
                      pag[ucp[j], ucp[j + 1]] <- pag[ucp[j + 1], ucp[j]] <- 3
                  }
                }
              }
            }
          }
        }
      }
      ##-- R6 ------------------------------------------------------------------
      if (rules[6]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          if (any(pag[b, ] == 3 & pag[, b] == 3)) {
            pag[c, b] <- 3
            if (verbose)
              cat("\nRule 6",
                  "\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
          }
        }
      }
      ##-- R7 ------------------------------------------------------------------
      if (rules[7]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          indA <- which((pag[b, ] == 3 & pag[, b] == 1) & (pag[c, ] == 0 & pag[, c] == 0))
          indA <- setdiff(indA, c)
          if (length(indA) > 0) {
            if (length(unfVect) == 0) {
              pag[c, b] <- 3
              if (verbose)
                cat("\nRule 7",
                    "\nOrient:", indA, "-o", b, "o-*",
                    c, "as", b, "-*", c, "\n")
            }
            else for (a in indA)
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pag[c, b] <- 3
                if (verbose)
                  cat("\nRule 7",
                      "\nConservatively orient:", a, "-o", b, "o-*",
                      c, "as", b, "-*", c, "\n")
              }
          }
        }
      }
    } ## R5-R7
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) { ## R8-R10
      old_pag1 <- pag
      ##-- R8 ------------------------------------------------------------------
      if (rules[8]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          # pag[a,c] == 2, pag[c,a] == 1
          indB <- which(pag[, a] == 3 & (pag[a, ] == 2 | pag[a, ] == 1) &
                        pag[c, ] == 3 & pag[, c] == 2)
          if (length(indB) > 0) {
            pag[c, a] <- 3
            if (verbose)
              cat("\nRule 8",
                  "\nOrient:", a, "->", indB, "->", c,
                  "or", a, "-o", indB, "->", c, "with", a,
                  "o->", c, "as", a, "->", c, "\n")
          }
        }
      }
      ##-- R9 ------------------------------------------------------------------
      if (rules[9]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          # pag[a,c] == 2, pag[c,a] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. a (o-)--(o>) b and b and c are not connected
          indB <- which((pag[a, ] == 2 | pag[a, ] == 1) &
                        (pag[, a] == 1 | pag[, a] == 3) &
                        (pag[c, ] == 0 & pag[, c] == 0))
          ## delete c from indB since it is surely inside
          indB <- setdiff(indB, c)
          ## chose one b s.t. the initial structure exists and the edge hasn't been oriented yet
          while ((length(indB) > 0) && (pag[c,a] == 1)) {
            b <- indB[1]
            indB <- indB[-1]
            ## find a minimal uncovered pd path from initial (a,b,c) :
            upd <- minUncovPdPath(p, pag, a,b,c,
                unfVect = unfVect, verbose = FALSE)
            ## there is a path ---> orient it
            if (length(upd) > 1) {
              pag[c, a] <- 3
              if (verbose)
                cat("\nRule 9",
                    "\nThere exists an uncovered potentially directed path between", a, "and", c,
                    ". Orient:", a, " ->",c, "\n")
            }
          }
        }
      }
      ##-- R10 ------------------------------------------------------------------
      if (rules[10]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          # pag[a,c] == 2, pac[c,a] == 1
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. b --> c
          indB <- which((pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (pag[c, a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (pag[c, a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if ((pag[a, b] == 1 || pag[a, b] == 2) &&
                    (pag[b, a] == 1 || pag[b, a] == 3) &&
                    (pag[a, d] == 1 || pag[a, d] == 2) &&
                    (pag[d, a] == 1 || pag[d, a] == 3) && pag[d, b] == 0 && pag[b, d] == 0) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[c, a] <- 3
                    if (verbose)
                      cat("\nRule 10 [easy]",
                          "\nOrient:", a, "->", c, "\n")
                  }
                  else ## conservative version: check faithfulness of b-a-d
                    if (!any(unfVect == triple2numb(p,b,a,d), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,d,a,b), na.rm = TRUE)) {
                      pag[c, a] <- 3
                      if (verbose)
                        cat("\nRule 10 [easy]",
                            "\nConservatively orient:", a, "->", c, "\n")
                    }
                }
                ## search with a breitensuche two minimal uncovered circle paths
                else {
                  ## find all x s.t. a (o-)--(o>) x
                  indX <- which((pag[a, ] == 1 | pag[a, ] == 2) &
                                (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    counterX1 <- 0
                    while (counterX1 < length(indX) && pag[c, a] == 1) {
                      counterX1 <- counterX1 + 1
                      first.pos <- indX[counterX1]
                      indX2 <- setdiff(indX, first.pos)
                      counterX2 <- 0
                      while (counterX2 < length(indX2) && pag[c, a] == 1) {
                        counterX2 <- counterX2 + 1
                        sec.pos <- indX2[counterX2]
                        if (first.pos == b) { t1 <- c(a,b)
                  } else {
                    t1 <- minUncovPdPath(p, pag, a, first.pos, b, unfVect = unfVect, verbose = verbose) ## is first X mu? meaning is there a path from a to b through X?
                  }
                  if (length(t1) > 1) {     
                    if (sec.pos == d) { t2 <- c(a,d) 
                    } else { 
                      t2 <- minUncovPdPath(p, pag, a, sec.pos, d, unfVect = unfVect, verbose = verbose)
                    }
                        ## start cut
                        ## t1 <- minUncovPdPath(p, pag, a, first.pos, b,
                           ##                  unfVect = unfVect, verbose = verbose)
                        ## if (length(t1) > 1) { # otherwise, can skip next minUnc..()
                          ## t2 <- minUncovPdPath(p, pag, a, sec.pos, d,
                             ##                  unfVect = unfVect, verbose = verbose)
                          ## end cut
                          if (length(t2) > 1 &&
                              first.pos != sec.pos && pag[first.pos, sec.pos] == 0) {
                            ## we found 2 uncovered pd paths
                            if (length(unfVect) == 0) { ## normal version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10", "\nOrient:", a, "->", c, "\n")
                            }
                            else if(!any(unfVect == triple2numb(p,first.pos, a, sec.pos), na.rm = TRUE) &&
                                    !any(unfVect == triple2numb(p,sec.pos, a, first.pos), na.rm = TRUE)) {
                              ## conservative version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10",
                                    "\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                      } #  # while ( counterX2 .. )
                    }
                  }
                } # else
              } # while ( counterD .. )
            } # while ( counterB .. )
          } # if (length(indB) .)
        }
      } ## if (rules[10] ..)
    } ## R8-R10
  }
  pag
} ## udag2pag()

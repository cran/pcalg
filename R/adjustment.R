adjustment <- function(amat, amat.type, x, y, set.type)
{
  ## Purpose: Compute adjustment sets
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## amat: input adjacency matrix encoded as in ?amatType
  ## amat.type: type of adjacency matrix (dag, cpdag, pdag, mag, pag)
  ## x: col number of exposure variable
  ## y: col number of outcome variable
  ## set.type: min/all (can?)
  ## No direct effect
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Apr 2017, 09:47
  stopifnot(set.type %in% c("minimal", "canonical", "all"))
  stopifnot(ncol(amat) == nrow(amat))
  stopifnot(is.numeric(x), is.numeric(y),
            x > 0, x <= ncol(amat), y > 0, y <= ncol(amat))
  stopifnot(amat.type %in% c("dag", "cpdag", "pdag", "mag", "pag"))
  
  if (amat.type %in% c("dag", "cpdag", "mag", "pag")) {
    ## convert to dagitty object
    if (length(colnames(amat)) > 0) {
      lb <- colnames(amat) ## use existing labels
    } else {
      lb <- as.character(seq(1, ncol(amat)))
    }
    x.lb <- lb[x]
    y.lb <- lb[y]
    cpdag.dagitty <- pcalg2dagitty(amat = amat,
                                   labels = lb,
                                   type = amat.type)
    
    ## compute adjustment sets in dagitty
    sets <- adjustmentSets(cpdag.dagitty,
                           exposure = x.lb, outcome = y.lb,
                           type= set.type, effect="total")
    
    ## convert to conventional list of character vectors
    ## (potentially of length zero)
    tmp <- unclass(sets) ## vectors of node labels
    
    ## if any adjustment sets exist, convert to col positions
    if (length(tmp) > 0) {
      res <- lapply(tmp, function(x) which(lb %in% x))
    } else {
      res <- tmp ## empty list if not adj set
    }
  } else { ## type pdag
    ## step 1
    isA <- isAmenable(m = amat, x = x, y = y, type = amat.type)
    if (!isA) {
      message("Graph is not amenable\n")
      return(vector(mode = "list", length = 0))
    }
    f <- bforbiddenNodes(m = amat, x = x, y = y)
        ##############################################################
    ### EMA : The changes start in the folowing line, added if (set.type == "canonical") {new code } else { old code }
    if (set.type == "canonical"){
      
      ##form canonical set adjustb
      possanx <- possany <- adjustb <-  c()
      for(i in 1:max(length(x),length(y)))
      {
        if (i <= length(x)){
          possanx <- union(possanx,possAn(m = amat, x = x[i], type = amat.type))
        }
        if (i <= length(y)){
          possany <- union(possany,possAn(m = amat, x = y[i], type = amat.type))
        }
      }
      adjustb <- union(possanx,possany)
      adjustb <- setdiff(adjustb,union(f,union(x,y)))
      adjustb <- sort(adjustb)
      
      result.gac <- gac(amat = amat, x = x, y = y, z = adjustb, type = amat.type)$gac
      
      ##output the canonical set if it is valid
      if (result.gac)
      {
        nn <- length(adjustb)
        if (nn > 0) {
            res <- list('1' = adjustb)
        } else {
            ## res <- vector("list", 0L)
            res <- list('1' = integer(0))
        }
      } else {
        res <-  vector(mode = "list", length = 0)
      }   
    }      ### EMA : changes end here, in terms of the if statement 
    else { ###  EMA : new line
    ## step 2
    oneDag <- pdag2dag(as(t(amat),"graphNEL"))    
    amatD <- t(as(oneDag$graph,"matrix"))
    lb <- colnames(amatD)
    x.lb <- lb[x]
    y.lb <- lb[y]
    
    ## step 3
    daggityD <- pcalg2dagitty(amat = amatD,
                              labels = lb,
                              type = "dag")
    
    ## step 4
    sets <- unclass(adjustmentSets(x = daggityD,
                                   exposure = x.lb, outcome = y.lb,
                                   type= set.type, effect="total") )
    
    ## step 5: convert node labels to col positions
    nn <- length(sets)
    if (nn > 0) {
      tmp <- lapply(sets, function(x) which(lb %in% x))
    } else {
      tmp <- sets ## empty list if not adj set
    }
    
    ## step 6:
    ## output all sets in tmp that do not contain forbidden nodes
    if (nn > 0) { ## at least one adj set
      idxSel <- sapply(tmp, function(x) any(x %in% f))
      res <- tmp[!idxSel] ## only sets WITHOUT forbidden nodes
    } else {
      res <- vector("list", 0L)
    }
    } ### EMA : new line
  }
  res
}

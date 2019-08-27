optAdjSet <- function(graphEst, x.pos, y.pos){
  ## graphEst    GraphNel object
  ## x.pos  Node positions of treatment nodes.
  ## y.pos  Node positions of outome nodes.
  
  if (inherits(graphEst, "graphNEL")){
     m <- ad.g <- wgtMatrix(graphEst)
     m[which(m != 0)] <- 1
  } else if (is.matrix(graphEst)) m <- graphEst
  
  # stopifnot(x.pos==(x <- as.integer(x.pos)), y.pos==(y <- as.integer(y.pos)), length(x)==1, length(y)==1)
  
    if (!isValidGraph(amat=m, type="pdag")) {
      stop("The input graph is not a valid DAG, CPDAG or maximal PDAG. See function isValidGraph() for details.\n")
    }
  

  
  possDeX <- unlist(lapply(x.pos,possDe,m=m,y=c(),type="pdag"))
  
  noDescY <- unlist(lapply(y.pos, function(y) length(intersect(y,possDeX))))
  
  if(sum(noDescY)==0){
    stop("\nNo node in Y is a descendant of X \n")
  }
  
  y.pruned <- y.pos[which(noDescY!=0)]
  
  if(sum(noDescY) < length(y.pos)){
    cat("\n\nNot all nodes in Y are descendants of the nodes in X \n")
    cat("The algorithm will continue with the pruned Y=",y.pruned,".\n") 
  }

  
  possAnY <- unlist(lapply(y.pruned,possAn,m=m,
                           y=x.pos,type="pdag"))
  
  cn <- intersect(possAnY,possDeX)
  forb <- unlist(lapply(cn,possDe,m=m,
                        y=c(),type="pdag"))
  
  if(length(intersect(forb,x.pos))!=0 | !isAmenable(m, x = x.pos, y = y.pruned, type = "pdag")) {
    stop("There is no valid adjustment set relative to (X,Y) in the given graph.\n")
  }
  
  
  Oset <- unlist(lapply(cn,function(x,m) setdiff(which(m[x,]==1),which(m[,x]==1)),m=m))
  Oset <- setdiff(Oset, c(forb,x.pos))
  return(Oset)
}



isAmenable <- function(m,x,y, type = "pag") {
  ## INPUT: adj.matrix m; sets of node positions x and y; type in DAG, CPDAG,
  ## MAG or PAG
  ## OUTPUT: TRUE if m is amenabel wrt x,y; o/w FALSE
  found <- FALSE ## if found == TRUE at any time, graph is not amenable wrt x,y
  ## DAG is always amenable
  if (type %in% c("pdag", "cpdag","dag","pag","mag")) { ##changed added dag just in case, makes no difference
    if (type == "dag")  ## added the if case for dags
      return(!found)
    i <- 0
    p <- length(x)
    
    ## for all nodes in x, if amenability is still possible
    while ( (i<p) & !found) {
      i <- i+1
      ## posDesc of x[i] without going through any other x node
      ## posDesc <- possibleDeProper(m,x[i],x[-i])
      posDesc <- possDe(m = m, x = x[i], y = x[-i], possible = TRUE,
                        ds = FALSE, type = type)
      ## potential problem for amenability only if there is a
      ## pdp from x[i] to y
      if ( length(intersect(y, posDesc)) != 0 ) {
        nb <- as.vector(which(m[x[i],]!=0 | m[,x[i]]!=0)) ## nbrs of x[i]
        ## potentially first node on pdp from x[i] to y; however, not yet sure
        cand <- intersect(nb, posDesc)
        j <- 0
        ## for all candidate nodes, if amenability is still possible
        ## (also covers case if cand is empty)
        while ( (j<length(cand)) & !found ) {
          j <- j+1
          ## check if there is a pdp from cand[j] to y without going through x[i]
          ## cand could already be in y
          ## pathOK <- ( length(intersect(y, possibleDeProper(m,cand[j],x[i]))) != 0 )
          pdpTemp <- possDe(m = m, x = cand[j], y = x[i],
                            possible = TRUE, ds = FALSE,
                            type = type)
          pathOK <- ( length(intersect(y, pdpTemp)) != 0 )
          if (pathOK) {
            isPDAG <- (type == "pdag" | type == "cpdag") ##changed
            PDAGproblem <- (isPDAG & (m[x[i],cand[j]] == 1)) ##changed
            PAGproblem1 <- (!isPDAG & (m[x[i], cand[j]] != 2) & (m[cand[j], x[i]] != 3)) ##changed
            isDirEdge <- ((m[x[i], cand[j]] == 2) & (m[cand[j], x[i]] == 3))
            PAGproblem2 <- (!isPDAG & isDirEdge & !visibleEdge(m, x[i], cand[j])) ##changed
            found <- (PDAGproblem | PAGproblem1 | PAGproblem2)
          } ## if pathOK
        } ## while cand
      } ## if path from x[i] to y
    } ## while x
    return(!found)
  } else { ## if graph type not known
    cat("Not a valid graph type! Should be written in lowercase! \n")
    return(NULL)
  }
}

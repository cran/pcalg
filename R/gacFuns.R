cond2node <- function(s, x,y,z,m,type) {
    ## INPUT: s in format for newStackEls
    ## OUTPUT: TRUE if pth is found; o/w FALSE

    ## Is there a path of length one ?
    tmp <- newStackEls(s, x,y,z,m,type)
    if (tmp$suc) {
        suc <- TRUE
    } else {
        ## Is there a longer path ?
        suc <- FALSE
        if (length(tmp$res) > 0) {
            ## nbrs to consider
            i <- 0
            while ( (!suc) & (i < length(tmp$res)) ) {
                i <- i+1
                s2 <- tmp$res[[i]]
                suc <- cond2node(s2, x,y,z,m,type)
            }
        }
    }
    suc
}
cond2 <- function(x,y,z,m,type) {
    ## OUTPUT: True if cond 2 is true; o/w false
    lx <- length(x)
    pthFound <- rep(FALSE, lx)
    for (i in 1:lx) {
        s <- list(pth=x[i], ncp=FALSE)
        ## cond2node: TRUE if path is found
        pthFound[i] <- cond2node(s,x,y,z,m,type)
    }
    ## if no path is found, cond 2 is true
    all(pthFound == FALSE)
}


defStat <- function(ln, cn, nn, m, type = "pag") {
  ## INPUT: Node positions ln, cn, nn; adj.mat m (for DAG/CPDAG/MAG/PAG)
  ## OUTPUT: TRUE if path ln-cn-nn in m is def.stat. o/w FALSE
  res <- FALSE
  stopifnot(type %in% c("dag", "cpdag", "mag", "pag"))

  if (type %in% c("mag", "pag")) {
      isCollider <- (m[ln,cn] == 2) & (m[nn,cn] == 2)
      hasTail <- (m[ln,cn] == 3) | (m[nn,cn] == 3)
      hasCircles <- (m[ln,cn] == 1) & (m[nn,cn] == 1)
      isUnshielded <- (m[ln,nn] == 0) & (m[nn,ln] == 0)
      isDefNonCollider <- ( (hasCircles & isUnshielded) | hasTail )
  } else {
      ## DAG, CPDAG
      isCollider <- (m[ln,cn]==0 & m[cn,ln]==1) & (m[cn,nn]==1 & m[nn,cn]==0)
      hasTail <- (m[ln,cn]==1 & m[cn,ln]==0) | (m[cn,nn]==0 & m[nn,cn]==1)
      isUndirected <- (m[ln,cn]==1 & m[cn,ln]==1) | (m[cn,nn]==1 & m[nn,cn]==1)
      isUnshielded <- (m[ln,nn]==0 & m[nn,ln]==0)
      isDefNonCollider <- ( (isUndirected & isUnshielded) | hasTail )
  }

  if (isCollider) {
    res <- TRUE
  } else {
    if (isDefNonCollider) res <- TRUE
  }

  res
}
desc <- function(m, possible = FALSE, type = type) {
  ## computes all (possible) descendants of every node
  ## INPUT: adj.matrix m in MAG/PAG or DAG/CPDAG coding; if possible=TRUE, possible desc are found; o/w descendants
  ## OUPUT: list of vectors; vector at list position i contains node positions of (possible) descendants of node i
  p <- ncol(m)
  res <- vector(mode = "list", length = p)
  for (i in 1:p) {
      ## res[[i]] <- possibleDeProper(m=m, x=i, y=NULL, possible=possible)
      res[[i]] <- possDe(m = m, x = i, y = NULL,
                               possible = possible,
                               ds = FALSE, type = type)
  }
  res
}
forbiddenNodes <- function(m,x,y,type)
{
    ## INPUT: adj.matrix in DAG,CPDAG,MAG,PAG coding; sets of node positions x and y
    ## OUTPUT: set of node positions of nodes in the forbidden set (sorted)
  n1 <- length(x)
  n2 <- length(y)

  possDeX <- possAnY <- c()
  for(i in 1:max(n1,n2))
  {
    if (i <= n1) {
      #find all possible descendants of a node x[i] that are on a proper path
      #relative to x (exclude x[i] because descendants of x[i] that are not
      #on a proper path are allowed)
        ## possDeX <- union(possDeX, setdiff(possibleDeProper(m,x[i],x),x[i]))
        pdpTmp <- possDe(m = m, x = x[i], y = x, possible = TRUE,
                         ds = FALSE, type = type)
        possDeX <- union(possDeX, setdiff(pdpTmp, x[i]))
    }
    if (i <= n2) {
      #find all possible ancestors of a node y[i] that are
      #on a proper path relative to x
        ## possAnY <- union(possAnY, possibleAnProper(m,y[i],x))
        papTmp <- possAn(m = m, x = y[i], y = x, possible = TRUE,
                         ds = FALSE, type = type)
        possAnY <- union(possAnY, papTmp)
    }

  }

  #a set of all nodes on a proper possibly directed path from X to Y
  pdp <- intersect(possDeX,possAnY)

  #the forbiden node set are all possible descendants of nodes in pdp
  fbnodes <- c()
  if (length(pdp) > 0) {
      for(j in 1:length(pdp))
      {
          ## fbnodes <- union(fbnodes,possibleDeProper(m,pdp[j],c()))
          pdpTmp2 <- possDe(m = m, x = pdp[j], y = c(), possible = TRUE,
                            ds = FALSE, type = type)
          fbnodes <- union(fbnodes, pdpTmp2)
      }
  }

  if (length(fbnodes) > 0) {
      return(sort(fbnodes))
  } else {
      return(fbnodes)
  }

}

##gac for pdags
##
## this function is the same as its equivalent in pcalg except it calls
## the pdag fucntions given above
##
##currrently separated dags,pdags,cpdags from mags,pags but could be merged

##however the res[3] must be calculated separately for dags,pdags,cpdags and mags,pags

## uses amat.cpdag encoding amat[i,j]=0,amat[j,i]=1 <=> i ->j
gac <- function (amat, x, y, z, type = "pag")
{
    if (inherits(amat, "amat")) amat <- as(amat, "matrix")
    if (!is.null(dimnames(amat))) dimnames(amat) <- NULL
    if (type %in% c("dag", "cpdag", "pdag")) {
          if (!isValidGraph(amat = amat, type = type)) {
              message("The input graph is not a valid ",type,
                      ". See function isValidGraph() for details.")
          }
    }
  res <- c(isAmenable(amat, x = x, y = y, type = type),
           NA, NA)
  if (type %in% c("dag","pdag","cpdag")){
    f <- bforbiddenNodes(m = amat, x = x, y = y)
    res[2] <- (length(intersect(f, z)) == 0)
    res[3] <- cond3fast(x = x, y = y, z = z, m = amat)
  } else { ##the code for MAGs, PAGs is still the same
    f <- forbiddenNodes(amat, x = x, y = y)
    res[2] <- (length(intersect(f, z)) == 0)
    res[3] <- cond2(x = x, y = y, z = z, m = amat, type = type)
  }
  list(gac = all(res), res = res, f = f)
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
mcon <- function(ln,cn,nn,m,z,descList) {
  ## INPUT: node positions of last, current and next node on path; adj.mat. m
    ## in DAG/CPDAGMAG/PAG format;
  ## set z; descList as returned from desc(m)
  ## OUTPUT: TRUE if ln and nn are mcon given z on path ln-cn-nn
  res <- FALSE
  isColliderPAG <- ( (m[ln,cn] == 2) & (m[nn,cn]==2) )
  isColliderDAG <- ( (m[ln,cn]==0 & m[cn,ln]==1) & (m[cn,nn]==1 & m[nn,cn]==0) )
  isCollider <- ( isColliderPAG | isColliderDAG)
  if (isCollider) {
    if (any(descList[[cn]] %in% z)) res <- TRUE
  } else {
    if ( !(cn %in% z) ) res <- TRUE
  }
  res
}
ncEdge <- function(cn, nn, m, type = "pag") {
    stopifnot(type %in% c("dag", "cpdag", "mag", "pag"))
    res <- FALSE
    if (type %in% c("mag", "pag")) {
        if (m[nn,cn]==2) res <- TRUE
    } else { ## must be DAG or CPDAG
        if ( (m[cn,nn] == 1) & (m[nn,cn] == 0)) res <- TRUE
    }
    res
}
newStackEls <- function(s,x,y,z,m,type) {
    ## INPUT: previous stack element of type list: {pth (num. vec.), ncp (bool);
    ## ncp is TRUE, if pth is non-causal; o/w FALSE
    ## OUTPUT: List 'res' of new stack elements (poss. empty); suc = TRUE, if an
    ## open path to y was found
    suc <- FALSE
    pth <- s$pth; ncp <- s$ncp
    lp <- length(pth)
    if (lp == 1) {
        ## Starting node on search
        cn <- pth
        nb <- setdiff(as.vector(which(m[cn,]!=0 | m[,cn]!=0)), x)
        ## nb does not contain x -> proper path
        res <- vector("list", length(nb))
        i <- 0
        while ( (!suc) & (i < length(nb)) ) {
            ## Search through possible neighbors
            i <- i + 1
            nn <- nb[i]
            ## m-con, def.stat. are for free -> check non-causal path
            ncpTmp <- ( ncp | ncEdge(cn,nn,m,type=type) )
            ## suc=T: Exists proper, m-con, d.s., non-c path from x to y
            suc <- ( (nn %in% y) & ncpTmp )
            res[[i]] <- list(pth = c(cn,nn), ncp = ncpTmp)
        }
    } else {
        ## Path has at least 2 nodes
        descList <- desc(m, type = type)
        ln <- pth[lp-1]; cn <- pth[lp]
        nb <- setdiff(setdiff(as.vector(which(m[cn,]!=0 | m[,cn]!=0)), x), pth)
        res <- list()
        j <- 0
        jj <- 0
        while( (!suc) & (j<length(nb)) ) {
            j <- j+1
            nn <- nb[j]
            mc <- mcon(ln=ln, cn=cn, nn=nn, m=m, z=z, descList=descList)
            ds <- defStat(ln=ln, cn=cn, nn=nn, m=m, type = type)
            ncpTmp <- ( ncp | ncEdge(cn,nn,m,type=type) )
            if (mc & ds) {
                jj <- jj + 1
                ## proper, m-connecting & def.stat. path
                if ( (nn %in% y) & ncpTmp ) suc <- TRUE
                res[[jj]] <- list(pth=c(pth,nn), ncp=ncpTmp)
            }
        }
    }
    list(res=res, suc=suc)
}

## ## CPDAG
## cp <- function(m) {
##     require(Rgraphviz)
##     plot(as(t(m),"graphNEL"))
## }

## ## MAG/PAG
## mp <- function(m) {
##     colnames(m) <- rownames(m) <- as.character(1:ncol(m))
##     require(Rgraphviz)
##     plotAG(m)
## }
possibleAnProper <- function(m,x,y=NULL)
{
    ## INPUT
    ## m: Adjacency matrix (coding for MAG/PAG)
    ## x: Starting node (col.position in m)
    ## y: Set Y (col.positions in m)
    ## OUTPUT
    ## All nodes with a possibly directed path from a to x not going through y
    ## including x itself

    p <- length(m[,1]) ## nmb of nodes
    q <- v <- rep(0,p) ## queue
    ## q has col.pos. of unvisited nodes = queue
    ## v has col.pos. of visited nodes
    i <- k <-  1 ## i: end of queue; k: current point in queue
    q[1] <- x ## x is first node in queue
    tmp <- m ## ???

    while(q[k]!=0 & k<=i) {## queue is not empty & current pos is within queue
        t <- q[k] ## take new node from queue
        v[k] <- t ## mark t as visited (if t is in v, it is visited)
        k <- k+1 ## increase current position in queue
        for(j in 1:p) { ## check if j is a possible parent of t
            ## check if edgemark at j is circle (1) or tail (3) but not head (2) or empty (0)
            if (tmp[t,j] %in% c(1,3)) { ## works for CPDAG, too
                ## only add nodes that:
                ## aren't already scheduled for a visit
                ## are not in y
                if (!(j %in% q) & !(j %in% y)) {##
                    i <- i+1 ## incearse size of queue by one
                    q[i] <- j ## add node to queue
                }
            }
        }
    }
    sort(setdiff(v,c(0))) ## remove trailing zeros from initial vector
}
#finds all possible descendants of a node x in a graph
#that are on a proper path relative to nodeset Y (that is, that don't go through Y)
#m is the adjacency matrix
possibleDeProper <- function(m,x,y=NULL,possible = TRUE)
{
  ## INPUT: adj.mat. m in MAG/PAG or DAG/CPDAG coding; node pos x;
  ## set of node pos y;
  ## If possible == TRUE, possible Desc. are found; o/w descendents are found
  ## OUTPUT: Node positions of (possible) descendents (ignoring path trough y); sorted
  #q denotes unvisited nodes/ nodes in queue
  #v denotes visited nodes
  q <- v <- rep(0,length(m[,1]))
  i <- k <-  1
  q[i] <- x

  while(q[k]!=0 & k<=i)
  {
    t <- q[k]
    #mark t as visited
    v[k] <- t
    k <- k+1
    #in this for cycle
    #add all nodes that have a possibly directed
    #edge with node t and all parents of node t to queue
    for(j in 1:length(m[1,])) {
      if (possible) {
        ## find possible descendents
        collect <- ( ((m[t,j]!=0) && (m[j,t] %in% c(1,3))) || ## PAG; CPDAG t --- j
                     ((m[t,j]==0) && (m[j,t]==1)) ) ## CPDAG: t ---> j
      } else {
        ## find descendents
        collect <- ( ((m[j,t]==3) && (m[t,j]==2)) ||
                     ((m[j,t]==1) && (m[t,j]==0)) )
      }
      if (collect) {
        #only add nodes that haven't been added
        #and that are on a proper path
        if (!(j %in% q) & !(j %in% y))
        {
          i <- i+1
          q[i] <- j
        }
       }
    }
  }
  sort(setdiff(v,c(0)))
}

##The following function forms a proper back-door graph G_{XY}^{pbd}
##given the adjacency matrix, node position vectors x and y
##
gbg <-  function(m,x,y)
{
  tmp <- m
  for (i in seq_along(x)) {
    Desc <- bPossibleDeProper(m,x[i],x[-i])
    if (length(intersect(y, Desc)) != 0) {
      ch <- as.vector(which(m[x[i], ] == 0 & m[, x[i]] == 1))
      cand <- intersect(ch, Desc)
      for(j in seq_along(cand)) {
        ##pathOK <- (length(intersect(y, bpossibleDeProper(m,cand[j], x[i]))) != 0)
        ## what if x - y? MISSING CASE!
        pathOK <- (length(intersect(y, bPossibleDeProper(m,cand[j], x))) != 0) ||
          (cand[j] %in% y) # added x - y case
        if (pathOK) {
          tmp[cand[j],x[i]] <- 0
        }
      }
    }
  }
  tmp
}

## the following function tests the separation condition for gac-pdags
## it is a pseudo check for this condition, because provided
## that amenability and forbidden set are satisfied
## this will exactly check the separation condition

## however if amenability or forb are not satisfied
## this check is not necessarily accurate

## given a pdag adjacency matrix m it checks whether node position vector z
## satisfies the separation condition relative to (x,y)

## m encoding amat.cpdag as above m[i,j]=0, m[j,i]=1 <=> i -> j
cond3fast <- function(x, y, z, m)
{
  ## the idea is to use the lemma which says
  ## that if m is amenable and z satisfies the forbidden cond
  ## then z satisfies separation in the pdag iff
  ## z satisfies separation in at least one dag

  ## First find one dag D in the equivalence class represented by m
  oneDag <- pdag2dag(as(t(m),"graphNEL"))
  dagAmat <- t(as(oneDag$graph,"matrix"))

  ## We instead of checking the separation crit in m
  ## w can check the d-separation criterion in D_{XY}^{pbd}
  gb <- gbg(dagAmat,x,y)
  msep(t(gb),x,y,z)  ## uses the transposed of amat.cpdag encoding
}

## basically only added b's

## note that provided that possibleDeProper is implemented in
## pcalg as looking at only def. stat. paths
## then the already implemented forbiddenNodes can be used directly
## (with the addition of type: pdag)

## m - adjacency matrix amat.cpdag encoding m[i,j]=0, m[j,i] = 1 <=> i -> j
## x - vector of node position
## y - vector of node positions

## as it is now this function can only be used for pdags, dags, cpdags
## because the bPossibleDeProper etc. functions are only defined for these
## adjacency matrices
bforbiddenNodes <- function (m, x, y)
{
  n1 <- length(x)
  n2 <- length(y)
  possDeX <- possAnY <- c()
  for (i in 1:max(n1, n2)) {
    if (i <= n1)
      possDeX <- union(possDeX, setdiff(bPossibleDeProper(m, x[i], x), x[i])) ##changed
    if (i <= n2)
      possAnY <- union(possAnY, bPossibleAnProper(m, y[i], x)) ##changed
  }
  pdp <- intersect(possDeX, possAnY)
  fbnodes <- c()
  if (length(pdp) > 0) {
    for (j in 1:length(pdp)) {
      fbnodes <- union(fbnodes, bPossibleDeProper(m, pdp[j], c())) ##changed
    }
  }
  if (length(fbnodes) > 0) {
    return(sort(fbnodes))
  }
  else {
    return(fbnodes)
  }
}

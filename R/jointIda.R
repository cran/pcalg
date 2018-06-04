b.extract.parent.sets <- function (x.pos, amat.cpdag) 
{ 
  ## if not working with amat.cpdag type then:
  ## amat.cpdag <- t(amat.cpdag)
  
  amat.cpdag[which(amat.cpdag != 0)] <- 1  ##just in case make all the non-zero's 1
  amat.undir <- amat.cpdag * t(amat.cpdag) ## amat.undir - all undirected edges
  ## the undirected edge i - j has amat[i,j] = amat[j,i] =1
  amat.dir <- amat.cpdag - amat.undir      ## amat.dir - all directed edges
  pasets.dir <- lapply(x.pos, function(x) which(amat.dir[x,] != 0))   ## find all parents of x.pos in the PDAG
  
  ## sibs will be a vector containing all undirected edges connected with x.pos
  ## if for example: i is in x.pos and i-j is in G
  ## then add j;i to sibs 
  sibs <- c()
  for (i in 1:length(x.pos))
  {
    tmp.sib <- which(amat.undir[x.pos[i],]!=0)  ## tmp.sib contains all siblings of x.pos[i]
    if (length(tmp.sib)!=0){  ## if x.pos[i] has a sibling, then add all those sibling edges to sibs
      for (j in 1:length(tmp.sib))
      {
        sibs <- c(sibs,paste(tmp.sib[j],x.pos[i], sep = ";")) ## sibs is a vector of type a;b c;d 
        ## where b and d are in x.pos and b-a d-c are in G
        ## note that if a,b are in x.pos and a-b is in G
        ## then both a;b and b;a are in sibs
      }
    }
  }
  
  ## if there are no undirected edges connected to x.pos, that is
  ## if sibs is empty, so return pasets.dir as the only parent set
  if (length(sibs)==0){         
    pasets <- list(pasets.dir) 
    return(pasets)
  } else {   ## if sibs is not empty, find all possible joint parent sets
    
    pasets <- list()   ## this is the object we return, containing a list of all valid joint parent sets
    count <- 1         ## this counter will contain the current number (+1) of valid joint parent sets
    
    ## first check the no additional parents option
    ## meaning that it is possible to orient all sibling edges out of x.pos
    toAdd <- makeBgKnowledge(sibs)
    ## require that no two nodes in x.pos are siblings
    if (length(intersect(x.pos,toAdd$x))==0){
      ## try to orient everything out of x.pos
      if (!is.null(addBgKnowledge(amat.cpdag,x=toAdd$y,y=toAdd$x))){  
        ## if it is possible to orient everyting out of x.pos add pasets.dir to valid joint parent sets
        pasets[[count]] <- pasets.dir
        count <- count + 1
      }
    }
    
    ## now for all subsets of possible parents hat is
    ## all possible subsets of sibs union pasets.dir 
    ## check if they work
    size <- length(sibs)
    for (k in 1:size)
    {
      possPa <- combn(sibs,k)  ## form all combinations of sibs of size k
      for (r in 1:length(possPa[1,]))
      {
        s <- possPa[,r]   ## get one subset of sibs which we will try to orient into x.pos
        toAdd1 <- makeBgKnowledge(s)  ## transform it from sibs format: a;b c;d into bgKnowledge format that is
        ## into a data.frame with x=c(a,c) y=c(b,d) so that a -> b, c-> d is bgKnowledge
        
        sbar <- setdiff(sibs,s)       ## the complement of our subset s should be oriented out of x.pos
        toAdd2 <- makeBgKnowledge(sbar)
        
        ## the following 2 lines define the bg Knowledge that we try to add
        addFromThis <- c(toAdd1$x,toAdd2$y)  
        addToThis <- c(toAdd1$y,toAdd2$x)
        
        ## only try to add this bg knowledge if its consistent within itself
        ## meaning if it does not contain contradictory edges
        ## for example addFromThis =c(1,2) addToThis = c(2,1)
        check <- FALSE  ## if check is true it will indicate that the background knowledge contradicts itself
        for (i in 1: length(addFromThis)){
          if (addToThis[i] %in% addFromThis[which(addToThis == addFromThis[i])]){# %in% addToThis)]){
            check <- TRUE
          }
        }
        if (!check){  ##only try to add background knowledge that does not contradict itself
          if (!is.null(addBgKnowledge(amat.cpdag,x=addFromThis,y=addToThis))){ ## if addBgKnowledge does not return NULL
            ## then this possible parent set is valid
            pasets[[count]] <- formPaSets(x.pos,s,pasets.dir) ## add s union pasets.dir as a valid joint parent set
            count <- count + 1                                
          }
        }
      }
    }
    return(pasets)
  }
}

## the following functions takes character vector s of type s=c("1;2","3;4")
## and transforms it into a data frame containing vectors x=c(1,3) and y=c(2,4)
## so that one can add x -> y (that is 1 -> 2, 3 -> 4) as bg knowledge easier
makeBgKnowledge <- function(s)
{
  x <- y <- c()
  
  if (length(s)==0){   ## if s is empty, return empty vectors
    df.final <- data.frame(x=x, y=y)
    return(df.final)
  } else {            ## otherwise transform the charactor vector s into a 
    ## bg knowledge data frame
    for (i in 1:length(s))   
    {
      addFromTo <- as.numeric(unlist(strsplit(x = s[i],split = ";")))
      x <- c(x,addFromTo[1])
      y <- c(y,addFromTo[2])
    }
    df.final <- data.frame(x=x, y=y)
    return(df.final)
  }
}

## the following function takes as arguments:
## x.pos - a vector of node positions (e.g., x = c(a,b))
## s - sibling edges character vector   (e.g, s = c("c;a","d;b"))
## pasets.dir  - a list of parents of x.pos in the PDAG (e.g. [[1]] e [[2]] integer(0))

## the function returns a list containing the parent set s union pasets.dir 
## e.g., the list
## [[1]] e c
## [[2]] d
##this function is written so that the ending format of the parent sets 
## matches the parent sets returned by the existing extract.parent.sets in pcalg 
formPaSets <- function(x.pos, s, pasets.dir){
  
  dfPossPa <- makeBgKnowledge(s)
  AllPa <- list()
  for(i in 1:length(x.pos)){
    AllPa[[i]] <- pasets.dir[[i]]
    if (x.pos[i] %in% dfPossPa$y){    
      AllPa[[i]] <- sort(as.integer(union(AllPa[[i]], dfPossPa$x[which(dfPossPa$y[] ==x.pos[i])])))
      names(AllPa[[i]]) <-rep("", length(AllPa[[i]])) ## as.character(AllPa[[i]]) ##EMA
    }
  }
  return(AllPa)
} 

### Preetam's old function for cpdags!
##
extract.parent.sets <- function (x.pos, amat.cpdag, isCPDAG = FALSE) 
{
  amat.cpdag[which(amat.cpdag != 0)] <- 1
  amat.undir <- amat.cpdag * t(amat.cpdag)
  amat.dir <- amat.cpdag - amat.undir
  pasets.dir <- lapply(x.pos, function(x) which(amat.dir[, 
                                                         x] != 0))
  conn.comp.imp <- NULL
  x.temp <- x.pos
  while (length(x.temp) > 0) {
    comp.temp <- graph.dfs(graph = graph.adjacency(amat.undir, 
                                                   mode = "undirected"), root = x.temp[1], unreachable = FALSE)$order
    comp.temp <- comp.temp[!is.na(comp.temp)]
    x.temp <- setdiff(x.temp, comp.temp)
    conn.comp.imp <- c(conn.comp.imp, list(comp.temp))
  }
  conn.comp.imp <- lapply(conn.comp.imp, as.integer)
  chordal <- if (!isCPDAG) {
    vapply(conn.comp.imp, function(i) is.chordal(graph.adjacency(amat.undir[i, 
                                                                            i], mode = "undirected"), fillin = FALSE)$chordal, 
           NA)
  }
  else {
    rep(TRUE, length(conn.comp.imp))
  }
  all.locally.valid.parents.undir <- function(amat, x) {
    amat.V <- as.integer(rownames(amat))
    pa.dir <- pasets.dir[[x.pos == amat.V[x]]]
    paset <- list(pa.dir)
    pa <- which(amat[, x] != 0)
    if (length(pa) == 1) {
      paset <- c(paset, list(c(amat.V[pa], pa.dir)))
    }
    else {
      for (i in 1:length(pa)) {
        pa.tmp <- combn(pa, i, simplify = TRUE)
        n.comb <- ncol(pa.tmp)
        for (j in 1:n.comb) {
          pa.t <- pa.tmp[, j]
          new.coll <- if (length(pa.t) > 1) {
            tmp <- amat[pa.t, pa.t]
            diag(tmp) <- 1
            (min(tmp) == 0)
          }
          else FALSE
          if (!new.coll) 
            paset <- c(paset, list(c(amat.V[pa.t], pa.dir)))
        }
      }
    }
    return(paset)
  }
  extract.parent.sets.from.conn.comp <- function(i) {
    all.nodes <- conn.comp.imp[[i]]
    nvar <- length(all.nodes)
    if (nvar == 1) {
      pasets.comp <- list(pasets.dir[match(all.nodes, x.pos)])
    }
    else {
      conn.comp.mat <- amat.undir[all.nodes, all.nodes]
      rownames(conn.comp.mat) <- all.nodes
      x.. <- intersect(all.nodes, x.pos)
      m.x. <- match(x.., all.nodes)
      ii.x <- seq_along(x..)
      if (chordal[i] & nvar <= 12) {
        rownames(conn.comp.mat) <- colnames(conn.comp.mat) <- 1:nvar
        all.extensions <- pdag2allDags(conn.comp.mat)$dags
        pa.fun <- function(amat, j) c(all.nodes[which(amat[, 
                                                           m.x.[j]] != 0)], pasets.dir[[match(x..[j], 
                                                                                              x.pos)]])
        parent.sets.fun <- function(r) lapply(ii.x, pa.fun, 
                                              amat = matrix(all.extensions[r, ], nrow = nvar))
        pasets.comp <- lapply(1:nrow(all.extensions), 
                              parent.sets.fun)
      }
      else {
        pasets.comp <- lapply(m.x., all.locally.valid.parents.undir, 
                              amat = conn.comp.mat)
        idx <- expand.grid(lapply(ii.x, function(j) seq_along(pasets.comp[[j]])))
        pasets.comp <- lapply(1:nrow(idx), function(r) lapply(ii.x, 
                                                              function(j) pasets.comp[[j]][[idx[r, j]]]))
      }
    }
    return(pasets.comp)
  }
  all.pasets <- lapply(1:length(conn.comp.imp), extract.parent.sets.from.conn.comp)
  idx <- expand.grid(lapply(1:length(all.pasets), function(i) 1:length(all.pasets[[i]])))
  x.conn.comp <- unlist(lapply(1:length(conn.comp.imp), function(i) intersect(conn.comp.imp[[i]], 
                                                                              x.pos)))
  i.match <- match(x.pos, x.conn.comp)
  lapply(1:nrow(idx), function(i) unlist(lapply(1:length(conn.comp.imp), 
                                                function(j) all.pasets[[j]][[idx[i, j]]]), recursive = FALSE)[i.match])
}

## Main Function ----------------------------------------------------------
##
jointIda <- function(x.pos, y.pos, mcov, graphEst = NULL,
                      all.pasets = NULL, technique = c("RRC","MCD"), type = c("pdag", "cpdag", "dag"))
{
  nx <- length(x.pos)
  if (length(type) > 1) type <- "pdag"
  if (is.null(all.pasets)) {
    amat <- as(graphEst,"matrix")
    amat[which(amat != 0)] <- 1

    ## check if valid input amat
    if (!isValidGraph(amat = t(amat), type = type)) {
      message("The input graph is not a valid ",type,". See function isValidGraph() for details.\n")
    }

##########################
    ##EMA changes below
    if (type == "pdag"){
      all.pasets <- b.extract.parent.sets(x.pos,t(amat)) 
    } else { ##asumes the graph is a cpdag or dag
      all.pasets <- extract.parent.sets(x.pos,amat) 
    }
    ####################
  } else { ## check format of all.pasets :
    if(!is.list(all.pasets) || vapply(all.pasets, length, 1L) != nx)
      stop("all.pasets is not given in an appropriate format.")
  }
  
  if (length(y.pos) > 1) { ## call myself for each y in y.pos :
    lapply(y.pos, function(y) jointIda(x.pos, y, mcov=mcov,
                                       all.pasets=all.pasets, technique=technique))
  } else {
    if (is.element(y.pos, x.pos))
      matrix(0, nrow = nx, ncol = length(all.pasets))
    else { ## return joint.effects
      technique <- match.arg(technique)
      switch(technique,
             RRC = matrix(unlist(lapply(all.pasets,function(pasets) RRC(mcov,x.pos,y.pos,pasets))),
                          nrow = nx),
             MCD = matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                          nrow = nx))
    }
  }
}

##' MCD :=  Modifying the Cholesky Decomposition
MCD <- function(cov.mat, intervention.set, var.y, pasets, return.modified.cov.mat = FALSE) {
  if (is.element(var.y,intervention.set) & !return.modified.cov.mat)
    return(rep(0,length(intervention.set)))
  if (length(intervention.set) == 1 & is.element(var.y,unlist(pasets)) & !return.modified.cov.mat)
    return(0)
  
  if (!return.modified.cov.mat) {
    imp.var <- unique(c(intervention.set,unlist(pasets),var.y))
    if (length(imp.var) < nrow(cov.mat)) {
      cov.mat <- cov.mat[imp.var,imp.var]
      intervention.set <- match(intervention.set,imp.var)
      var.y <- match(var.y,imp.var)
      pasets <- if(length(intervention.set) > 1)
        lapply(pasets, function(x) match(unlist(x), imp.var))
      else match(unlist(pasets), imp.var)
    }
  }
  
  do.Cholesky.modification <- function(x) {
    pa.x <- if(length(intervention.set) > 1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if (length(pa.x) == 0) return(cov.mat)
    ind <- c(pa.x,x,(1:nrow(cov.mat))[-c(pa.x,x)])
    cov.mat <- cov.mat[ind,ind]
    x <- match(x,ind)
    pa.x <- match(pa.x,ind)
    temp <- gchol(cov.mat)
    Lower.tri.mat <- solve(as.matrix(temp))
    ## tmp1 <- bdsmatrix::diag(temp)  ## bug in diag function ?
    #################################################
    ## temp solution until bug in bdsmatrix is fixed
    stopifnot(length(temp@Dim)==2, temp@Dim[1] == temp@Dim[2]) ## double check input for temp solution
    tmpSeq <- seq(from = 1, to = prod(temp@Dim), by = (temp@Dim[1]+1)) ## index of diagonal
    tmp1 <- temp@.Data[tmpSeq]
    ##################################################
    Diag.mat <- base::diag(tmp1)
    Lower.tri.mat[x,pa.x] <- 0
    ## MM{FIXME !} :
    cov.mat <- solve(Lower.tri.mat) %*% Diag.mat %*% t(solve(Lower.tri.mat))
    ## return
    cov.mat[order(ind),order(ind)]
  }
  
  for (i in 1:length(intervention.set)) {
    cov.mat <- do.Cholesky.modification(intervention.set[i])
  }
  
  if(return.modified.cov.mat)
    cov.mat
  else {
    MCD.estimate <- function(x) {
      if (is.element(var.y, unlist(pasets[match(x,intervention.set)])))
        0
      else
        cov.mat[var.y,x]/cov.mat[x,x]
    }
    vapply(intervention.set, MCD.estimate, double(1))
  }
} ## {MCD}


##' RRC := Recursive Regressions for Causal effects
RRC <- function(cov.mat, intervention.set, var.y, pasets) {
  
  adjusted.regression <- function(x,y) {
    if (x == y) return(0)
    pa.x <- if(length(intervention.set) > 1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if(is.element(y,pa.x)) 0 else
      solve(cov.mat[c(x,pa.x), c(x,pa.x)],
            cov.mat[c(x,pa.x), y])[1]
  }
  
  ## Define the vector of causal effects of intervention variables on var.y
  ## ISo   := intervention.set.on
  ## ISoIS := ISo.intervention.set := intervention.set.on.intervention.set
  ISo.var.y <- sapply(intervention.set, adjusted.regression, y = var.y)
  
  ## Define the required matrix of single intervention effects
  if (length(intervention.set) > 1) {
    ISoIS <-
      matrix(apply(expand.grid(intervention.set, intervention.set),
                   1L, function(x) adjusted.regression(x[1],x[2])),
             nrow = length(intervention.set))
  } else {
    return(ISo.var.y)
  }
  
  joint.effect.fun  <- function(x) {
    if(is.element(var.y, unlist(pasets[match(x, intervention.set)])))
      return(0)
    x.temp <- match(x, intervention.set)
    ## Intervention Set without x --> I.S. w/o x  --> IS.wo.x
    IS.wo.x <- intervention.set[-x.temp]
    ## Intitialize the RR estimate as the single intervention effect of intervention.set on var.y
    RR.estimate <- ISo.var.y[x.temp]
    ## Define the vector of causal effects of intervention.set on other intervention variables
    x.t.oIS <- ISoIS[x.temp,-x.temp]
    ## Define the vector of causal effects of "other" intervention variables on var.y
    ISo.var.y.temp <- ISo.var.y[-x.temp]
    ## Define the required matrix of single intervention effects
    ISoIS.temp <- ISoIS[-x.temp,-x.temp]
    
    while(length(IS.wo.x) > 1) {
      ## update RR.estimate and the other things accounting for
      ## the elimination of the first entry of the current intervention set
      RR.estimate <- RR.estimate - x.t.oIS[1]*ISo.var.y.temp[1]
      ISo.var.y.temp <- ISo.var.y.temp[-1] - ISoIS.temp[-1,1] * ISo.var.y.temp[1]
      x.t.oIS <- x.t.oIS[-1] - x.t.oIS[1]*ISoIS.temp[1,-1]
      if (length(IS.wo.x) > 2)
        ISoIS.temp <- ISoIS.temp[-1,-1] - tcrossprod(ISoIS.temp[-1,1],
                                                     ISoIS.temp[1,-1])
      IS.wo.x <- IS.wo.x[-1]
    }
    ## return
    RR.estimate - x.t.oIS * ISo.var.y.temp
  }
  
  sapply(intervention.set, joint.effect.fun)
}## {RRC}


###  MM: (ess-set-style 'DEFAULT) : we have much nesting ==> only indent by 2
## Local Variables:
## eval: (ess-set-style 'DEFAULT 'quiet)
## delete-old-versions: never
## End:

## the following function returns TRUE if there are no partially directed or directed cycles in G
## where G is a pdag (cpdag) represented by it's adjacency matrix: mat (which is in the form of amat.cpdag)
noCycles <- function(mat)
{
  ok <- TRUE
  ## for every node i in G
  for(i in 1: length(mat[1,]))
  {
    ##parents of i in the dag/cpdag/pdag if any
    pa.i <- which( mat[i,]!=0 & mat[,i]==0)
    
    ##now find all possible ancestors of the parents
    ## if there are any parents
    if (length(pa.i) !=0){
      for (j in 1: length(pa.i)){
          pos.anc <-  setdiff(possAn(m=mat, x=pa.i[j], ds=FALSE, possible=TRUE),
                              pa.i[j])
          if (i %in% pos.anc)    ## if i is a possible ancestor of a parent of i in G
              ## then there is a partially directed (or directed)
              ## cycle in G
              ok <- FALSE
      }
    }
  }
  ok
}


####################################

##the following function is adapted from the code of jointIda in pcalg
## it uses some functions of the igraph package
## given a pdag (cpdag) this function finds all the connected components of the undirected part
## and checks if each of those is chordal
## it returns a list (of connected components), where each element of the list (are two sublists):
## [[1]] the elements of a connected component [[2]] TRUE/FALSE: depending on whether this component is chordal
chordalComponents <- function(amat.cpdag){
  
  amat.undir <- amat.cpdag * t(amat.cpdag)    ##get the undirected part of the graph

  conn.comp.imp <- NULL   ## conn.comp.imp will contain all connected undirected components as a list

  x.temp <- which(apply(amat.undir,1,sum)!=0)  ## which nodes in the undirect. part are not singletons
  
  #the following while loop finds all connected components of the undirected part of the graph  
  while (length(x.temp) > 0) {
    comp.temp <- dfs(graph = graph_from_adjacency_matrix(amat.undir, mode = "undirected"), root = x.temp[1], unreachable = FALSE)$order
    comp.temp <- comp.temp[!is.na(comp.temp)]
    x.temp <- setdiff(x.temp, comp.temp)
    conn.comp.imp <- c(conn.comp.imp, list(comp.temp))
  }
  conn.comp.imp <- lapply(conn.comp.imp, as.integer, conn.comp.imp)
  
  ## chordal is a vector of booleans, the length of the number of connected components in conn.comp.imp
  ## TRUE or FALSE value depends on whether each conn comp is chordal
  chordal <-  vapply(conn.comp.imp, function(i) is_chordal(graph_from_adjacency_matrix(amat.undir[i,i], mode = "undirected"), 
                                                           fillin = FALSE)$chordal, NA)
  
  
  ##forming the final list:
  chordalComp <- list()
  if (length(conn.comp.imp)!=0){
    for (i in 1:length(conn.comp.imp)){
      chordalComp[[i]] <- list(component = conn.comp.imp[[i]],chordal = chordal[i])
    }
  }
 chordalComp
}

#########
## this function takes as parameter the adjacency matrix of a pdag (or cpdag)
## and returns the pattern of this pdag in the Meek sense, that is,
## it returns the adjacency matrix of the graph with the same skeleton where the only oriented
## edges are the v-structures (can be easily modified to work for MAGs/PAGs)
getPattern <- function(amat){
  
  ## makes the whole graph undirected
  tmp <- amat + t(amat)
  tmp[tmp == 2] <- 1
  
  ## find all v-structures i -> k <- j s.t. i not adj to k
  ## and make only those edges directed
  for (i in 1: (length(tmp[1,])-1))
    for (j in (i+1): length(tmp[1,])){
     if ((amat[j,i] ==0) & (amat[i,j] ==0) & (i!=j)){ ## if i no adjacent with j in G
        
        possible.k <- which(amat[,i]!= 0 & amat[i,]==0) ## finds all k such that i -> k is in G
        
        if (length(possible.k)!=0){    ## if there are any such k's then check whether j -> k for any of them
          for (k in 1: length(possible.k)){
            if ((amat[possible.k[k],j] ==1) & (amat[j,possible.k[k]]==0)) { ## if j -> k add the v-struc orientation to tmp
              tmp[i,possible.k[k]] <- 0   
              tmp[j,possible.k[k]] <- 0
            }
          }
        }
      }
     }
  tmp
}


####
## given the adjacency matrix of a pdag (or cpdag) - amat
## the following function obtains the corresponding cpdag of this pdag
## that is the cpdag that contains the same skeleton and v-structures as the given pdag
correspondingCpdag <- function(amat){
  patt <- getPattern(amat)
  corresp.cpdag <- addBgKnowledge(patt, checkInput = FALSE)
  corresp.cpdag
}

####
  
## amat - adjacency matrix of a dag/pdag/cpdag
## type - a string that is either dag, pdag or cpdag
## the following function returns TRUE if amat is inded of the specified type
##otherwise returns FALSE
isValidGraph <- function(amat, type = c("pdag", "cpdag", "dag"), verbose = FALSE) {
  
    ## we will first check that the amat contains only 0 and 1 as entries  
   tmp.amat <- amat
   tmp.amat[amat ==1] <- 0
   if (sum(tmp.amat) !=0){
       if (verbose) message("This is not a valid adjacency matrix.\n Your matrix must be of the type amat.cpdag\n Check: ?amatType for more info. \n\n")
     return(FALSE)
   }
  
   if (type == "dag"){ ## for a dag, we check for no cycles and no undirected edges
    okDir <- all(amat + t(amat) <= 1)
    if (!okDir) {
        if (verbose) message("There are undirected edges in this DAG. \n")
    }
    okCycle <- noCycles(amat)  
    if (!okCycle) {
        if (verbose) message("There is a cycle in this DAG. \n")
    }
    ok <- (okDir & okCycle)
    
  } else {
    
    if (type == "cpdag"){ 
      ## for a cpdag, we check that there are no cycles or partially directed cycles
      ## we check that the undirected part is made up of chordal components
      ## and we check that the cpdag is maximally oriented 
      
      ## check for cycles: TRUE if no cycles, FALSE otherwise
      ok1 <- noCycles(amat)
      
      ## check for chordality, ok2 is TRUE if chordality is satisfied, otherwise FALSE
      chordalComps <- chordalComponents(amat)
      chordal <- unlist(sapply(chordalComps,function(x) x[2]))
      ok2 <- !(FALSE %in% chordal)
      
      ## check that we cannot orient more edges by applying Meek rules
      ## ok3 is TRUE if the cpdag is maximally oriented, otherwise FALSE
      tmp <- t(applyOrientationRules(t(amat),verbose))
      isMaximal <- (sum(tmp-amat)==0)
      ok3 <- isMaximal
      
      ## the corresponding cpdag of this graph must have the same skeleton and orientations!
      corresp.cpdag <- correspondingCpdag(amat)
      diff.mat <- corresp.cpdag-amat
      ok4 <- (sum(abs(diff.mat))==0)
      
      ## we return the following value:
      ok <- (ok1 & ok2) & ok3 & ok4  
      
      if (!ok1){
          if (verbose) message("There is a directed or partially directed cycle in this CPDAG! \n\n")}
      
      if (!ok2){
          if (verbose) message("The undirected component of this CPDAG is not made up of chordal components! \nComponents: ",
                               which(chordal == "FALSE"),
                               " are not chordal. \n")
      }
        ## You can check this by calling: chordalComponents(Input)\n Be sure to replace Input with your adjacency matrix.\n\n")}
       
        if (!ok3){
          if (verbose) message("Your input CPDAG is not maximally oriented! \n You can obtain a maximally oriented CPDAG by calling: addBgKnowledge(amat)\n Be sure to replace amat with your adjacency matrix.\n\n")
      }
        if (!ok4){
            if (verbose) message("Your input CPDAG has more orientations then it should and/or contradicting orientations! Maybe you meant to specify a pdag? \n")}
    } else { 
      if (type=="pdag"){
        ## for a pdag we check whether all orientations in the corresponding cpdag are also in our pdag 
        ## and wehther the undirected part of the corresponding cpdag consists of chordal graphs
        ## we also check that the directed part of the pdag does not contain cycles
        ## and that the pdag is maximally oriented
        
        ## we first find the corresponding cpdag of our pdag
        ## and check that the undirected part of it is made up of chordal graphs
        ## ok1 is TRUE if this is the case and FALSE otherwise
        corresp.cpdag <- correspondingCpdag(amat)
        
        chordalComps <- chordalComponents(corresp.cpdag)
        chordal <- unlist(sapply(chordalComps,function(x) x[2]))
        
        ok1 <- !(FALSE %in% chordal)
        
        ## then we check that all orientations in the corresponding cpdag are also in our pdag
        ## ok4 is TRUE if there are no contradicting orientations between our pdag and the corresp. cpdag
        ## and FALSE otherwise
        dif.mat <- corresp.cpdag - amat ## should contain only positive values as 
                                        ## corresp.cpdag has the same skeleton as amat
                                        ## and should contain less orientations then amat
        dif.mat[dif.mat>0] <- 0
      
        ok4 <- (sum(dif.mat)>=0) ## MK: replaced '>' by '>='
        
        ## we then check whether the directed component of our pdag contains cycles
        ## ok2 is TRUE if there are no cycles and false otherwise
        amat.undir <- amat * t(amat) ## amat.undir - all undirected edges
        amat.dir <- amat - amat.undir      ## amat.dir - all directed edges
        
        ok2 <- noCycles(amat.dir)
        
        ## lastly we check whether by applying the Meek rules we can get more orientations
        ## ok3 is TRUE if our pdag is maximally oriented and FALSE otherwise
        tmp <- t(applyOrientationRules(t(amat),verbose))
        isMaximal <- (sum(tmp-amat)==0)
        
        ok3 <- isMaximal
        
        ## final result
        ok <- ok1 & ok2 & ok3 & ok4
        
        if (!ok1){
          if (verbose) message("The corresponding CPDAG of this PDAG is invalid. The undirected component of the CPDAG is not made up of chordal components!\n")
                  ## You can check this by calling: correspondingCpdag(chordalComponents(amat))\n Be sure to replace amat with your adjacency matrix.\n\n")
        }
        if (!ok4){
          if (verbose) message("Your input PDAG and its corresponding CPDAG have a mismatch in orientations!\n 
                  Check: graph::plot(as(t(amat),\"graphNEL\") \n
                  Be sure to replace amat with your adjacency matrix. \n\n")
        }
        
        if (!ok2){
          if (verbose) message("There is a directed cycle in this PDAG! \n\n")
        }
        
       if (!ok3){
          if (verbose) message("Your input PDAG is not maximally oriented! You can obtain a maximally oriented PDAG by calling: addBgKnowledge(amat)\n
                  Be sure to replace amat with your adjacency matrix.\n\n")
        }
        
      } else {
        if (verbose) message("Your input graph is not a DAG, PDAG, or CPDAG and so cannot be evaluated.\n")
        return(FALSE)
      }
     }
  }
  ok 
}


## PART 1: Unifying functions for possible De/An ####
possDe <- function(m, x, y = NULL, possible = TRUE, ds = TRUE, ## EXTERNAL
                   type = c("cpdag", "pdag", "dag", "mag", "pag") ) {
  ## only via def status paths
  if (ds) {
    if (type %in% c("dag", "cpdag", "pdag")) {
      res <- bPossibleDeProper(m=m, x=x, y=y, possible=possible)
      ## allows possible TRUE/FALSE
    } else { ## MAG, PAG
      ## !! only amat.pag coding => ok for this else-branch
      ## !! only single x
      if (length(x) == 1) {
        if (possible) {
          res <- possibleDe(amat=m, x=x)
        } else { ## possible = FALSE
          stop("Not yet implemented") 
        }
      } else { ## x is a set
        stop("Not yet implemented") 
      }
    }
  } else { ## via all paths
    if (length(x) == 1) {
      ## covers all types of graphs
      res <- possibleDeProper(m=m, x=x, y=y, possible=possible)
    } else { ## x is a set
      stop("Not yet implemented")
    } 
  }
  res
}

possAn <- function(m, x, y = NULL, possible = TRUE, ds = TRUE, ## EXTERNAL
                   type = c("cpdag", "pdag", "dag", "mag", "pag") ) {
  ## only PROPRER paths are considered
  ## m: adj matrix in coding according to type
  ## x: node positions
  ## y: paths must not go through y
  ## possible: find possible paths
  ## ds: find def status paths
  ## Value: 
  
  ## only via def status paths
  if (ds) {
    if (type %in% c("dag", "cpdag", "pdag")) {
      res <- bPossibleAnProper(m=m, x=x, y=y, possible=possible)
      ## allows possible TRUE/FALSE
    } else { ## MAG, PAG
      stop("Not yet implemented")
    }
  } else { ## via all paths
    if (length(x) == 1) {
      if (possible) {
        ## covers all types of graphs
        res <- possibleAnProper(m=m, x=x, y=y) 
      } else { ## possible = FALSE
        stop("Not yet implemented")
      }
    } else { ## x is a set
      stop("Not yet implemented")
    }
  }
  res
}

bPossibleDeProper <- function(m,x,y=NULL,possible = TRUE)
{
  #I will use depth first search
  #q denotes unvisited nodes/ nodes in queue
  #v denotes visited nodes
  q <- v <- previous <-  rep(0,length(m[,1])) 
  i <- k <-  1     
  if(length(x)>1){
    cat("Need to do this node by node!\n")
    return(NULL)
  }
  q <- sort(x)           
  tmp <- m
  ## previous will remember the previous node
  ## on the path, so we can check for definite status
  previous[1] <- q[1]
  
  while(q[k]!=0 & k<=i)
  {
    t <- q[k]
    #mark t as visited
    v[k] <- t       
    k <- k+1
    #in this for cycle adds all children of t and all nodes j  
    # such that t-j is in the pdag and  <previous[k-1],t,j> is of def. status
    ##i'm using the amat.cpdag encoding: amat[i,j] = 0, amat[j,i]=1 iff i -> j
    for(j in 1:length(tmp[1,])) 
      if (tmp[j,t] != 0  & tmp[t,j] != 2 ){#cat(previous[k-1],t,j,"\n")
        if ((tmp[j,t] ==1 & tmp[t,j] == 0) | (previous[k-1]==t) | (tmp[j,previous[k-1]] ==0 & tmp[previous[k-1],j] ==0))
          #only add nodes that haven't been added
          if (!(j %in% q)& !(j %in% y))   
          {
            i <- i+1
            previous[i] <- t
            q[i] <- j
          }}
  }
  ## remove all leftover zeros from initialization
  bpossDes <-setdiff(v,c(0))   
  
  return(sort(bpossDes))          
  
}

bPossibleAnProper <- function(m,x,y=NULL,possible = TRUE)
{
  #q denotes unvisited nodes/ nodes in queue
  #v denotes visited nodes
  q <- v <- previous <-  rep(0,length(m[,1])) 
  i <- k <-  1     
  if(length(x)>1){
    cat("Need to do this node by node!\n")
    return(NULL)
  }
  q <- sort(x)           
  tmp <- m
  previous[1] <- q[1]
  
  while(q[k]!=0 & k<=i)
  {
    t <- q[k]
    #mark t as visited
    v[k] <- t       
    k <- k+1
    #in this for cycle adds all parents of t and all nodes j  
    # such that t-j is in the pdag and  <previous[k-1],t,j> is of def. status
    ##i'm using the amat.cpdag encoding: amat[i,j] = 0, amat[j,i]=1 iff i -> j
    for(j in 1:length(tmp[1,])) 
      if (tmp[t,j] != 0){
        if ((tmp[j,t] ==0 & tmp[t,j] == 1) | (previous[k-1]==t) | ( tmp[j,previous[k-1]] ==0 & tmp[previous[k-1],j] ==0))
          #only add nodes that haven't been added
          # and that are b-possAn of t along a proper path w.r.t. to y! (i.e., not in y)
          if (!(j %in% q)& !(j %in% y))   
          {
            i <- i+1
            previous[i] <- t
            q[i] <- j
          }}
  }
  ## remove all leftover zeros from initialization
  bpossAnc <-setdiff(v,c(0))   
  
  return(sort(bpossAnc))          
  
}

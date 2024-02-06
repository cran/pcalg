## Adapting the dagitty functions to our problem

## need   pcalg2dagitty()   to utilize the functionality of dagitty

## as input this function takes the adjacency matrix of a graph, 
## labels=list of (observed) variable names and 
## type which can take values DAG, CPDAG, MAG, PAG
## the function outputs an object of class dagitty that corresponds
## to the same graph.

## the adjacency matrix encoding is amat.cpdag
pcalg2dagitty <- function(amat,labels,type="cpdag")
{
  if (type == "cpdag" | type=="pdag") {
      result <- "pdag{"
  }  else if (type=="dag") {
      result <- "dag{"  
  } else if (type=="mag") {
      result <- "mag{"
  } else if (type=="pag") {
      result <- "pag{"
  } 
  
  n <- length(amat[1,]) # =? ncol(amat) == #{nodes}
  
  ## if type of graph is DAG/CPDAG then edges are denoted with 0 and 1
  if (type %in% c("pdag","dag","cpdag")){
    
    for (i in 1:(n-1))
    {
      ## for single (unconnected) nodes
      if (sum(amat[,i]+amat[i,])==0){
        result <- paste0(result," ",labels[i]," ")
      }
      for(j in (i+1) : n)
      {
        if (amat[i,j] & amat[j,i]) {
          result <- paste(result,labels[i],"--",labels[j],sep=" ")
        } else if (amat[i,j] & !amat[j,i]) {
          result <- paste(result,labels[i],"<-",labels[j],sep=" ")
        } else if (!amat[i,j] & amat[j,i]) {
          result <- paste(result,labels[i],"->",labels[j],sep=" ")
        }
      }
    }
    
  } else if (type %in% c("mag","pag")) {
    
    ## if the type of graph is MAG/PAG then edgemarks are encoded with
    ## 1,2 and 3 in the adjacency matrix
      
      for (i in 1:(n-1))
      {
        ##for single nodes
        if (sum(amat[,i]+amat[i,])==0){
          result <- paste0(result," ",labels[i]," ")
        }
        for(j in (i+1) : n)
        {
          if (amat[i,j]==1 & amat[j,i]==1){
            result <- paste(result,labels[i],"@-@",labels[j],sep=" ")
          } else if (amat[i,j]==2 & amat[j,i]==3) {
            result <- paste(result,labels[i],"->",labels[j],sep=" ")
          } else if (amat[i,j]==3 & amat[j,i]==2) {
            result <- paste(result,labels[i],"<-",labels[j],sep=" ")
          } else if (amat[i,j]==3 & amat[j,i]==3) {
            result <- paste(result,labels[i],"--",labels[j],sep=" ")
          } else if (amat[i,j]==2 & amat[j,i]==2) {
            result <- paste(result,labels[i],"<->",labels[j],sep=" ")
          } else if (amat[i,j]==2 & amat[j,i]==1) {
            result <- paste(result,labels[i],"@->",labels[j],sep=" ")
          } else if (amat[i,j]==1 & amat[j,i]==2) {
            result <- paste(result,labels[i],"<-@",labels[j],sep=" ")
          }
        }
      }
  }
  ## if last node is unconnected
  if (sum(amat[,n]+amat[n,])==0){
    result <- paste0(result," ",labels[n]," ")
  }
  result <- paste0(result,"}")
  
  ## transform the resulting string into a daggity object
  dagitty::dagitty(result)
}


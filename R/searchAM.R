searchAM <- function(amat,x,type = c("an","de","ant","sp","nb","pa","ch","pde"))
{
    ## Purpose: in a DAG, CPDAG, MAG, or PAG determine the 
    ## ancestors/descendants/anteriors/spouses/neighbors/parents/children
    ## or possible descendants of a node x (or set of nodes x)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## - amat: adjacency matrix of type amat.pag
    ##           amat[i,j] = 0 iff no edge btw i,j
    ##           amat[i,j] = 1 iff i *-o j
    ##           amat[i,j] = 2 iff i *-> j
    ##           amat[i,j] = 3 iff i *-- j
    ## - x: node index or set of node indices of interest
    ## - type: "an" (ancestors, i.e., nodes y s.t. y-->...-->x), 
    ##         "de" (descendants, i.e., nodes y s.t. x-->...-->y), 
    ##         "ant" (anteriors, i.e., nodes y s.t. y---...---z-->...-->x,
    ##                i.e. there is an undirected path from y to a node z
    ##                followed by a directed path from z to x),
    ##         "sp" (spouses, i.e., nodes y s.t. y<->x),
    ##         "nb" (neighbors, i.e., nodes y s.t. y---x),
    ##         "pa" (parents, i.e., nodes y s.t. y-->x),
    ##         "ch" (children, i.e., nodes y s.t. y<--x), or
    ##         "pde" (possible descendants: nodes y s.t. there is a possibly
    ##                directed path from y to x: y {o,-}--{o,>} ... {o,-}--{o,>} x)
    ## ----------------------------------------------------------------------
    ## Value:
    ## - nodes: array containing the node indices related to x in the specified way
    ## ----------------------------------------------------------------------
    ## Author: Joris Mooij, Date: Oct 2020

    ## check inputs
    if( type != "de" && type != "pde" && type != "an" && type != "ant" && type != "sp" && type != "nb" && type != "pa" && type != "ch" ) 
      stop("Unknown type")
    stopifnot(is.matrix(amat))
    stopifnot(ncol(amat) == nrow(amat))

    ## initialize
    p <- nrow(amat)
    isType <- rep.int(FALSE, p)

    # start search at x
    indD <- c(x)
    if( type == "sp" || type == "nb" || type == "pa" || type == "ch" ) {
      while (length(indD) > 0) {
        d <- indD[1]
        indD <- indD[-1]
        if( type == "sp" ) {
          ## find all spouses of d
          indR <- which(amat[d,] == 2 & amat[,d] == 2) ## d <-> r
        } else if( type == "nb" ) {
          ## find all neighbors of d
          indR <- which(amat[d,] == 3 & amat[,d] == 3) ## d --- r
        } else if( type == "pa" ) {
          ## find all parents of d
          indR <- which(amat[d,] == 3 & amat[,d] == 2) ## d <-- r
        } else if( type == "ch" ) {
          ## find all children of d
          indR <- which(amat[d,] == 2 & amat[,d] == 3) ## d --> r
        }
        isType[indR] <- TRUE
      }
    } else {
      while (length(indD) > 0) {
        ## next element in the queue
        d <- indD[1]
        indD <- indD[-1]
        isType[d] <- TRUE
        if( type == "de" ) {
          ## find all children of d not visited yet
          indR <- which(amat[d,] == 2 & amat[,d] == 3 & !isType) ## d --> r
        } else if( type == "pde" ) {
          ## find all x s.t. d {o,-}--{o,>} x
          indR <- which((amat[d, ] == 1 | amat[d, ] == 2) &
                        (amat[, d] == 1 | amat[, d] == 3) & !isType)
        } else if( type == "an" ) {
          ## find all parents of d not visited yet
          indR <- which(amat[d,] == 3 & amat[,d] == 2 & !isType) ## d <-- ... <-- r
        } else if( type == "ant" ) {
          ## find all parents of d not visited yet
          indR <- which(amat[d,] == 3 & amat[,d] == 2 & !isType) ## d <-- ... <-- r
        }
        indD <- c(indD, indR)
      }
      if( type == "ant" ) {
        indD <- which(isType)
        while (length(indD) > 0) {
          ## next element in the queue
          d <- indD[1]
          indD <- indD[-1]
          isType[d] <- TRUE
          ## find all neighbors of d not visited yet
          indR <- which(amat[d,] == 3 & amat[,d] == 3 & !isType) ## d --- r
          indD <- c(indD, indR)
        }
      }
    }
    nodes<-which(isType)
    return(nodes)
} ## {searchAM}

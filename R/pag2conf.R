pag2conf <- function(P) {
  ## Purpose: Constructs a matrix which contains identifiably unconfounded
  ## node pairs in the Markov equivalence class represented
  ## by directed partial ancestral graph P.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - P: adjacency matrix of type amat.pag
  ##        P[i,j] = 0 iff no edge btw i,j
  ##        P[i,j] = 1 iff i *-o j
  ##        P[i,j] = 2 iff i *-> j
  ##        P[i,j] = 3 iff i *-- j
  ##      P has to be directed, i.e., it should not contain undirected
  ##      (x --- y) or circle-tail (x o-- y) edges
  ## ----------------------------------------------------------------------
  ## Value:
  ## - conf: matrix containing for each ordered pair of nodes whether the
  ##   absence of a bidiricted edge was identified
  ## ----------------------------------------------------------------------
  ## Author: Joris Mooij, Date: Nov 2020

  p<-dim(P)[1]
  conf<-matrix(0,p,p)
  rownames(conf)<-rownames(P)
  colnames(conf)<-colnames(P)
  for( i in seq_len(p) ) {
    for( j in seq_len(p) ) { 
      ## Proposition 6 in Mooij & Claassen (2020)
      if( P[i,j] == 0 ) { # nonadjacent ==> nonconfounded
        conf[i,j] = -1
        conf[j,i] = -1
      } else if( P[i,j] == 2 && P[j,i] == 3 ) { # i->j
        if( visibleEdge(P,i,j) ) { # not confounded
          conf[i,j] = -1
          conf[j,i] = -1
        }
      }
    }
  }

  return(conf)
} ## {pag2conf}
